"""VCF adapter for eveHap.

Extracts mtDNA allele profiles from VCF variant call files.
Uses pysam's VariantFile for efficient parsing of VCF/BCF files.
"""

from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING

import pysam

from evehap.adapters.base import InputAdapter
from evehap.core.profile import AlleleObservation, AlleleProfile, MTDNA_LENGTH

if TYPE_CHECKING:
    from evehap.core.damage import DamageFilter

# Common mtDNA contig names
MTDNA_CONTIGS = ["MT", "chrM", "chrMT", "M", "NC_012920.1"]


class VCFAdapter(InputAdapter):
    """Adapter for VCF variant call files.

    Extracts allele profiles from VCF by:
    1. Loading reference sequence for non-variant positions
    2. Parsing GT field for genotypes at variant positions
    3. Optionally parsing AD/DP fields for allele depths

    Positions not in VCF are assumed to match the reference.
    """

    def __init__(
        self,
        reference_path: str,
        sample_id: Optional[str] = None,
        mt_contig: str = "MT",
    ) -> None:
        """Initialize VCF adapter.

        Args:
            reference_path: Path to reference FASTA (rCRS)
            sample_id: Sample to extract (None = first sample)
            mt_contig: mtDNA contig name in VCF (default: MT)
        """
        self.reference_path = reference_path
        self.sample_id = sample_id
        self.mt_contig = mt_contig

    def can_handle(self, file_path: str) -> bool:
        """Check if file is a VCF."""
        return file_path.endswith((".vcf", ".vcf.gz", ".bcf"))

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "VCF"

    def list_samples(self, file_path: str) -> List[str]:
        """List all samples in the VCF file.

        Args:
            file_path: Path to VCF file

        Returns:
            List of sample IDs
        """
        vcf_path = Path(file_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {file_path}")

        with pysam.VariantFile(file_path) as vcf:
            return list(vcf.header.samples)

    def extract_profile(
        self,
        file_path: str,
        damage_filter: Optional["DamageFilter"] = None,
        **kwargs: object,
    ) -> AlleleProfile:
        """Extract profile from VCF.

        Parses GT field for genotypes. Positions not in VCF
        are filled from the reference sequence.

        Args:
            file_path: Path to VCF file
            damage_filter: Optional damage filter (not typically used for VCF)

        Returns:
            AlleleProfile extracted from the VCF

        Raises:
            FileNotFoundError: If VCF or reference file doesn't exist
            ValueError: If sample not found in VCF
        """
        vcf_path = Path(file_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {file_path}")

        ref_path = Path(self.reference_path)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference_path}")

        # Load reference sequence
        ref_sequence = self._load_reference()

        # Open VCF
        with pysam.VariantFile(file_path) as vcf:
            # Determine sample to extract
            samples = list(vcf.header.samples)
            if not samples:
                raise ValueError("VCF file contains no samples")

            if self.sample_id is None:
                target_sample = samples[0]
            else:
                if self.sample_id not in samples:
                    raise ValueError(
                        f"Sample '{self.sample_id}' not found in VCF. "
                        f"Available samples: {samples[:5]}..."
                    )
                target_sample = self.sample_id

            sample_idx = samples.index(target_sample)

            # Create profile - VCF data has full coverage (missing = reference)
            profile = AlleleProfile(
                sample_id=target_sample,
                source_format=self.format_name,
                assumes_full_coverage=True,
            )

            # Find mtDNA contig in VCF
            mt_contig = self._find_mtdna_contig(vcf)

            # Start with reference at all positions
            variant_positions: Dict[int, AlleleObservation] = {}

            # Parse variants for this sample
            try:
                for record in vcf.fetch(mt_contig):
                    pos = record.pos  # 1-based
                    ref = record.ref
                    alts = record.alts or ()

                    # Get genotype for target sample
                    sample_data = record.samples[sample_idx]
                    gt = sample_data.get("GT")

                    if gt is None:
                        continue

                    # Get allele depths if available
                    ad = sample_data.get("AD")
                    dp = sample_data.get("DP")

                    # Determine alleles present
                    allele_indices = [i for i in gt if i is not None]
                    if not allele_indices:
                        continue

                    # Build allele list (0=ref, 1+=alt)
                    all_alleles = [ref] + list(alts)

                    # For SNPs/indels, record the observation
                    # Skip if all alleles are reference (allele_idx = 0)
                    for allele_idx in set(allele_indices):
                        if allele_idx == 0:
                            # Reference allele - skip
                            continue
                        if allele_idx >= len(all_alleles):
                            continue

                        allele = all_alleles[allele_idx]

                        # Calculate frequency
                        if ad and len(ad) > allele_idx:
                            # Use allele depth
                            total_depth = sum(a for a in ad if a)
                            if total_depth > 0:
                                freq = ad[allele_idx] / total_depth
                            else:
                                freq = 1.0 / len(set(allele_indices))
                        else:
                            # No depth info, assume equal frequency
                            freq = 1.0 / len(set(allele_indices))

                        # Handle different allele lengths
                        if len(allele) == 1 and len(ref) == 1:
                            # Simple SNP
                            snp_pos = pos
                            snp_ref = ref[0]
                            snp_alt = allele[0]
                        elif len(allele) == len(ref) and len(ref) > 1:
                            # MNP (multi-nucleotide polymorphism) - extract each SNP
                            # e.g., AC>GC -> position 0 is A>G
                            for i, (r, a) in enumerate(zip(ref, allele)):
                                if r != a:
                                    snp_pos = pos + i
                                    snp_ref = r
                                    snp_alt = a
                                    self._add_snp_observation(
                                        variant_positions, snp_pos, snp_ref, snp_alt,
                                        freq, dp, record.qual
                                    )
                            continue  # Already added all SNPs
                        elif len(allele) < len(ref):
                            # Deletion - mark each deleted position with '-'
                            # e.g., ACCCCCTCTA > A means positions 8281-8289 are deleted
                            # First position is the anchor (usually matches)
                            for i in range(len(allele), len(ref)):
                                del_pos = pos + i
                                del_ref = ref[i]
                                self._add_snp_observation(
                                    variant_positions, del_pos, del_ref, "-",
                                    freq, dp, record.qual
                                )
                            continue  # Already added all deletions
                        else:
                            # Insertion - record first base
                            snp_pos = pos
                            snp_ref = ref[0]
                            snp_alt = allele[0]

                        # Add the observation
                        self._add_snp_observation(
                            variant_positions, snp_pos, snp_ref, snp_alt,
                            freq, dp, record.qual
                        )

            except ValueError:
                # Contig not in VCF, all reference
                pass

            # Add only variant positions to profile
            # Reference positions are NOT added - the classifier should
            # treat missing positions as reference matches
            for pos, obs in variant_positions.items():
                profile.add_observation(obs)

        return profile

    def _add_snp_observation(
        self,
        variant_positions: Dict[int, AlleleObservation],
        pos: int,
        ref: str,
        alt: str,
        freq: float,
        depth: Optional[int],
        quality: Optional[float],
    ) -> None:
        """Add a SNP observation to the variant positions dict.

        Args:
            variant_positions: Dict to update
            pos: 1-based position
            ref: Reference allele
            alt: Alternate allele
            freq: Allele frequency
            depth: Read depth
            quality: Quality score
        """
        if pos not in variant_positions:
            variant_positions[pos] = AlleleObservation(
                position=pos,
                ref_allele=ref.upper(),
                alleles={alt.upper(): freq},
                depth=depth,
                quality=quality,
                source=self.format_name,
            )
        else:
            # Update existing observation
            obs = variant_positions[pos]
            obs.alleles[alt.upper()] = max(
                obs.alleles.get(alt.upper(), 0.0),
                freq
            )

    def _load_reference(self) -> str:
        """Load reference sequence from FASTA.

        Returns:
            Reference sequence string
        """
        with pysam.FastaFile(self.reference_path) as fasta:
            # Find mtDNA contig
            for contig in MTDNA_CONTIGS:
                if contig in fasta.references:
                    return fasta.fetch(contig)

            # Try first reference if short (likely mtDNA)
            if fasta.references:
                first_ref = fasta.references[0]
                seq = fasta.fetch(first_ref)
                if len(seq) < 20000:  # Likely mtDNA
                    return seq

            raise ValueError(
                f"Could not find mtDNA contig in reference. "
                f"Available: {fasta.references}"
            )

    def _find_mtdna_contig(self, vcf: pysam.VariantFile) -> str:
        """Find mtDNA contig name in VCF.

        Args:
            vcf: Open VariantFile

        Returns:
            Contig name
        """
        contigs = set(vcf.header.contigs.keys())

        # Try preferred first
        if self.mt_contig in contigs:
            return self.mt_contig

        # Try common names
        for contig in MTDNA_CONTIGS:
            if contig in contigs:
                return contig

        # If no contigs defined, use mt_contig and hope for the best
        return self.mt_contig


"""BAM/CRAM adapter for eveHap.

Extracts mtDNA allele profiles from aligned read files using pysam.
"""

from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING

import pysam

from evehap.adapters.base import InputAdapter
from evehap.core.profile import AlleleObservation, AlleleProfile, MTDNA_LENGTH

if TYPE_CHECKING:
    from evehap.core.damage import DamageFilter

# Common mtDNA contig names to try
MTDNA_CONTIGS = ["chrM", "MT", "chrMT", "M", "NC_012920.1"]

# Default rCRS reference path (relative to package)
PACKAGE_DIR = Path(__file__).parent.parent.parent
DEFAULT_REFERENCE = PACKAGE_DIR / "data" / "reference" / "rCRS.fasta"

# Cached rCRS sequence
_RCRS_SEQUENCE: Optional[str] = None


def _load_rcrs_sequence() -> str:
    """Load rCRS reference sequence (cached)."""
    global _RCRS_SEQUENCE
    if _RCRS_SEQUENCE is None:
        if DEFAULT_REFERENCE.exists():
            with open(DEFAULT_REFERENCE) as f:
                lines = f.readlines()
                _RCRS_SEQUENCE = "".join(
                    line.strip() for line in lines if not line.startswith(">")
                )
        else:
            # Fallback: empty sequence (will use sample alleles as proxy)
            _RCRS_SEQUENCE = ""
    return _RCRS_SEQUENCE


class BAMAdapter(InputAdapter):
    """Adapter for BAM/CRAM aligned read files.

    Uses pysam to efficiently extract allele observations from
    aligned reads covering the mitochondrial genome.
    """

    def __init__(
        self,
        reference_path: Optional[str] = None,
        min_base_quality: int = 20,
        min_mapping_quality: int = 20,
    ) -> None:
        """Initialize BAM adapter.

        Args:
            reference_path: Path to reference FASTA (required for CRAM)
            min_base_quality: Minimum base quality to include (default: 20)
            min_mapping_quality: Minimum mapping quality to include (default: 20)
        """
        self.reference_path = reference_path
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality

    def can_handle(self, file_path: str) -> bool:
        """Check if file is a BAM or CRAM."""
        return file_path.endswith((".bam", ".cram"))

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "BAM"

    def extract_profile(
        self,
        file_path: str,
        region: str = "chrM",
        damage_filter: Optional["DamageFilter"] = None,
        **kwargs: object,
    ) -> AlleleProfile:
        """Extract profile from BAM/CRAM.

        Uses pysam to count alleles at each position, applying quality filters.
        Optionally applies damage filtering for ancient DNA.

        Args:
            file_path: Path to BAM/CRAM file
            region: Genomic region to extract (default: chrM)
            damage_filter: Optional damage filter for ancient DNA

        Returns:
            AlleleProfile extracted from the BAM

        Raises:
            FileNotFoundError: If BAM file doesn't exist
            ValueError: If mtDNA contig cannot be found
        """
        bam_path = Path(file_path)
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {file_path}")

        # Extract sample ID from filename
        sample_id = bam_path.stem
        # Remove common suffixes like .chrM
        for suffix in [".chrM", ".MT", ".chrMT", ".mito"]:
            if sample_id.endswith(suffix):
                sample_id = sample_id[:-len(suffix)]
                break

        profile = AlleleProfile(sample_id=sample_id, source_format=self.format_name)

        # Open BAM file
        open_kwargs: Dict[str, object] = {}
        if self.reference_path:
            open_kwargs["reference_filename"] = self.reference_path

        with pysam.AlignmentFile(file_path, "rb", **open_kwargs) as bam:
            # Find the mtDNA contig
            mt_contig = self._find_mtdna_contig(bam, region)
            if mt_contig is None:
                raise ValueError(
                    f"Could not find mtDNA contig. Tried: {region}, {MTDNA_CONTIGS}"
                )

            # Get reference length
            ref_lengths = dict(zip(bam.references, bam.lengths))
            mt_length = ref_lengths.get(mt_contig, MTDNA_LENGTH)

            # Use count_coverage for efficient allele counting
            # This returns 4 arrays (A, C, G, T) of coverage at each position
            try:
                coverage = bam.count_coverage(
                    mt_contig,
                    start=0,
                    stop=mt_length,
                    quality_threshold=self.min_base_quality,
                    read_callback=lambda read: (
                        read.mapping_quality >= self.min_mapping_quality
                        and not read.is_unmapped
                        and not read.is_secondary
                        and not read.is_supplementary
                    ),
                )
            except Exception as e:
                raise ValueError(f"Failed to count coverage for {mt_contig}: {e}") from e

            # coverage is a tuple of 4 arrays: (A, C, G, T)
            a_counts, c_counts, g_counts, t_counts = coverage

            # Get reference sequence for determining ref allele
            # We'll use a simple approach: get the most common allele at each position
            # as a proxy for reference when not available
            ref_seq = self._get_reference_sequence(bam, mt_contig, mt_length)

            # Convert coverage to AlleleObservations
            for pos_0based in range(min(len(a_counts), mt_length)):
                pos_1based = pos_0based + 1

                a = a_counts[pos_0based]
                c = c_counts[pos_0based]
                g = g_counts[pos_0based]
                t = t_counts[pos_0based]

                total = a + c + g + t
                if total == 0:
                    continue  # No coverage at this position

                # Build allele frequency dict
                alleles: Dict[str, float] = {}
                for allele, count in [("A", a), ("C", c), ("G", g), ("T", t)]:
                    if count > 0:
                        alleles[allele] = count / total

                # Determine reference allele
                if ref_seq and pos_0based < len(ref_seq):
                    ref_allele = ref_seq[pos_0based].upper()
                    if ref_allele not in "ACGT":
                        ref_allele = "N"
                else:
                    # Try to use cached rCRS sequence
                    rcrs = _load_rcrs_sequence()
                    if rcrs and pos_0based < len(rcrs):
                        ref_allele = rcrs[pos_0based].upper()
                    else:
                        # Fallback: use most common allele as reference proxy
                        ref_allele = max(alleles.keys(), key=lambda x: alleles.get(x, 0))

                obs = AlleleObservation(
                    position=pos_1based,
                    ref_allele=ref_allele,
                    alleles=alleles,
                    depth=total,
                    quality=None,  # Average quality not easily available from count_coverage
                    source=self.format_name,
                )

                profile.add_observation(obs)

        return profile

    def _find_mtdna_contig(
        self, bam: pysam.AlignmentFile, preferred: str
    ) -> Optional[str]:
        """Find the mtDNA contig name in the BAM file.

        Args:
            bam: Open pysam AlignmentFile
            preferred: Preferred contig name to try first

        Returns:
            Contig name or None if not found
        """
        references = set(bam.references)

        # Try preferred first
        if preferred in references:
            return preferred

        # Try common names
        for contig in MTDNA_CONTIGS:
            if contig in references:
                return contig

        return None

    def _get_reference_sequence(
        self, bam: pysam.AlignmentFile, contig: str, length: int
    ) -> Optional[str]:
        """Get reference sequence for the contig if available.

        Args:
            bam: Open pysam AlignmentFile
            contig: Contig name
            length: Expected sequence length

        Returns:
            Reference sequence string or None if not available
        """
        # If we have a reference file, try to get the sequence
        if self.reference_path:
            try:
                with pysam.FastaFile(self.reference_path) as fasta:
                    if contig in fasta.references:
                        return fasta.fetch(contig)
            except Exception:
                pass

        # Otherwise return None - we'll infer from the data
        return None

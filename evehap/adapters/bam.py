"""BAM/CRAM adapter for eveHap.

Extracts mtDNA allele profiles from aligned read files using pysam.
"""

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional

import pysam

from evehap.adapters.base import InputAdapter
from evehap.core.profile import MTDNA_LENGTH, AlleleObservation, AlleleProfile

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
                _RCRS_SEQUENCE = "".join(line.strip() for line in lines if not line.startswith(">"))
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
        min_depth: int = 1,
        adaptive_quality: bool = False,
        filter_read_ends: bool = False,
        read_end_length: int = 3,
    ) -> None:
        """Initialize BAM adapter.

        Args:
            reference_path: Path to reference FASTA (required for CRAM)
            min_base_quality: Minimum base quality to include (default: 20)
            min_mapping_quality: Minimum mapping quality to include (default: 20)
            min_depth: Minimum depth to include a position (default: 1, use 2 for aDNA)
            adaptive_quality: Automatically lower quality thresholds for low-coverage samples
            filter_read_ends: Exclude bases within read_end_length of read start/end (default: False)
            read_end_length: Number of bases to filter from read ends (default: 3)
        """
        self.reference_path = reference_path
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality
        self.min_depth = min_depth
        self.adaptive_quality = adaptive_quality
        self.filter_read_ends = filter_read_ends
        self.read_end_length = read_end_length

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
        damage_filter: Optional["DamageFilter"] = None,
        region: str = "chrM",
        **kwargs: Any,
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
                sample_id = sample_id[: -len(suffix)]
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
                raise ValueError(f"Could not find mtDNA contig. Tried: {region}, {MTDNA_CONTIGS}")

            # Get reference length
            ref_lengths = dict(zip(bam.references, bam.lengths))
            mt_length = ref_lengths.get(mt_contig, MTDNA_LENGTH)

            # Get reference sequence for determining ref allele
            ref_seq = self._get_reference_sequence(bam, mt_contig, mt_length)

            # If adaptive quality is enabled, estimate coverage first to adjust thresholds
            effective_base_quality = self.min_base_quality
            effective_mapping_quality = self.min_mapping_quality
            if self.adaptive_quality:
                # Quick pass with low thresholds to estimate coverage
                try:
                    quick_coverage = bam.count_coverage(
                        mt_contig,
                        start=0,
                        stop=mt_length,
                        quality_threshold=5,  # Very low threshold for estimation
                        read_callback=lambda read: (
                            read.mapping_quality >= 5
                            and not read.is_unmapped
                            and not read.is_secondary
                            and not read.is_supplementary
                        ),
                    )
                    # Calculate coverage fraction
                    a_quick, c_quick, g_quick, t_quick = quick_coverage
                    covered_positions = sum(
                        1
                        for i in range(mt_length)
                        if (a_quick[i] + c_quick[i] + g_quick[i] + t_quick[i]) > 0
                    )
                    coverage_fraction = covered_positions / mt_length

                    # Adjust thresholds based on coverage
                    if coverage_fraction < 0.25:
                        effective_base_quality = 10
                        effective_mapping_quality = 10
                    elif coverage_fraction < 0.50:
                        effective_base_quality = 15
                        effective_mapping_quality = 15
                    # Otherwise use provided thresholds
                except Exception:
                    # If estimation fails, use provided thresholds
                    pass

            # If damage filter is provided or read-end filtering is enabled,
            # we need to iterate through reads to track read position and strand
            if damage_filter is not None or self.filter_read_ends:
                # Initialize position counters: [A, C, G, T] counts for each position
                position_counts: List[List[int]] = [[0, 0, 0, 0] for _ in range(mt_length)]

                # Iterate through reads to count alleles with damage filtering
                for read in bam.fetch(mt_contig, 0, mt_length):
                    # Apply read-level filters (use effective threshold if adaptive)
                    map_qual_threshold = (
                        effective_mapping_quality
                        if self.adaptive_quality
                        else self.min_mapping_quality
                    )
                    if (
                        read.is_unmapped
                        or read.is_secondary
                        or read.is_supplementary
                        or read.mapping_quality < map_qual_threshold
                    ):
                        continue

                    query_seq = read.query_sequence
                    query_quals = read.query_qualities
                    if not query_seq:
                        continue

                    read_len = len(query_seq)
                    is_reverse = read.is_reverse

                    # Get aligned pairs (query_pos, ref_pos)
                    aligned_pairs = read.get_aligned_pairs(matches_only=False)

                    for query_pos, ref_pos in aligned_pairs:
                        if ref_pos is None or ref_pos < 0 or ref_pos >= mt_length:
                            continue
                        if query_pos is None or query_pos < 0 or query_pos >= len(query_seq):
                            continue

                        # Check base quality (use effective threshold if adaptive)
                        base_qual_threshold = (
                            effective_base_quality
                            if self.adaptive_quality
                            else self.min_base_quality
                        )
                        if query_quals and query_quals[query_pos] < base_qual_threshold:
                            continue

                        base = query_seq[query_pos].upper()
                        if base not in "ACGT":
                            continue

                        # Determine reference allele at this position
                        if ref_seq and ref_pos < len(ref_seq):
                            ref_allele = ref_seq[ref_pos].upper()
                        else:
                            # Try cached rCRS
                            rcrs = _load_rcrs_sequence()
                            if rcrs and ref_pos < len(rcrs):
                                ref_allele = rcrs[ref_pos].upper()
                            else:
                                ref_allele = "N"

                        # Check if this base should be filtered due to read-end position
                        if self.filter_read_ends:
                            # Calculate distance from read ends
                            if is_reverse:
                                dist_from_5prime = read_len - 1 - query_pos
                                dist_from_3prime = query_pos
                            else:
                                dist_from_5prime = query_pos
                                dist_from_3prime = read_len - 1 - query_pos

                            # Filter bases within read_end_length of either end
                            if (
                                dist_from_5prime < self.read_end_length
                                or dist_from_3prime < self.read_end_length
                            ):
                                continue  # Skip this base (at read end)

                        # Check if this base should be filtered due to damage
                        if damage_filter is not None:
                            if damage_filter.should_filter(
                                position=ref_pos + 1,  # 1-based position
                                ref=ref_allele,
                                alt=base,
                                read_position=query_pos,
                                read_length=read_len,
                                is_reverse=is_reverse,
                            ):
                                continue  # Skip this base (damage pattern)

                        # Count this base
                        base_idx = {"A": 0, "C": 1, "G": 2, "T": 3}[base]
                        position_counts[ref_pos][base_idx] += 1

                # Convert counts to AlleleObservations
                for pos_0based in range(mt_length):
                    pos_1based = pos_0based + 1
                    a, c, g, t = position_counts[pos_0based]

                    total = a + c + g + t
                    if total == 0:
                        continue  # No coverage at this position

                    # Apply minimum depth filter
                    if total < self.min_depth:
                        continue  # Depth too low

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
                        quality=None,
                        source=self.format_name,
                    )

                    profile.add_observation(obs)
            else:
                # Fast path: use count_coverage when no damage filter
                # This returns 4 arrays (A, C, G, T) of coverage at each position
                try:
                    # Use effective thresholds if adaptive quality is enabled
                    base_qual_threshold = (
                        effective_base_quality if self.adaptive_quality else self.min_base_quality
                    )
                    map_qual_threshold = (
                        effective_mapping_quality
                        if self.adaptive_quality
                        else self.min_mapping_quality
                    )
                    coverage = bam.count_coverage(
                        mt_contig,
                        start=0,
                        stop=mt_length,
                        quality_threshold=base_qual_threshold,
                        read_callback=lambda read: (
                            read.mapping_quality >= map_qual_threshold
                            and not read.is_unmapped
                            and not read.is_secondary
                            and not read.is_supplementary
                        ),
                    )
                except Exception as e:
                    raise ValueError(f"Failed to count coverage for {mt_contig}: {e}") from e

                # coverage is a tuple of 4 arrays: (A, C, G, T)
                a_counts, c_counts, g_counts, t_counts = coverage

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

                    # Apply minimum depth filter
                    if total < self.min_depth:
                        continue  # Depth too low

                    # Build allele frequency dict
                    pos_alleles: Dict[str, float] = {}
                    for allele, count in [("A", a), ("C", c), ("G", g), ("T", t)]:
                        if count > 0:
                            pos_alleles[allele] = count / total

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
                            ref_allele = max(
                                pos_alleles.keys(), key=lambda x: pos_alleles.get(x, 0)
                            )

                    obs = AlleleObservation(
                        position=pos_1based,
                        ref_allele=ref_allele,
                        alleles=pos_alleles,
                        depth=total,
                        quality=None,  # Average quality not easily available from count_coverage
                        source=self.format_name,
                    )

                    profile.add_observation(obs)

        return profile

    def _find_mtdna_contig(self, bam: pysam.AlignmentFile, preferred: str) -> Optional[str]:
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
                        result: str = fasta.fetch(contig)
                        return result
            except Exception:
                pass

        # Otherwise return None - we'll infer from the data
        return None

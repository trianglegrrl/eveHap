"""AlleleProfile and AlleleObservation data structures for eveHap.

These are the core data structures that represent mtDNA allele observations
from any input source (BAM, VCF, FASTA, etc.).
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set, Tuple

# mtDNA genome length (rCRS)
MTDNA_LENGTH = 16569

# Heteroplasmy threshold (minor allele frequency to consider heteroplasmic)
HETEROPLASMY_THRESHOLD = 0.10


@dataclass
class AlleleObservation:
    """Allele observation at a single mtDNA position.

    Represents the observed allele(s) at a position in the mtDNA genome,
    along with associated quality metrics.

    Attributes:
        position: 1-based position in mtDNA (1-16569)
        ref_allele: rCRS reference allele at this position
        alleles: Dictionary mapping allele to frequency (e.g., {'G': 0.98, 'A': 0.02})
        depth: Read depth at position (None for FASTA/microarray)
        quality: Mean base quality (None if unavailable)
        source: Source format identifier ('bam', 'vcf', 'fasta', etc.)
    """

    position: int
    ref_allele: str
    alleles: Dict[str, float]
    depth: Optional[int] = None
    quality: Optional[float] = None
    source: str = "unknown"

    def __post_init__(self) -> None:
        """Validate the observation."""
        if not 1 <= self.position <= MTDNA_LENGTH:
            raise ValueError(
                f"Position {self.position} out of range [1, {MTDNA_LENGTH}]"
            )
        if not self.alleles:
            raise ValueError("At least one allele must be provided")
        if not self.ref_allele:
            raise ValueError("Reference allele must be provided")

    @property
    def major_allele(self) -> str:
        """Return the most frequent allele.

        In case of a tie, returns the alphabetically first allele.
        """
        if not self.alleles:
            raise ValueError("No alleles present")

        max_freq = max(self.alleles.values())
        # Get all alleles with max frequency (handle ties)
        candidates = [a for a, f in self.alleles.items() if f == max_freq]
        # Return alphabetically first to ensure consistency
        return sorted(candidates)[0]

    @property
    def is_variant(self) -> bool:
        """True if major allele differs from reference allele."""
        return self.major_allele != self.ref_allele

    @property
    def is_heteroplasmic(self) -> bool:
        """True if multiple alleles present above heteroplasmy threshold.

        Default threshold is 10% minor allele frequency.
        """
        if len(self.alleles) <= 1:
            return False

        # Count alleles above threshold
        above_threshold = sum(
            1 for freq in self.alleles.values()
            if freq >= HETEROPLASMY_THRESHOLD
        )
        return above_threshold > 1

    def allele_frequency(self, allele: str) -> float:
        """Get frequency of specific allele, 0.0 if absent.

        Args:
            allele: The allele to query (A, C, G, T, or other)

        Returns:
            Frequency of the allele (0.0-1.0), or 0.0 if not present
        """
        return self.alleles.get(allele, 0.0)


class AlleleProfile:
    """Complete mtDNA allele profile for a sample.

    Stores allele observations for all covered positions and provides
    methods for querying the profile.

    Attributes:
        sample_id: Sample identifier
        source_format: Source format ('bam', 'vcf', 'fasta', etc.)
        observations: Dictionary mapping position to AlleleObservation
        covered_positions: Set of all covered positions
        metadata: Additional sample metadata
    """

    def __init__(
        self,
        sample_id: str,
        source_format: str,
        assumes_full_coverage: bool = False,
    ) -> None:
        """Initialize an AlleleProfile.

        Args:
            sample_id: Sample identifier
            source_format: Source format identifier
            assumes_full_coverage: If True, missing positions are assumed to be
                covered with the reference allele. This is typical for VCF data
                where only variants are recorded.
        """
        self.sample_id = sample_id
        self.source_format = source_format
        self.assumes_full_coverage = assumes_full_coverage
        self.observations: Dict[int, AlleleObservation] = {}
        self.covered_positions: Set[int] = set()
        self.metadata: Dict[str, Any] = {}

    def add_observation(self, obs: AlleleObservation) -> None:
        """Add an allele observation at a position.

        If an observation already exists at this position, it is replaced.

        Args:
            obs: The AlleleObservation to add
        """
        self.observations[obs.position] = obs
        self.covered_positions.add(obs.position)

    def get_allele(self, position: int) -> Optional[str]:
        """Get major allele at position, None if not covered.

        Args:
            position: 1-based mtDNA position

        Returns:
            Major allele at position, or None if position is not covered
        """
        obs = self.observations.get(position)
        if obs is None:
            return None
        return obs.major_allele

    def get_observation(self, position: int) -> Optional[AlleleObservation]:
        """Get full observation at position.

        Args:
            position: 1-based mtDNA position

        Returns:
            AlleleObservation at position, or None if not covered
        """
        return self.observations.get(position)

    def is_derived(self, position: int, derived_allele: str) -> Optional[bool]:
        """Check if position has derived allele.

        Args:
            position: 1-based mtDNA position
            derived_allele: The allele to check for

        Returns:
            True if position has derived allele, False if not,
            None if position is not covered
        """
        allele = self.get_allele(position)
        if allele is None:
            return None
        return allele == derived_allele

    def coverage_fraction(self) -> float:
        """Fraction of mtDNA genome covered (0.0-1.0).

        If assumes_full_coverage is True (e.g., VCF data), returns 1.0.

        Returns:
            Fraction of the 16569 bp mtDNA genome that has coverage
        """
        if self.assumes_full_coverage:
            return 1.0
        return len(self.covered_positions) / MTDNA_LENGTH

    def mean_depth(self) -> Optional[float]:
        """Mean read depth across covered positions.

        Only considers observations that have depth information.

        Returns:
            Mean depth, or None if no observations have depth
        """
        depths = [
            obs.depth for obs in self.observations.values()
            if obs.depth is not None
        ]
        if not depths:
            return None
        return sum(depths) / len(depths)

    def get_variants(self) -> List[Tuple[int, str, str]]:
        """Return list of (position, ref, alt) variants vs reference.

        Returns:
            List of tuples (position, ref_allele, alt_allele) for all
            positions where the major allele differs from reference
        """
        variants = []
        for pos, obs in self.observations.items():
            if obs.is_variant:
                variants.append((pos, obs.ref_allele, obs.major_allele))
        return variants

    def __repr__(self) -> str:
        """Return string representation."""
        return (
            f"AlleleProfile(sample_id={self.sample_id!r}, "
            f"source={self.source_format!r}, "
            f"positions={len(self.covered_positions)}, "
            f"coverage={self.coverage_fraction():.1%})"
        )


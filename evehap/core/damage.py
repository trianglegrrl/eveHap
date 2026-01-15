"""Ancient DNA damage filtering for eveHap.

Ancient DNA exhibits characteristic damage patterns:
- C→T substitutions at 5' ends (deamination of cytosine to uracil)
- G→A substitutions at 3' ends (same damage on opposite strand)

This module provides tools to:
1. Detect damage patterns in sequencing data
2. Filter potentially damaged bases from profiles
3. Estimate damage rates for QC
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional

import pysam

if TYPE_CHECKING:
    from evehap.core.profile import AlleleProfile


@dataclass
class DamageStats:
    """Statistics about damage patterns in a sample.

    Attributes:
        ct_5prime_rate: C→T rate at 5' positions
        ga_3prime_rate: G→A rate at 3' positions
        total_reads: Total reads analyzed
        damaged_reads: Reads with damage signature
        position_rates: Position-specific damage rates
    """

    ct_5prime_rate: float = 0.0
    ga_3prime_rate: float = 0.0
    total_reads: int = 0
    damaged_reads: int = 0
    position_rates: Dict[int, float] = field(default_factory=dict)

    @property
    def is_ancient(self) -> bool:
        """True if damage patterns suggest ancient DNA.

        Typically, ancient DNA has >3% C→T at 5' ends.
        """
        return self.ct_5prime_rate > 0.03 or self.ga_3prime_rate > 0.03

    @property
    def damage_rate(self) -> float:
        """Overall damage rate (average of 5' and 3' rates)."""
        return (self.ct_5prime_rate + self.ga_3prime_rate) / 2


class DamageFilter:
    """Filter for ancient DNA damage patterns.

    Ancient DNA exhibits characteristic damage:
    - C→T substitutions at 5' ends (deamination of cytosine)
    - G→A substitutions at 3' ends (deamination on opposite strand)

    This filter can:
    1. Identify potentially damaged bases
    2. Remove or downweight damaged observations
    3. Estimate damage rates from BAM files
    """

    def __init__(
        self,
        filter_ct_5prime: bool = True,
        filter_ga_3prime: bool = True,
        end_positions: int = 3,
        min_damage_rate: float = 0.1,
    ) -> None:
        """Initialize damage filter.

        Args:
            filter_ct_5prime: Filter C→T at 5' ends
            filter_ga_3prime: Filter G→A at 3' ends
            end_positions: Number of positions from read ends to filter
            min_damage_rate: Minimum damage rate to apply filtering
        """
        self.filter_ct_5prime = filter_ct_5prime
        self.filter_ga_3prime = filter_ga_3prime
        self.end_positions = end_positions
        self.min_damage_rate = min_damage_rate

    def should_filter(
        self,
        position: int,
        ref: str,
        alt: str,
        read_position: int,
        read_length: int,
        is_reverse: bool,
    ) -> bool:
        """Determine if this observation should be filtered.

        Args:
            position: mtDNA position (not used, for context)
            ref: Reference allele
            alt: Observed allele
            read_position: Position within read (0-based)
            read_length: Total read length
            is_reverse: True if read is reverse strand

        Returns:
            True if observation matches damage pattern and should be filtered
        """
        # Calculate distance from ends
        if is_reverse:
            # For reverse reads, coordinates are flipped
            dist_from_5prime = read_length - 1 - read_position
            dist_from_3prime = read_position
        else:
            dist_from_5prime = read_position
            dist_from_3prime = read_length - 1 - read_position

        # Check C→T at 5' end
        if self.filter_ct_5prime:
            if ref == "C" and alt == "T":
                if dist_from_5prime < self.end_positions:
                    return True

        # Check G→A at 3' end
        if self.filter_ga_3prime:
            if ref == "G" and alt == "A":
                if dist_from_3prime < self.end_positions:
                    return True

        return False

    def filter_profile(
        self,
        profile: "AlleleProfile",
        remove_potential_damage: bool = False,
        damage_stats: Optional[DamageStats] = None,
        reference: Optional[str] = None,
    ) -> "AlleleProfile":
        """Filter a profile to remove or mark potentially damaged observations.

        For ancient DNA, this filters:
        - C→T variants (most common damage pattern)
        - G→A variants (complementary damage on opposite strand)

        When reference is provided, checks ref→alt transitions.
        Otherwise, filters all T variants (conservative approach for aDNA).

        Args:
            profile: AlleleProfile to filter
            remove_potential_damage: If True, remove damaged observations
            damage_stats: Optional damage stats for context
            reference: Optional reference sequence for precise filtering

        Returns:
            Filtered AlleleProfile (may be same object if no filtering)
        """
        from evehap.core.profile import AlleleProfile

        if not remove_potential_damage:
            return profile

        # Create filtered profile
        filtered = AlleleProfile(
            sample_id=profile.sample_id,
            source_format=profile.source_format,
        )
        filtered.metadata = profile.metadata.copy()
        filtered.metadata["damage_filtered"] = True

        damage_removed = 0
        total_variants = 0

        for pos, obs in profile.observations.items():
            # Check if marked as potential damage
            if hasattr(obs, "metadata") and obs.metadata:
                if obs.metadata.get("potential_damage", False):
                    damage_removed += 1
                    continue  # Skip this observation

            # Check for C→T / G→A damage pattern
            if obs.is_variant and obs.major_allele:
                total_variants += 1
                is_damage = False

                if reference and 1 <= pos <= len(reference):
                    ref_base = reference[pos - 1].upper()
                    alt_base = obs.major_allele.upper()
                    # C→T damage
                    if ref_base == "C" and alt_base == "T":
                        is_damage = True
                    # G→A damage (complement)
                    elif ref_base == "G" and alt_base == "A":
                        is_damage = True
                else:
                    # Without reference, filter all T transitions
                    # (conservative for aDNA)
                    if obs.major_allele.upper() == "T":
                        is_damage = True

                if is_damage:
                    damage_removed += 1
                    continue  # Skip damaged observation

            filtered.add_observation(obs)

        filtered.metadata["damage_removed"] = damage_removed
        filtered.metadata["original_variants"] = total_variants

        return filtered

    def filter_profile_strict(
        self,
        profile: "AlleleProfile",
        reference: str,
    ) -> "AlleleProfile":
        """Strict damage filtering: remove all C→T and G→A variants.

        This is aggressive but effective for heavily damaged aDNA.

        Args:
            profile: AlleleProfile to filter
            reference: Reference sequence (required)

        Returns:
            Filtered AlleleProfile with damage removed
        """
        return self.filter_profile(
            profile,
            remove_potential_damage=True,
            reference=reference,
        )

    def estimate_damage(
        self,
        bam_path: str,
        max_reads: int = 10000,
        mt_contig: str = "chrM",
    ) -> DamageStats:
        """Estimate damage rates from BAM file.

        Samples reads to estimate C→T and G→A rates at read ends.

        Args:
            bam_path: Path to BAM file
            max_reads: Maximum number of reads to sample
            mt_contig: mtDNA contig name

        Returns:
            DamageStats with estimated damage rates
        """
        stats = DamageStats()

        bam_file = Path(bam_path)
        if not bam_file.exists():
            return stats

        # Count damage events
        ct_5prime = 0  # C→T at 5' end
        c_5prime = 0  # Total C at 5' end
        ga_3prime = 0  # G→A at 3' end
        g_3prime = 0  # Total G at 3' end

        reads_analyzed = 0
        damaged_reads = 0

        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                # Find mtDNA contig
                contigs = list(bam.references)
                mt = None
                for c in [mt_contig, "MT", "chrMT", "M"]:
                    if c in contigs:
                        mt = c
                        break

                if mt is None:
                    return stats

                # Sample reads
                for read in bam.fetch(mt):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    if reads_analyzed >= max_reads:
                        break

                    reads_analyzed += 1
                    read_damaged = False

                    # Get aligned pairs (query_pos, ref_pos, ref_base)
                    try:
                        aligned_pairs = read.get_aligned_pairs(with_seq=True)
                    except Exception:
                        continue

                    query_seq = read.query_sequence
                    if not query_seq:
                        continue

                    read_len = len(query_seq)
                    is_reverse = read.is_reverse

                    for query_pos, _ref_pos, ref_base in aligned_pairs:
                        if query_pos is None or ref_base is None:
                            continue

                        ref_base = ref_base.upper()
                        query_base = query_seq[query_pos].upper()

                        if ref_base not in "ACGT" or query_base not in "ACGT":
                            continue

                        # Calculate distance from ends
                        if is_reverse:
                            dist_from_5prime = read_len - 1 - query_pos
                            dist_from_3prime = query_pos
                        else:
                            dist_from_5prime = query_pos
                            dist_from_3prime = read_len - 1 - query_pos

                        # Count 5' C and C→T
                        if dist_from_5prime < self.end_positions:
                            if ref_base == "C":
                                c_5prime += 1
                                if query_base == "T":
                                    ct_5prime += 1
                                    read_damaged = True

                        # Count 3' G and G→A
                        if dist_from_3prime < self.end_positions:
                            if ref_base == "G":
                                g_3prime += 1
                                if query_base == "A":
                                    ga_3prime += 1
                                    read_damaged = True

                    if read_damaged:
                        damaged_reads += 1

        except Exception:
            pass

        # Calculate rates
        stats.total_reads = reads_analyzed
        stats.damaged_reads = damaged_reads
        stats.ct_5prime_rate = ct_5prime / c_5prime if c_5prime > 0 else 0.0
        stats.ga_3prime_rate = ga_3prime / g_3prime if g_3prime > 0 else 0.0

        return stats

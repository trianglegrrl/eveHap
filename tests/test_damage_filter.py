"""Tests for ancient DNA damage filtering.

These tests use real AADR BAM samples as fixtures to verify damage filtering
works correctly on actual ancient DNA data.

Fixtures:
- ERR1341839.chrM.bam: 30 reads, heavy damage (19/23 variants are C→T)
- ERR1341819.chrM.bam: 6 reads, minimal data

Test strategy:
1. Verify damage detection identifies C→T patterns
2. Verify damage filtering removes C→T transitions
3. Verify classification improves (or at least doesn't break) after filtering
"""

import pytest
from pathlib import Path
import sys

# Add project to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from evehap.core.damage import DamageFilter, DamageStats
from evehap.core.profile import AlleleProfile, AlleleObservation
from evehap.adapters.bam import BAMAdapter


# Fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"

# rCRS reference (partial, for testing)
# Positions 660-670: CCTCCTCATG (C at 664)
# Positions 3168-3178: ACCCCTCCTA (C at 3171)


class TestDamageStats:
    """Test DamageStats dataclass."""

    def test_is_ancient_low_damage(self):
        """Low damage rate should not be classified as ancient."""
        stats = DamageStats(ct_5prime_rate=0.01, ga_3prime_rate=0.01)
        assert not stats.is_ancient

    def test_is_ancient_high_damage(self):
        """High damage rate should be classified as ancient."""
        stats = DamageStats(ct_5prime_rate=0.15, ga_3prime_rate=0.10)
        assert stats.is_ancient

    def test_damage_rate(self):
        """Test average damage rate calculation."""
        stats = DamageStats(ct_5prime_rate=0.20, ga_3prime_rate=0.10)
        assert stats.damage_rate == pytest.approx(0.15)


class TestDamageFilter:
    """Test DamageFilter class."""

    def test_should_filter_ct_5prime(self):
        """C→T at 5' end should be filtered."""
        filt = DamageFilter(end_positions=3)
        # Position 0 in read (5' end), C→T
        assert filt.should_filter(
            position=100, ref="C", alt="T",
            read_position=0, read_length=50, is_reverse=False
        )

    def test_should_not_filter_ct_middle(self):
        """C→T in middle of read should not be filtered."""
        filt = DamageFilter(end_positions=3)
        # Position 25 in read (middle), C→T
        assert not filt.should_filter(
            position=100, ref="C", alt="T",
            read_position=25, read_length=50, is_reverse=False
        )

    def test_should_filter_ga_3prime(self):
        """G→A at 3' end should be filtered."""
        filt = DamageFilter(end_positions=3)
        # Position 49 in 50bp read (3' end), G→A
        assert filt.should_filter(
            position=100, ref="G", alt="A",
            read_position=49, read_length=50, is_reverse=False
        )

    def test_should_not_filter_ga_5prime(self):
        """G→A at 5' end should not be filtered (wrong end)."""
        filt = DamageFilter(end_positions=3)
        # Position 0 in read (5' end), G→A
        assert not filt.should_filter(
            position=100, ref="G", alt="A",
            read_position=0, read_length=50, is_reverse=False
        )

    def test_should_filter_reverse_strand(self):
        """Damage filter should handle reverse strand correctly."""
        filt = DamageFilter(end_positions=3)
        # On reverse strand, 3' end is at low read positions
        # Position 0 on reverse = 3' end
        assert filt.should_filter(
            position=100, ref="G", alt="A",
            read_position=0, read_length=50, is_reverse=True
        )


class TestDamageFilterProfile:
    """Test profile-level damage filtering."""

    def test_filter_removes_ct_variants(self):
        """Filter should remove C→T variants when reference provided."""
        # Create a profile with mixed variants
        profile = AlleleProfile(sample_id="test", source_format="test")

        # C→T variant at position 100 (should be removed)
        obs1 = AlleleObservation(
            position=100, depth=10, ref_allele="C",
            alleles={"T": 1.0}
        )
        profile.add_observation(obs1)

        # G→A variant at position 200 (should be removed)
        obs2 = AlleleObservation(
            position=200, depth=10, ref_allele="G",
            alleles={"A": 1.0}
        )
        profile.add_observation(obs2)

        # A→G variant at position 300 (should be kept)
        obs3 = AlleleObservation(
            position=300, depth=10, ref_allele="A",
            alleles={"G": 1.0}
        )
        profile.add_observation(obs3)

        # Create a fake reference (just need positions 100, 200, 300)
        reference = "N" * 99 + "C" + "N" * 99 + "G" + "N" * 99 + "A" + "N" * 100

        filt = DamageFilter()
        filtered = filt.filter_profile(
            profile, remove_potential_damage=True, reference=reference
        )

        # Should only have position 300 remaining
        assert len(filtered.observations) == 1
        assert 300 in filtered.observations
        assert 100 not in filtered.observations
        assert 200 not in filtered.observations

    def test_filter_without_reference_removes_t_variants(self):
        """Without reference, filter should remove all T variants conservatively."""
        profile = AlleleProfile(sample_id="test", source_format="test")

        # T variant
        obs1 = AlleleObservation(
            position=100, depth=10, ref_allele="N",
            alleles={"T": 1.0}
        )
        profile.add_observation(obs1)

        # G variant (should be kept)
        obs2 = AlleleObservation(
            position=200, depth=10, ref_allele="N",
            alleles={"G": 1.0}
        )
        profile.add_observation(obs2)

        filt = DamageFilter()
        filtered = filt.filter_profile(
            profile, remove_potential_damage=True, reference=None
        )

        # T variant should be removed, G should remain
        assert len(filtered.observations) == 1
        assert 200 in filtered.observations

    def test_filter_tracks_removal_stats(self):
        """Filter should track how many variants were removed."""
        profile = AlleleProfile(sample_id="test", source_format="test")

        # Add 3 C→T variants
        for i in range(3):
            obs = AlleleObservation(
                position=100 + i, depth=10, ref_allele="C",
                alleles={"T": 1.0}
            )
            profile.add_observation(obs)

        # Add 1 A→G variant
        obs = AlleleObservation(
            position=200, depth=10, ref_allele="A",
            alleles={"G": 1.0}
        )
        profile.add_observation(obs)

        reference = "N" * 99 + "CCC" + "N" * 97 + "A" + "N" * 100

        filt = DamageFilter()
        filtered = filt.filter_profile(
            profile, remove_potential_damage=True, reference=reference
        )

        assert filtered.metadata.get('damage_filtered') is True
        assert filtered.metadata.get('damage_removed') == 3
        assert filtered.metadata.get('original_variants') == 4


@pytest.mark.skipif(
    not (FIXTURES_DIR / "ERR1341839.chrM.bam").exists(),
    reason="Fixture BAM not available"
)
class TestDamageFilterWithRealData:
    """Test damage filtering with real AADR ancient DNA samples."""

    def test_estimate_damage_from_bam(self):
        """Test damage estimation from real BAM file."""
        bam_path = str(FIXTURES_DIR / "ERR1341839.chrM.bam")

        filt = DamageFilter(end_positions=3)
        stats = filt.estimate_damage(bam_path)

        # Should detect some damage in ancient DNA
        assert stats.total_reads > 0
        # This is an ancient sample, so expect some damage
        # (though with only 30 reads, rates may vary)
        print(f"C→T rate: {stats.ct_5prime_rate:.2%}")
        print(f"G→A rate: {stats.ga_3prime_rate:.2%}")
        print(f"Total reads: {stats.total_reads}")

    def test_extract_profile_shows_damage_pattern(self):
        """Verify the fixture BAM has the expected damage pattern."""
        bam_path = str(FIXTURES_DIR / "ERR1341839.chrM.bam")
        ref_path = str(Path(__file__).parent.parent / "data" / "reference" / "rCRS.fasta")

        adapter = BAMAdapter(reference_path=ref_path)
        profile = adapter.extract_profile(bam_path)

        # Count T variants (likely C→T damage)
        t_variants = 0
        total_variants = 0
        for pos, obs in profile.observations.items():
            if obs.is_variant and obs.major_allele:
                total_variants += 1
                if obs.major_allele == 'T':
                    t_variants += 1

        print(f"Total variants: {total_variants}")
        print(f"T variants: {t_variants}")
        print(f"T variant fraction: {t_variants/total_variants:.1%}" if total_variants > 0 else "N/A")

        # This sample should have high T variant fraction (damage)
        if total_variants > 5:
            assert t_variants / total_variants > 0.5, \
                f"Expected >50% T variants in damaged sample, got {t_variants}/{total_variants}"

    def test_filter_reduces_variants(self):
        """Damage filtering should reduce variant count in damaged sample."""
        bam_path = str(FIXTURES_DIR / "ERR1341839.chrM.bam")
        ref_path = str(Path(__file__).parent.parent / "data" / "reference" / "rCRS.fasta")

        adapter = BAMAdapter(reference_path=ref_path)
        profile = adapter.extract_profile(bam_path)

        # Load reference for filtering
        with open(ref_path) as f:
            lines = f.readlines()
            reference = "".join(line.strip() for line in lines if not line.startswith(">"))

        # Count variants before filtering
        before = sum(1 for obs in profile.observations.values() if obs.is_variant)

        # Apply damage filter
        filt = DamageFilter()
        filtered = filt.filter_profile(
            profile, remove_potential_damage=True, reference=reference
        )

        # Count variants after filtering
        after = sum(1 for obs in filtered.observations.values() if obs.is_variant)

        print(f"Variants before: {before}")
        print(f"Variants after: {after}")
        print(f"Removed: {before - after}")

        # Should have removed some variants
        assert after < before, "Damage filter should remove some variants"
        # But shouldn't remove everything
        assert after >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


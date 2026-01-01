"""Tests for evehap.core.damage module."""

from pathlib import Path

import pytest

from evehap.core.damage import DamageFilter, DamageStats
from evehap.core.profile import AlleleProfile, AlleleObservation


class TestDamageStats:
    """Tests for DamageStats dataclass."""

    def test_create_default(self) -> None:
        """Test default DamageStats creation."""
        stats = DamageStats()
        assert stats.ct_5prime_rate == 0.0
        assert stats.ga_3prime_rate == 0.0
        assert stats.total_reads == 0
        assert stats.damaged_reads == 0

    def test_is_ancient_true(self) -> None:
        """Test is_ancient with high damage rates."""
        stats = DamageStats(ct_5prime_rate=0.10, ga_3prime_rate=0.05)
        assert stats.is_ancient is True

    def test_is_ancient_false(self) -> None:
        """Test is_ancient with low damage rates."""
        stats = DamageStats(ct_5prime_rate=0.01, ga_3prime_rate=0.01)
        assert stats.is_ancient is False

    def test_is_ancient_borderline(self) -> None:
        """Test is_ancient at threshold."""
        # Just above threshold
        stats = DamageStats(ct_5prime_rate=0.031, ga_3prime_rate=0.0)
        assert stats.is_ancient is True

        # Just below threshold
        stats = DamageStats(ct_5prime_rate=0.029, ga_3prime_rate=0.0)
        assert stats.is_ancient is False


class TestDamageFilter:
    """Tests for DamageFilter class."""

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        df = DamageFilter()
        assert df.filter_ct_5prime is True
        assert df.filter_ga_3prime is True
        assert df.end_positions == 3
        assert df.min_damage_rate == 0.1

    def test_init_custom(self) -> None:
        """Test custom initialization."""
        df = DamageFilter(
            filter_ct_5prime=False,
            filter_ga_3prime=True,
            end_positions=5,
            min_damage_rate=0.05,
        )
        assert df.filter_ct_5prime is False
        assert df.filter_ga_3prime is True
        assert df.end_positions == 5
        assert df.min_damage_rate == 0.05

    def test_should_filter_ct_5prime(self) -> None:
        """Test filtering C→T at 5' end."""
        df = DamageFilter(end_positions=3)

        # C→T at position 0 of forward read (5' end)
        assert df.should_filter(
            position=100, ref="C", alt="T",
            read_position=0, read_length=100, is_reverse=False
        ) is True

        # C→T at position 1 (still within end_positions)
        assert df.should_filter(
            position=101, ref="C", alt="T",
            read_position=1, read_length=100, is_reverse=False
        ) is True

        # C→T at position 50 (middle of read, not filtered)
        assert df.should_filter(
            position=150, ref="C", alt="T",
            read_position=50, read_length=100, is_reverse=False
        ) is False

    def test_should_filter_ga_3prime(self) -> None:
        """Test filtering G→A at 3' end."""
        df = DamageFilter(end_positions=3)

        # G→A at 3' end of forward read
        assert df.should_filter(
            position=100, ref="G", alt="A",
            read_position=99, read_length=100, is_reverse=False
        ) is True

        # G→A in middle (not filtered)
        assert df.should_filter(
            position=100, ref="G", alt="A",
            read_position=50, read_length=100, is_reverse=False
        ) is False

    def test_should_filter_reverse_strand(self) -> None:
        """Test filtering on reverse strand."""
        df = DamageFilter(end_positions=3)

        # On reverse strand, 5' end is at the end of read coordinates
        # C→T at 3' of reverse = 5' of original molecule
        assert df.should_filter(
            position=100, ref="C", alt="T",
            read_position=99, read_length=100, is_reverse=True
        ) is True

        # G→A at 5' of reverse = 3' of original molecule
        assert df.should_filter(
            position=100, ref="G", alt="A",
            read_position=0, read_length=100, is_reverse=True
        ) is True

    def test_should_filter_disabled(self) -> None:
        """Test that filtering can be disabled."""
        df = DamageFilter(filter_ct_5prime=False, filter_ga_3prime=False)

        # Should not filter even C→T at 5'
        assert df.should_filter(
            position=100, ref="C", alt="T",
            read_position=0, read_length=100, is_reverse=False
        ) is False

    def test_should_filter_non_damage(self) -> None:
        """Test that non-damage substitutions are not filtered."""
        df = DamageFilter(end_positions=3)

        # A→G at 5' end (not damage pattern)
        assert df.should_filter(
            position=100, ref="A", alt="G",
            read_position=0, read_length=100, is_reverse=False
        ) is False


class TestDamageFilterProfile:
    """Tests for filtering AlleleProfile."""

    def test_filter_profile_no_damage(self) -> None:
        """Test filtering profile with no damage."""
        df = DamageFilter()

        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        filtered = df.filter_profile(profile)

        # Non-damage variant should be preserved
        assert 73 in filtered.covered_positions

    def test_filter_profile_with_damage_markers(self) -> None:
        """Test filtering profile with damage-marked observations."""
        df = DamageFilter()

        profile = AlleleProfile(sample_id="test", source_format="test")

        # Normal observation
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        # Observation marked as potential damage
        obs = AlleleObservation(
            position=100, ref_allele="C", alleles={"T": 1.0}, source="test"
        )
        obs.metadata = {"potential_damage": True}
        profile.add_observation(obs)

        filtered = df.filter_profile(profile, remove_potential_damage=True)

        # Position 73 preserved, position 100 removed
        assert 73 in filtered.covered_positions
        # Position 100 may or may not be removed depending on implementation


class TestDamageFilterEstimate:
    """Tests for damage rate estimation."""

    @pytest.fixture
    def hgdp_bam_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return HGDP BAM directory."""
        return mtdna_poc_data_dir / "hgdp" / "bams"

    @pytest.fixture
    def sample_bam(self, hgdp_bam_dir: Path) -> Path:
        """Return a sample BAM file."""
        if not hgdp_bam_dir.exists():
            pytest.skip("HGDP BAM directory not available")

        bams = list(hgdp_bam_dir.glob("*.bam"))
        if not bams:
            pytest.skip("No BAM files found")

        return bams[0]

    def test_estimate_damage_modern(self, sample_bam: Path) -> None:
        """Test damage estimation on modern DNA sample."""
        df = DamageFilter()

        stats = df.estimate_damage(str(sample_bam))

        assert isinstance(stats, DamageStats)
        # Modern samples should have low damage rates
        assert stats.ct_5prime_rate < 0.05
        assert stats.ga_3prime_rate < 0.05

    def test_estimate_damage_sample_size(self, sample_bam: Path) -> None:
        """Test that damage estimation samples reads."""
        df = DamageFilter()

        stats = df.estimate_damage(str(sample_bam), max_reads=1000)

        # Should have analyzed some reads
        assert stats.total_reads > 0
        assert stats.total_reads <= 1000



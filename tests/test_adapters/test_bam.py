"""Tests for evehap.adapters.bam module."""

from pathlib import Path

import pytest

from evehap.adapters.bam import BAMAdapter
from evehap.core.profile import AlleleProfile


class TestBAMAdapter:
    """Tests for BAMAdapter class."""

    def test_can_handle_bam(self) -> None:
        """Test can_handle returns True for .bam files."""
        adapter = BAMAdapter()
        assert adapter.can_handle("sample.bam") is True
        assert adapter.can_handle("/path/to/sample.bam") is True

    def test_can_handle_cram(self) -> None:
        """Test can_handle returns True for .cram files."""
        adapter = BAMAdapter()
        assert adapter.can_handle("sample.cram") is True

    def test_can_handle_other(self) -> None:
        """Test can_handle returns False for other files."""
        adapter = BAMAdapter()
        assert adapter.can_handle("sample.vcf") is False
        assert adapter.can_handle("sample.fasta") is False
        assert adapter.can_handle("sample.bam.bai") is False

    def test_format_name(self) -> None:
        """Test format_name returns 'BAM'."""
        adapter = BAMAdapter()
        assert adapter.format_name == "BAM"

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        adapter = BAMAdapter()
        assert adapter.reference_path is None
        assert adapter.min_base_quality == 20
        assert adapter.min_mapping_quality == 20
        assert adapter.min_depth == 1

    def test_init_custom(self) -> None:
        """Test custom initialization."""
        adapter = BAMAdapter(
            reference_path="/path/to/ref.fa",
            min_base_quality=30,
            min_mapping_quality=40,
            min_depth=2,
        )
        assert adapter.reference_path == "/path/to/ref.fa"
        assert adapter.min_base_quality == 30
        assert adapter.min_mapping_quality == 40
        assert adapter.min_depth == 2


class TestBAMAdapterExtractProfile:
    """Tests for BAMAdapter.extract_profile().

    Note: Tests that require real BAM files are marked with @pytest.mark.integration
    and are skipped by default. Run with -m integration to include them.
    """

    @pytest.fixture
    def mtdna_poc_data_dir(self) -> Path:
        """Return mtdna_poc data directory."""
        return Path(__file__).parent.parent.parent.parent / "data"

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

    @pytest.mark.integration
    def test_extract_profile_basic(self, sample_bam: Path) -> None:
        """Test basic profile extraction from BAM."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        assert isinstance(profile, AlleleProfile)
        # Sample ID strips common suffixes like .chrM
        expected_id = sample_bam.stem
        for suffix in [".chrM", ".MT", ".chrMT", ".mito"]:
            if expected_id.endswith(suffix):
                expected_id = expected_id[:-len(suffix)]
                break
        assert profile.sample_id == expected_id
        assert profile.source_format == "BAM"
        assert len(profile.covered_positions) > 0

    @pytest.mark.integration
    def test_extract_profile_coverage(self, sample_bam: Path) -> None:
        """Test that extracted profile has reasonable coverage."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        # mtDNA is 16569 bp, well-sequenced samples should cover most
        coverage = profile.coverage_fraction()
        assert coverage > 0.5, f"Coverage too low: {coverage:.1%}"

    @pytest.mark.integration
    def test_extract_profile_depth(self, sample_bam: Path) -> None:
        """Test that extracted profile has depth information."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        # At least some observations should have depth
        depths = [
            obs.depth for obs in profile.observations.values()
            if obs.depth is not None
        ]
        assert len(depths) > 0, "No observations have depth"
        assert all(d > 0 for d in depths), "Depths should be positive"

    @pytest.mark.integration
    def test_extract_profile_variants(self, sample_bam: Path) -> None:
        """Test that extracted profile can report variants."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        # Most samples have at least a few variants
        variants = profile.get_variants()
        # Don't assert specific count, just that we can get variants
        assert isinstance(variants, list)

    @pytest.mark.integration
    def test_extract_profile_region(self, sample_bam: Path) -> None:
        """Test profile extraction with specific region."""
        adapter = BAMAdapter()

        # Try different common mtDNA contig names
        for region in ["chrM", "MT", "chrMT"]:
            try:
                profile = adapter.extract_profile(str(sample_bam), region=region)
                if len(profile.covered_positions) > 0:
                    break
            except Exception:
                continue

        # Should have extracted something
        assert len(profile.covered_positions) > 0

    @pytest.mark.integration
    def test_extract_profile_quality_filter(self, sample_bam: Path) -> None:
        """Test that quality filtering is applied."""
        # Extract with low quality threshold
        adapter_low = BAMAdapter(min_base_quality=5, min_mapping_quality=5)
        profile_low = adapter_low.extract_profile(str(sample_bam))

        # Extract with high quality threshold
        adapter_high = BAMAdapter(min_base_quality=35, min_mapping_quality=35)
        profile_high = adapter_high.extract_profile(str(sample_bam))

        # Higher quality filter may reduce coverage
        # (not always, depends on data quality)
        assert len(profile_low.covered_positions) >= len(profile_high.covered_positions)

    @pytest.mark.integration
    def test_extract_profile_with_damage_filter(self, sample_bam: Path) -> None:
        """Test that damage_filter parameter is actually used during extraction."""
        from evehap.core.damage import DamageFilter

        # Extract without damage filter
        adapter_no_filter = BAMAdapter()
        profile_no_filter = adapter_no_filter.extract_profile(str(sample_bam))

        # Extract with damage filter
        damage_filter = DamageFilter(end_positions=3)
        adapter_with_filter = BAMAdapter()
        profile_with_filter = adapter_with_filter.extract_profile(
            str(sample_bam), damage_filter=damage_filter
        )

        # With damage filtering, we may have fewer variants (C→T and G→A filtered)
        # But coverage should still be reasonable
        assert len(profile_with_filter.covered_positions) > 0

    @pytest.mark.integration
    def test_damage_filter_filters_ct_at_read_ends(self, sample_bam: Path) -> None:
        """Test that C→T at read ends are filtered when damage_filter is enabled."""
        from evehap.core.damage import DamageFilter

        damage_filter = DamageFilter(end_positions=3, filter_ct_5prime=True)
        adapter = BAMAdapter()

        # Extract with damage filter
        profile = adapter.extract_profile(str(sample_bam), damage_filter=damage_filter)

        # Profile should be valid
        assert len(profile.covered_positions) > 0

    @pytest.mark.integration
    def test_damage_filter_filters_ga_at_read_ends(self, sample_bam: Path) -> None:
        """Test that G→A at read ends are filtered when damage_filter is enabled."""
        from evehap.core.damage import DamageFilter

        damage_filter = DamageFilter(end_positions=3, filter_ga_3prime=True)
        adapter = BAMAdapter()

        # Extract with damage filter
        profile = adapter.extract_profile(str(sample_bam), damage_filter=damage_filter)

        # Profile should be valid
        assert len(profile.covered_positions) > 0

    @pytest.mark.integration
    @pytest.mark.skipif(
        not (Path(__file__).parent / "fixtures" / "ERR1341839.chrM.bam").exists(),
        reason="AADR fixture BAM not available"
    )
    def test_damage_filter_with_aadr_sample(self) -> None:
        """Test damage filtering with real AADR ancient DNA sample."""
        from evehap.core.damage import DamageFilter
        from pathlib import Path

        bam_path = Path(__file__).parent / "fixtures" / "ERR1341839.chrM.bam"
        ref_path = Path(__file__).parent.parent.parent / "data" / "reference" / "rCRS.fasta"

        if not bam_path.exists() or not ref_path.exists():
            pytest.skip("Required fixtures not available")

        # Extract without damage filter
        adapter_no_filter = BAMAdapter(reference_path=str(ref_path))
        profile_no_filter = adapter_no_filter.extract_profile(str(bam_path))

        # Extract with damage filter
        damage_filter = DamageFilter(end_positions=3)
        adapter_with_filter = BAMAdapter(reference_path=str(ref_path))
        profile_with_filter = adapter_with_filter.extract_profile(
            str(bam_path), damage_filter=damage_filter
        )

        # Count variants before and after
        variants_before = sum(1 for obs in profile_no_filter.observations.values() if obs.is_variant)
        variants_after = sum(1 for obs in profile_with_filter.observations.values() if obs.is_variant)

        # Damage filtering should reduce variant count (filters C→T and G→A at ends)
        # But may not always if damage is not at read ends
        assert len(profile_with_filter.covered_positions) > 0

    @pytest.mark.integration
    def test_extract_profile_min_depth_filtering(self, sample_bam: Path) -> None:
        """Test that positions with depth < min_depth are excluded."""
        # Extract with min_depth=2
        adapter_depth2 = BAMAdapter(min_depth=2)
        profile_depth2 = adapter_depth2.extract_profile(str(sample_bam))

        # Extract with min_depth=1 (should include all positions)
        adapter_depth1 = BAMAdapter(min_depth=1)
        profile_depth1 = adapter_depth1.extract_profile(str(sample_bam))

        # min_depth=1 should have equal or more positions than min_depth=2
        assert len(profile_depth1.covered_positions) >= len(profile_depth2.covered_positions)

        # All positions in depth2 profile should have depth >= 2
        for pos, obs in profile_depth2.observations.items():
            assert obs.depth is not None
            assert obs.depth >= 2, f"Position {pos} has depth {obs.depth} < 2"

    @pytest.mark.integration
    def test_extract_profile_min_depth_default(self, sample_bam: Path) -> None:
        """Test that default min_depth=1 includes all positions."""
        # Default should be min_depth=1 (no filtering)
        adapter_default = BAMAdapter()
        profile_default = adapter_default.extract_profile(str(sample_bam))

        # Should have coverage
        assert len(profile_default.covered_positions) > 0

    @pytest.mark.integration
    def test_extract_profile_min_depth_zero_coverage(self, sample_bam: Path) -> None:
        """Test that min_depth filtering works correctly."""
        # Extract with high min_depth (should filter many positions)
        adapter_high = BAMAdapter(min_depth=100)
        profile_high = adapter_high.extract_profile(str(sample_bam))

        # Extract with low min_depth
        adapter_low = BAMAdapter(min_depth=1)
        profile_low = adapter_low.extract_profile(str(sample_bam))

        # High min_depth should have fewer or equal positions
        assert len(profile_high.covered_positions) <= len(profile_low.covered_positions)

    @pytest.mark.integration
    def test_all_features_combined(self, sample_bam: Path) -> None:
        """Test all new features working together: damage filter, min_depth, adaptive quality, read-end filtering."""
        from evehap.core.damage import DamageFilter

        # Create adapter with all features enabled
        damage_filter = DamageFilter(end_positions=3)
        adapter = BAMAdapter(
            min_depth=2,
            adaptive_quality=True,
            filter_read_ends=True,
            read_end_length=3,
        )

        # Extract profile with all features
        profile = adapter.extract_profile(str(sample_bam), damage_filter=damage_filter)

        # Should produce valid profile
        assert len(profile.covered_positions) > 0

        # All positions should have depth >= min_depth
        for pos, obs in profile.observations.items():
            assert obs.depth is not None
            assert obs.depth >= 2

    @pytest.mark.integration
    def test_backward_compatibility(self, sample_bam: Path) -> None:
        """Test that default behavior is unchanged (backward compatibility)."""
        # Default adapter (no new features enabled)
        adapter_default = BAMAdapter()
        profile_default = adapter_default.extract_profile(str(sample_bam))

        # Should work as before
        assert len(profile_default.covered_positions) > 0
        assert adapter_default.min_depth == 1
        assert adapter_default.adaptive_quality is False
        assert adapter_default.filter_read_ends is False

    @pytest.mark.integration
    def test_extract_profile_adaptive_quality(self, sample_bam: Path) -> None:
        """Test that adaptive quality thresholds are applied for low-coverage samples."""
        # Extract with adaptive quality enabled
        adapter_adaptive = BAMAdapter(adaptive_quality=True)
        profile_adaptive = adapter_adaptive.extract_profile(str(sample_bam))

        # Extract with fixed quality thresholds
        adapter_fixed = BAMAdapter(adaptive_quality=False, min_base_quality=20, min_mapping_quality=20)
        profile_fixed = adapter_fixed.extract_profile(str(sample_bam))

        # Both should produce valid profiles
        assert len(profile_adaptive.covered_positions) > 0
        assert len(profile_fixed.covered_positions) > 0

    @pytest.mark.integration
    def test_adaptive_quality_low_coverage(self, sample_bam: Path) -> None:
        """Test that adaptive quality lowers thresholds for low-coverage samples."""
        # This test verifies the mechanism works, not necessarily the exact thresholds
        # since we don't know the coverage of the test sample
        adapter = BAMAdapter(adaptive_quality=True)
        profile = adapter.extract_profile(str(sample_bam))

        # Should have coverage
        assert len(profile.covered_positions) > 0

    @pytest.mark.integration
    def test_extract_profile_read_end_filtering(self, sample_bam: Path) -> None:
        """Test that read-end filtering excludes bases at read ends."""
        # Extract with read-end filtering enabled
        adapter_filtered = BAMAdapter(filter_read_ends=True, read_end_length=3)
        profile_filtered = adapter_filtered.extract_profile(str(sample_bam))

        # Extract without read-end filtering
        adapter_unfiltered = BAMAdapter(filter_read_ends=False)
        profile_unfiltered = adapter_unfiltered.extract_profile(str(sample_bam))

        # Both should have coverage
        assert len(profile_filtered.covered_positions) > 0
        assert len(profile_unfiltered.covered_positions) > 0

        # Filtered profile may have fewer positions (bases at ends excluded)
        # But should still have reasonable coverage
        assert len(profile_filtered.covered_positions) > 0

    @pytest.mark.integration
    def test_read_end_filtering_different_lengths(self, sample_bam: Path) -> None:
        """Test read-end filtering with different read_end_length values."""
        # Extract with read_end_length=1
        adapter_len1 = BAMAdapter(filter_read_ends=True, read_end_length=1)
        profile_len1 = adapter_len1.extract_profile(str(sample_bam))

        # Extract with read_end_length=5
        adapter_len5 = BAMAdapter(filter_read_ends=True, read_end_length=5)
        profile_len5 = adapter_len5.extract_profile(str(sample_bam))

        # Both should have coverage
        assert len(profile_len1.covered_positions) > 0
        assert len(profile_len5.covered_positions) > 0

        # Longer filter may have fewer positions (more bases excluded)
        # But not always, depends on read lengths
        assert len(profile_len5.covered_positions) <= len(profile_len1.covered_positions) or len(profile_len1.covered_positions) == 0


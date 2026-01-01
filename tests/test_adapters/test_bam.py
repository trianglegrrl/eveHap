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

    def test_init_custom(self) -> None:
        """Test custom initialization."""
        adapter = BAMAdapter(
            reference_path="/path/to/ref.fa",
            min_base_quality=30,
            min_mapping_quality=40,
        )
        assert adapter.reference_path == "/path/to/ref.fa"
        assert adapter.min_base_quality == 30
        assert adapter.min_mapping_quality == 40


class TestBAMAdapterExtractProfile:
    """Tests for BAMAdapter.extract_profile()."""

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

    def test_extract_profile_coverage(self, sample_bam: Path) -> None:
        """Test that extracted profile has reasonable coverage."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        # mtDNA is 16569 bp, well-sequenced samples should cover most
        coverage = profile.coverage_fraction()
        assert coverage > 0.5, f"Coverage too low: {coverage:.1%}"

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

    def test_extract_profile_variants(self, sample_bam: Path) -> None:
        """Test that extracted profile can report variants."""
        adapter = BAMAdapter()

        profile = adapter.extract_profile(str(sample_bam))

        # Most samples have at least a few variants
        variants = profile.get_variants()
        # Don't assert specific count, just that we can get variants
        assert isinstance(variants, list)

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


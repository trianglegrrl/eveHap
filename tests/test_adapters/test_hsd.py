"""Tests for evehap.adapters.hsd module."""

from pathlib import Path

import pytest

from evehap.adapters.hsd import HSDAdapter, parse_polymorphism
from evehap.core.profile import AlleleProfile


class TestParsePolymorphism:
    """Tests for polymorphism parsing."""

    def test_parse_substitution(self) -> None:
        """Test parsing simple substitution."""
        pos, ref, alt, ptype = parse_polymorphism("73G")
        assert pos == 73
        assert alt == "G"
        assert ptype == "substitution"

    def test_parse_substitution_with_ref(self) -> None:
        """Test parsing substitution with reference base."""
        pos, ref, alt, ptype = parse_polymorphism("A73G")
        assert pos == 73
        assert ref == "A"
        assert alt == "G"
        assert ptype == "substitution"

    def test_parse_insertion(self) -> None:
        """Test parsing insertion (e.g., 315.1C)."""
        pos, ref, alt, ptype = parse_polymorphism("315.1C")
        assert pos == 315
        assert alt == "C"
        assert ptype == "insertion"

    def test_parse_multi_insertion(self) -> None:
        """Test parsing multi-base insertion (e.g., 309.1CC)."""
        pos, ref, alt, ptype = parse_polymorphism("309.1CC")
        assert pos == 309
        assert alt == "CC"
        assert ptype == "insertion"

    def test_parse_deletion_lowercase(self) -> None:
        """Test parsing deletion with lowercase d."""
        pos, ref, alt, ptype = parse_polymorphism("523d")
        assert pos == 523
        assert ptype == "deletion"

    def test_parse_deletion_uppercase(self) -> None:
        """Test parsing deletion with DEL."""
        pos, ref, alt, ptype = parse_polymorphism("523DEL")
        assert pos == 523
        assert ptype == "deletion"

    def test_parse_deletion_hyphen(self) -> None:
        """Test parsing deletion with hyphen notation."""
        pos, ref, alt, ptype = parse_polymorphism("8281-8289d")
        assert pos == 8281
        assert ptype == "deletion"


class TestHSDAdapter:
    """Tests for HSDAdapter class."""

    def test_can_handle_hsd(self) -> None:
        """Test can_handle returns True for .hsd files."""
        adapter = HSDAdapter()
        assert adapter.can_handle("sample.hsd") is True
        assert adapter.can_handle("/path/to/sample.hsd") is True

    def test_can_handle_other(self) -> None:
        """Test can_handle returns False for other files."""
        adapter = HSDAdapter()
        assert adapter.can_handle("sample.bam") is False
        assert adapter.can_handle("sample.vcf") is False
        assert adapter.can_handle("sample.fasta") is False

    def test_format_name(self) -> None:
        """Test format_name returns 'HSD'."""
        adapter = HSDAdapter()
        assert adapter.format_name == "HSD"


class TestHSDAdapterExtractProfile:
    """Tests for HSDAdapter.extract_profile()."""

    @pytest.fixture
    def haplogrep_hsd(self) -> Path:
        """Return Haplogrep example HSD file."""
        # Use the example HSD from haplogrep3
        hsd_path = Path("/home/a/mtdna_poc/haplogrep3/data/examples/evaluation-data.hsd")
        if not hsd_path.exists():
            pytest.skip("Haplogrep HSD file not available")
        return hsd_path

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    def test_extract_profile_basic(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test basic profile extraction from HSD."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(haplogrep_hsd))

        assert isinstance(profile, AlleleProfile)
        assert profile.source_format == "HSD"
        # HSD only contains variants, so covered positions are variant positions
        assert len(profile.covered_positions) > 0

    def test_extract_profile_sample_id(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test that sample ID is extracted from HSD."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(haplogrep_hsd))

        assert profile.sample_id is not None
        # First sample in evaluation-data.hsd starts with "Africa"
        assert profile.sample_id.startswith("Africa")

    def test_extract_profile_specific_sample(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test extracting a specific sample by ID."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        # Extract Africa01 specifically
        profile = adapter.extract_profile(str(haplogrep_hsd), sample_id="Africa01")

        assert profile.sample_id == "Africa01"

    def test_extract_profile_variants(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test that variants are correctly extracted."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(haplogrep_hsd))

        variants = profile.get_variants()
        assert isinstance(variants, list)
        assert len(variants) > 0

    def test_extract_profile_file_not_found(self, reference_fasta: Path) -> None:
        """Test FileNotFoundError for missing file."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        with pytest.raises(FileNotFoundError):
            adapter.extract_profile("/nonexistent/file.hsd")

    def test_extract_profile_sample_not_found(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test ValueError for non-existent sample."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        with pytest.raises(ValueError, match="Sample .* not found"):
            adapter.extract_profile(str(haplogrep_hsd), sample_id="NONEXISTENT")

    def test_list_samples(
        self, haplogrep_hsd: Path, reference_fasta: Path
    ) -> None:
        """Test listing samples in HSD file."""
        adapter = HSDAdapter(reference_path=str(reference_fasta))

        samples = adapter.list_samples(str(haplogrep_hsd))

        assert isinstance(samples, list)
        assert len(samples) > 0
        # evaluation-data.hsd has Africa and Asia samples
        assert "Africa01" in samples
        assert "Asia01" in samples


class TestHSDAdapterWithReference:
    """Tests for HSD adapter with reference filling."""

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    def test_fill_reference(self, reference_fasta: Path, tmp_path: Path) -> None:
        """Test filling in reference positions."""
        # Create minimal HSD with just one variant
        hsd_file = tmp_path / "test.hsd"
        with open(hsd_file, "w") as f:
            f.write("ID\tRange\tHaplogroup\tPolymorphisms\n")
            f.write("test_sample\t1-16569\t?\t73G\n")

        adapter = HSDAdapter(reference_path=str(reference_fasta), fill_reference=True)
        profile = adapter.extract_profile(str(hsd_file))

        # Should have many positions (filled from reference)
        assert profile.coverage_fraction() > 0.95

        # Position 73 should be G (variant)
        obs = profile.observations.get(73)
        assert obs is not None
        assert "G" in obs.alleles



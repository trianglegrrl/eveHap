"""Tests for evehap.adapters.vcf module."""

from pathlib import Path
from typing import Optional

import pytest

from evehap.adapters.vcf import VCFAdapter
from evehap.core.profile import AlleleProfile


class TestVCFAdapter:
    """Tests for VCFAdapter class."""

    def test_can_handle_vcf(self) -> None:
        """Test can_handle returns True for .vcf files."""
        adapter = VCFAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.vcf") is True
        assert adapter.can_handle("/path/to/sample.vcf") is True

    def test_can_handle_vcf_gz(self) -> None:
        """Test can_handle returns True for .vcf.gz files."""
        adapter = VCFAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.vcf.gz") is True
        assert adapter.can_handle("/path/to/sample.vcf.gz") is True

    def test_can_handle_other(self) -> None:
        """Test can_handle returns False for other files."""
        adapter = VCFAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.bam") is False
        assert adapter.can_handle("sample.fasta") is False
        assert adapter.can_handle("sample.vcf.tbi") is False

    def test_format_name(self) -> None:
        """Test format_name returns 'VCF'."""
        adapter = VCFAdapter(reference_path="/path/to/ref.fa")
        assert adapter.format_name == "VCF"

    def test_init_defaults(self) -> None:
        """Test default initialization."""
        adapter = VCFAdapter(reference_path="/path/to/ref.fa")
        assert adapter.reference_path == "/path/to/ref.fa"
        assert adapter.sample_id is None
        assert adapter.mt_contig == "MT"

    def test_init_custom(self) -> None:
        """Test custom initialization."""
        adapter = VCFAdapter(
            reference_path="/path/to/ref.fa",
            sample_id="HG00096",
            mt_contig="chrM",
        )
        assert adapter.reference_path == "/path/to/ref.fa"
        assert adapter.sample_id == "HG00096"
        assert adapter.mt_contig == "chrM"


class TestVCFAdapterExtractProfile:
    """Tests for VCFAdapter.extract_profile()."""

    @pytest.fixture
    def onekg_vcf(self, mtdna_poc_data_dir: Path) -> Path:
        """Return 1000 Genomes mtDNA VCF."""
        vcf_path = mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"
        if not vcf_path.exists():
            pytest.skip("1000 Genomes VCF not available")
        return vcf_path

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    def test_extract_profile_basic(self, onekg_vcf: Path, reference_fasta: Path) -> None:
        """Test basic profile extraction from VCF."""
        adapter = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="HG00096",
        )

        profile = adapter.extract_profile(str(onekg_vcf))

        assert isinstance(profile, AlleleProfile)
        assert profile.sample_id == "HG00096"
        assert profile.source_format == "VCF"
        # VCF should have coverage at most mtDNA positions (from reference)
        assert len(profile.covered_positions) > 0

    def test_extract_profile_first_sample(self, onekg_vcf: Path, reference_fasta: Path) -> None:
        """Test extracting first sample when no sample_id specified."""
        adapter = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id=None,  # Use first sample
        )

        profile = adapter.extract_profile(str(onekg_vcf))

        assert isinstance(profile, AlleleProfile)
        # Sample ID should be the first sample in the VCF
        assert profile.sample_id == "HG00096"

    def test_extract_profile_variants(self, onekg_vcf: Path, reference_fasta: Path) -> None:
        """Test that extracted profile can report variants."""
        adapter = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="HG00096",
        )

        profile = adapter.extract_profile(str(onekg_vcf))

        variants = profile.get_variants()
        assert isinstance(variants, list)
        # Most samples have at least some variants

    def test_extract_profile_different_samples(
        self, onekg_vcf: Path, reference_fasta: Path
    ) -> None:
        """Test that different samples produce different profiles."""
        adapter1 = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="HG00096",
        )
        adapter2 = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="HG00097",
        )

        profile1 = adapter1.extract_profile(str(onekg_vcf))
        profile2 = adapter2.extract_profile(str(onekg_vcf))

        # Different samples may have different variants
        variants1 = set(
            (v["position"], v["alt_allele"])
            for v in profile1.get_variants()
        )
        variants2 = set(
            (v["position"], v["alt_allele"])
            for v in profile2.get_variants()
        )

        # At least samples should have different IDs
        assert profile1.sample_id != profile2.sample_id

    def test_extract_profile_file_not_found(self, reference_fasta: Path) -> None:
        """Test FileNotFoundError for missing file."""
        adapter = VCFAdapter(reference_path=str(reference_fasta))

        with pytest.raises(FileNotFoundError):
            adapter.extract_profile("/nonexistent/file.vcf")

    def test_extract_profile_sample_not_found(
        self, onekg_vcf: Path, reference_fasta: Path
    ) -> None:
        """Test ValueError for non-existent sample."""
        adapter = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="NONEXISTENT_SAMPLE",
        )

        with pytest.raises(ValueError, match="Sample .* not found"):
            adapter.extract_profile(str(onekg_vcf))

    def test_extract_profile_coverage(self, onekg_vcf: Path, reference_fasta: Path) -> None:
        """Test that VCF profile has full mtDNA coverage from reference."""
        adapter = VCFAdapter(
            reference_path=str(reference_fasta),
            sample_id="HG00096",
        )

        profile = adapter.extract_profile(str(onekg_vcf))

        # VCF + reference should give near-complete coverage
        coverage = profile.coverage_fraction()
        assert coverage > 0.95, f"Expected >95% coverage, got {coverage:.1%}"


class TestVCFAdapterMultiSample:
    """Tests for multi-sample VCF operations."""

    @pytest.fixture
    def onekg_vcf(self, mtdna_poc_data_dir: Path) -> Path:
        """Return 1000 Genomes mtDNA VCF."""
        vcf_path = mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"
        if not vcf_path.exists():
            pytest.skip("1000 Genomes VCF not available")
        return vcf_path

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    def test_list_samples(self, onekg_vcf: Path, reference_fasta: Path) -> None:
        """Test listing available samples in VCF."""
        adapter = VCFAdapter(reference_path=str(reference_fasta))

        samples = adapter.list_samples(str(onekg_vcf))

        assert isinstance(samples, list)
        assert len(samples) > 0
        assert "HG00096" in samples
        assert "HG00097" in samples



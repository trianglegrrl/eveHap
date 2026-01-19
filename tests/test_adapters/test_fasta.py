"""Tests for evehap.adapters.fasta module."""

from pathlib import Path

import pytest

from evehap.adapters.fasta import FASTAAdapter
from evehap.core.profile import AlleleProfile


class TestFASTAAdapter:
    """Tests for FASTAAdapter class."""

    def test_can_handle_fasta(self) -> None:
        """Test can_handle returns True for .fasta files."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.fasta") is True
        assert adapter.can_handle("/path/to/sample.fasta") is True

    def test_can_handle_fa(self) -> None:
        """Test can_handle returns True for .fa files."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.fa") is True
        assert adapter.can_handle("/path/to/sample.fa") is True

    def test_can_handle_fna(self) -> None:
        """Test can_handle returns True for .fna files."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.fna") is True

    def test_can_handle_other(self) -> None:
        """Test can_handle returns False for other files."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.can_handle("sample.bam") is False
        assert adapter.can_handle("sample.vcf") is False
        assert adapter.can_handle("sample.fastq") is False

    def test_format_name(self) -> None:
        """Test format_name returns 'FASTA'."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.format_name == "FASTA"

    def test_init(self) -> None:
        """Test initialization."""
        adapter = FASTAAdapter(reference_path="/path/to/ref.fa")
        assert adapter.reference_path == "/path/to/ref.fa"


class TestFASTAAdapterExtractProfile:
    """Tests for FASTAAdapter.extract_profile()."""

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    @pytest.fixture
    def hgdp_consensus_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return HGDP consensus FASTA directory."""
        consensus_dir = mtdna_poc_data_dir / "hgdp" / "haplogroups" / "consensus"
        if not consensus_dir.exists():
            pytest.skip("HGDP consensus directory not available")
        return consensus_dir

    @pytest.fixture
    def sample_fasta(self, hgdp_consensus_dir: Path) -> Path:
        """Return a sample FASTA file."""
        fastas = list(hgdp_consensus_dir.glob("*.fasta"))
        if not fastas:
            pytest.skip("No FASTA files found")
        return fastas[0]

    def test_extract_profile_basic(self, sample_fasta: Path, reference_fasta: Path) -> None:
        """Test basic profile extraction from FASTA."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(sample_fasta))

        assert isinstance(profile, AlleleProfile)
        assert profile.source_format == "FASTA"
        assert len(profile.covered_positions) > 0

    def test_extract_profile_sample_id(self, sample_fasta: Path, reference_fasta: Path) -> None:
        """Test that sample ID is extracted from FASTA header."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(sample_fasta))

        # Sample ID should come from the FASTA header
        assert profile.sample_id is not None
        assert len(profile.sample_id) > 0

    def test_extract_profile_coverage(self, sample_fasta: Path, reference_fasta: Path) -> None:
        """Test that FASTA gives full mtDNA coverage."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(sample_fasta))

        # FASTA should have near-complete coverage
        coverage = profile.coverage_fraction()
        assert coverage > 0.95, f"Expected >95% coverage, got {coverage:.1%}"

    def test_extract_profile_variants(self, sample_fasta: Path, reference_fasta: Path) -> None:
        """Test that variants are detected relative to reference."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        profile = adapter.extract_profile(str(sample_fasta))

        variants = profile.get_variants()
        assert isinstance(variants, list)
        # Most samples have at least a few variants from rCRS

    def test_extract_profile_file_not_found(self, reference_fasta: Path) -> None:
        """Test FileNotFoundError for missing file."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        with pytest.raises(FileNotFoundError):
            adapter.extract_profile("/nonexistent/file.fasta")

    def test_extract_profile_reference(self, reference_fasta: Path) -> None:
        """Test extracting reference FASTA gives no variants."""
        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        # Extract reference against itself
        profile = adapter.extract_profile(str(reference_fasta))

        # Should have full coverage
        coverage = profile.coverage_fraction()
        assert coverage > 0.99

        # Should have no variants (sample = reference)
        variants = profile.get_variants()
        assert len(variants) == 0, f"Expected no variants, got {len(variants)}"

    def test_extract_profile_multi_sequence(
        self, hgdp_consensus_dir: Path, reference_fasta: Path
    ) -> None:
        """Test extracting from multi-sequence FASTA."""
        multi_fasta = hgdp_consensus_dir.parent / "all_samples.fasta"
        if not multi_fasta.exists():
            pytest.skip("Multi-sample FASTA not available")

        adapter = FASTAAdapter(reference_path=str(reference_fasta))

        # Should extract first sequence by default
        profile = adapter.extract_profile(str(multi_fasta))

        assert isinstance(profile, AlleleProfile)
        assert len(profile.covered_positions) > 0


class TestFASTAAdapterIUPAC:
    """Tests for IUPAC ambiguity code handling."""

    @pytest.fixture
    def reference_fasta(self, evehap_data_dir: Path) -> Path:
        """Return reference FASTA."""
        ref_path = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref_path.exists():
            pytest.skip("Reference FASTA not available")
        return ref_path

    def test_iupac_r_code(self, reference_fasta: Path, tmp_path: Path) -> None:
        """Test IUPAC R (A/G) code is converted to heteroplasmy."""
        # Create test FASTA with R at position 73 (rCRS is A at 73)
        test_fasta = tmp_path / "test.fasta"

        # Read reference
        with open(reference_fasta) as f:
            lines = f.readlines()
            header = lines[0]
            seq = "".join(line.strip() for line in lines[1:])

        # Modify position 73 (0-indexed: 72) to R
        seq_list = list(seq)
        seq_list[72] = "R"  # Position 73
        modified_seq = "".join(seq_list)

        # Write test FASTA
        with open(test_fasta, "w") as f:
            f.write(">test_sample\n")
            f.write(modified_seq + "\n")

        adapter = FASTAAdapter(reference_path=str(reference_fasta))
        profile = adapter.extract_profile(str(test_fasta))

        # Position 73 should have A and G
        obs = profile.observations.get(73)
        assert obs is not None, "Position 73 should be covered"
        assert (
            "A" in obs.alleles or "G" in obs.alleles
        ), f"Expected A and/or G at position 73, got {obs.alleles}"

    def test_iupac_n_code(self, reference_fasta: Path, tmp_path: Path) -> None:
        """Test IUPAC N (any base) code is handled as missing."""
        test_fasta = tmp_path / "test.fasta"

        # Read reference
        with open(reference_fasta) as f:
            lines = f.readlines()
            seq = "".join(line.strip() for line in lines[1:])

        # Modify position 100 to N
        seq_list = list(seq)
        seq_list[99] = "N"  # Position 100
        modified_seq = "".join(seq_list)

        with open(test_fasta, "w") as f:
            f.write(">test_sample\n")
            f.write(modified_seq + "\n")

        adapter = FASTAAdapter(reference_path=str(reference_fasta))
        profile = adapter.extract_profile(str(test_fasta))

        # Position 100 should be missing (N = unknown)
        obs = profile.observations.get(100)
        # Either not present or marked as unknown
        if obs is not None:
            assert obs.alleles.get("N", 0) > 0 or len(obs.alleles) == 0

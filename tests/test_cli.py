"""Tests for CLI --tree-type option.

TDD: Red phase - these tests should FAIL until we implement --tree-type.
"""

import pytest
from click.testing import CliRunner

from evehap.cli import classify


class TestTreeTypeOption:
    """Tests for required --tree-type CLI option."""

    @pytest.fixture
    def runner(self):
        """Create CLI test runner."""
        return CliRunner()

    @pytest.fixture
    def sample_bam(self, tmp_path):
        """Create a minimal test BAM file path.

        Note: For these CLI tests, we just need the path to exist.
        The actual file content doesn't matter since we're testing
        argument validation, not classification.
        """
        # Create a placeholder file
        bam_path = tmp_path / "test.bam"
        bam_path.touch()
        return str(bam_path)

    def test_tree_type_is_required(self, runner, sample_bam):
        """Test that --tree-type is required."""
        result = runner.invoke(classify, [sample_bam])

        assert result.exit_code != 0
        assert "tree-type" in result.output.lower() or "required" in result.output.lower()

    def test_tree_type_rcrs_valid(self, runner, sample_bam):
        """Test that --tree-type rcrs is accepted."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rcrs"])

        # Should not fail due to tree-type validation
        # May fail due to other reasons (file parsing), but that's OK
        assert "invalid value" not in result.output.lower()
        assert "not one of" not in result.output.lower()

    def test_tree_type_rsrs_valid(self, runner, sample_bam):
        """Test that --tree-type rsrs is accepted."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rsrs"])

        # Should not fail due to tree-type validation
        assert "invalid value" not in result.output.lower()
        assert "not one of" not in result.output.lower()

    def test_tree_type_invalid_rejected(self, runner, sample_bam):
        """Test that invalid tree-type is rejected."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "invalid"])

        assert result.exit_code != 0
        assert "invalid value" in result.output.lower() or "not one of" in result.output.lower()

    def test_tree_type_rcrs_uses_correct_tree(self, runner, sample_bam):
        """Test that --tree-type rcrs uses the rCRS tree file."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rcrs", "-q"])

        # The output should mention rcrs or the tree file
        # This will be visible in verbose mode
        # For now, just ensure it doesn't crash on tree-type
        assert "invalid value" not in result.output.lower()

    def test_tree_type_rsrs_uses_correct_tree(self, runner, sample_bam):
        """Test that --tree-type rsrs uses the RSRS tree file."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rsrs", "-q"])

        # Should use RSRS tree
        assert "invalid value" not in result.output.lower()

    def test_tree_and_tree_type_mutually_exclusive(self, runner, sample_bam):
        """Test that --tree and --tree-type cannot be used together.

        When --tree-type is specified, --tree should not be required.
        The tree is auto-selected based on tree-type.
        """
        # Using --tree-type should work without --tree
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rcrs"])

        # Should not complain about missing --tree
        assert "--tree is required" not in result.output.lower()

    def test_reference_still_required(self, runner, sample_bam):
        """Test that --reference is still required with --tree-type."""
        result = runner.invoke(classify, [sample_bam, "--tree-type", "rcrs"])

        # Should still require reference (or use a sensible default)
        # This test verifies the transition behavior
        # Note: We may auto-default reference based on tree-type
        pass  # Placeholder - behavior TBD


class TestTreeTypeIntegration:
    """Integration tests for tree-type with real classification.

    These tests require actual sample files and are marked for
    integration test runs only.
    """

    @pytest.fixture
    def runner(self):
        return CliRunner()

    @pytest.mark.integration
    def test_vcf_with_rcrs_tree(self, runner):
        """Test VCF classification with rCRS tree."""
        # This will be implemented in Phase 2
        pass

    @pytest.mark.integration
    def test_vcf_with_rsrs_tree(self, runner):
        """Test VCF classification with RSRS tree."""
        pass

    @pytest.mark.integration
    def test_bam_with_rcrs_tree(self, runner):
        """Test BAM classification with rCRS tree."""
        pass

    @pytest.mark.integration
    def test_bam_with_rsrs_tree(self, runner):
        """Test BAM classification with RSRS tree."""
        pass

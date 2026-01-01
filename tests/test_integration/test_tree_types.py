"""Integration tests for tree-type x input format matrix.

Tests that both VCF and BAM inputs work correctly with both rcrs and rsrs trees,
and that results match Haplogrep3 when using the corresponding tree.

TDD: These tests define expected behavior for Phase 2.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pytest

# Test configuration
VCF_PATH = "/home/a/mtdna_poc/data/raw/1kg_chrMT.vcf.gz"
HGDP_BAM_DIR = "/home/a/mtdna_poc/data/hgdp/bams"
HAPLOGREP3_CMD = "/home/a/anaconda3/envs/haplogrep3/bin/haplogrep3"

# Skip if test data not available
pytestmark = pytest.mark.skipif(
    not Path(VCF_PATH).exists(),
    reason="Test data not available"
)


def run_haplogrep3(variants: list, tree_type: str) -> str:
    """Run Haplogrep3 on a list of variants and return haplogroup.

    Args:
        variants: List of variant strings like "73G", "263G"
        tree_type: 'rcrs' or 'rsrs'

    Returns:
        Haplogroup called by Haplogrep3
    """
    # Map tree type to Haplogrep3 tree name
    hp3_tree = {
        "rcrs": "phylotree-fu-rcrs@1.2",
        "rsrs": "phylotree-rsrs@17.1",
    }[tree_type]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.hsd', delete=False) as f:
        hsd_line = f"sample\t1-16569\t?\t" + "\t".join(variants)
        f.write(hsd_line + "\n")
        hsd_path = f.name

    result_path = hsd_path.replace('.hsd', '.txt')

    try:
        subprocess.run([
            HAPLOGREP3_CMD,
            "classify",
            "--tree", hp3_tree,
            "--input", hsd_path,
            "--output", result_path,
        ], capture_output=True, text=True, check=True)

        with open(result_path) as f:
            lines = f.readlines()
            if len(lines) > 1:
                parts = lines[1].strip().replace('"', '').split('\t')
                return parts[1] if len(parts) > 1 else "?"
    except Exception:
        return "error"
    finally:
        os.unlink(hsd_path)
        if os.path.exists(result_path):
            os.unlink(result_path)

    return "?"


class TestVCFWithTreeTypes:
    """Test VCF classification with both tree types."""

    @pytest.fixture
    def evehap_classify(self):
        """Return function to classify a sample using evehap."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        import pysam
        from evehap.adapters.vcf import VCFAdapter
        from evehap.core.phylotree import Phylotree
        from evehap.core.classifier import Classifier
        from evehap.cli import TREE_TYPES

        def classify(sample_id: str, tree_type: str) -> tuple:
            """Classify sample and return (haplogroup, variants)."""
            config = TREE_TYPES[tree_type]

            with pysam.FastaFile(str(config["reference"])) as fasta:
                reference = fasta.fetch(fasta.references[0])

            phylotree = Phylotree.load(str(config["tree"]), None)
            classifier = Classifier(
                phylotree,
                reference=reference,
                weights_path=str(config["weights"]) if config["weights"].exists() else None,
            )

            adapter = VCFAdapter(
                reference_path=str(config["reference"]),
                sample_id=sample_id,
            )
            profile = adapter.extract_profile(VCF_PATH)
            result = classifier.classify(profile)

            # Get variants for Haplogrep3 comparison
            variants = [
                f"{pos}{obs.major_allele}"
                for pos, obs in sorted(profile.observations.items())
                if obs.is_variant and obs.major_allele != "-"
            ]

            return result.haplogroup, variants

        return classify

    @pytest.mark.integration
    def test_vcf_rcrs_matches_haplogrep3(self, evehap_classify):
        """Test VCF with rCRS tree matches Haplogrep3."""
        samples = ["HG00096", "HG00097", "HG00099"]

        for sample_id in samples:
            evehap_hg, variants = evehap_classify(sample_id, "rcrs")
            hp3_hg = run_haplogrep3(variants, "rcrs")

            assert evehap_hg == hp3_hg, \
                f"Sample {sample_id}: evehap={evehap_hg}, haplogrep3={hp3_hg}"

    @pytest.mark.integration
    def test_vcf_rsrs_works(self, evehap_classify):
        """Test VCF with RSRS tree produces valid results."""
        samples = ["HG00096", "HG00097"]

        for sample_id in samples:
            evehap_hg, _ = evehap_classify(sample_id, "rsrs")

            # Just verify it produces a valid haplogroup
            assert evehap_hg is not None
            assert evehap_hg != "unknown"
            assert len(evehap_hg) > 0


class TestBAMWithTreeTypes:
    """Test BAM classification with both tree types."""

    @pytest.fixture
    def evehap_classify_bam(self):
        """Return function to classify a BAM file using evehap."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        import pysam
        from evehap.adapters.bam import BAMAdapter
        from evehap.core.phylotree import Phylotree
        from evehap.core.classifier import Classifier
        from evehap.cli import TREE_TYPES

        def classify(bam_path: str, tree_type: str) -> str:
            """Classify BAM file and return haplogroup."""
            config = TREE_TYPES[tree_type]

            with pysam.FastaFile(str(config["reference"])) as fasta:
                reference = fasta.fetch(fasta.references[0])

            phylotree = Phylotree.load(str(config["tree"]), None)
            classifier = Classifier(
                phylotree,
                reference=reference,
                weights_path=str(config["weights"]) if config["weights"].exists() else None,
            )

            adapter = BAMAdapter(reference_path=str(config["reference"]))
            profile = adapter.extract_profile(bam_path)
            result = classifier.classify(profile)

            return result.haplogroup

        return classify

    @pytest.mark.integration
    @pytest.mark.skipif(
        not Path(HGDP_BAM_DIR).exists(),
        reason="HGDP BAM files not available"
    )
    def test_bam_rcrs_works(self, evehap_classify_bam):
        """Test BAM with rCRS tree produces valid results."""
        import glob

        bam_files = glob.glob(os.path.join(HGDP_BAM_DIR, "*.bam"))[:3]

        for bam_path in bam_files:
            evehap_hg = evehap_classify_bam(bam_path, "rcrs")

            assert evehap_hg is not None
            assert evehap_hg != "unknown"
            assert len(evehap_hg) > 0

    @pytest.mark.integration
    @pytest.mark.skipif(
        not Path(HGDP_BAM_DIR).exists(),
        reason="HGDP BAM files not available"
    )
    def test_bam_rsrs_works(self, evehap_classify_bam):
        """Test BAM with RSRS tree produces valid results."""
        import glob

        bam_files = glob.glob(os.path.join(HGDP_BAM_DIR, "*.bam"))[:3]

        for bam_path in bam_files:
            evehap_hg = evehap_classify_bam(bam_path, "rsrs")

            assert evehap_hg is not None
            assert evehap_hg != "unknown"
            assert len(evehap_hg) > 0


class TestTreeTypeConsistency:
    """Test that tree types produce consistent results."""

    @pytest.fixture
    def classify_both(self):
        """Return function to classify with both tree types."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        import pysam
        from evehap.adapters.vcf import VCFAdapter
        from evehap.core.phylotree import Phylotree
        from evehap.core.classifier import Classifier
        from evehap.cli import TREE_TYPES

        def classify(sample_id: str) -> dict:
            """Classify with both tree types."""
            results = {}

            for tree_type in ["rcrs", "rsrs"]:
                config = TREE_TYPES[tree_type]

                with pysam.FastaFile(str(config["reference"])) as fasta:
                    reference = fasta.fetch(fasta.references[0])

                phylotree = Phylotree.load(str(config["tree"]), None)
                classifier = Classifier(
                    phylotree,
                    reference=reference,
                    weights_path=str(config["weights"]) if config["weights"].exists() else None,
                )

                adapter = VCFAdapter(
                    reference_path=str(config["reference"]),
                    sample_id=sample_id,
                )
                profile = adapter.extract_profile(VCF_PATH)
                result = classifier.classify(profile)
                results[tree_type] = result.haplogroup

            return results

        return classify

    @pytest.mark.integration
    def test_major_clade_consistent(self, classify_both):
        """Test that both trees agree on major clade (first letter)."""
        sample_id = "HG00096"

        results = classify_both(sample_id)

        # Both should start with 'H' for this sample
        assert results["rcrs"][0] == results["rsrs"][0], \
            f"Major clade mismatch: rcrs={results['rcrs']}, rsrs={results['rsrs']}"


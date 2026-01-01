"""Integration tests for microarray adapter with 23andMe format.

Tests the round-trip accuracy: VCF -> 23andMe -> classify
Verifies that classification results match between formats.
Also tests real 23andMe files from /home/a/mtdna_poc/data/23andme/
"""

import os
import subprocess
import tempfile
from pathlib import Path

import pytest

# Test configuration
VCF_PATH = "/home/a/mtdna_poc/data/raw/1kg_chrMT.vcf.gz"
RSID_LOOKUP_PATH = Path(__file__).parent.parent.parent / "data" / "dbsnp" / "mtdna_rsids.tsv"

# Real 23andMe test files
REAL_23ANDME_DIR = Path("/home/a/mtdna_poc/data/23andme")
REAL_23ANDME_FILES = [
    REAL_23ANDME_DIR / "a.txt",
    REAL_23ANDME_DIR / "b.txt",
]

# Skip if test data not available
pytestmark = pytest.mark.skipif(
    not Path(VCF_PATH).exists(),
    reason="1KG VCF test data not available"
)


class TestMicroarrayRoundTrip:
    """Test VCF -> 23andMe -> classify round-trip accuracy."""

    @pytest.fixture
    def vcf_classifier(self):
        """Create classifier for VCF data."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        import pysam
        from evehap.adapters.vcf import VCFAdapter
        from evehap.core.phylotree import Phylotree
        from evehap.core.classifier import Classifier
        from evehap.cli import TREE_TYPES

        tree_type = "rcrs"
        config = TREE_TYPES[tree_type]

        with pysam.FastaFile(str(config["reference"])) as fasta:
            reference = fasta.fetch(fasta.references[0])

        phylotree = Phylotree.load(str(config["tree"]), None)
        weights_path = str(config["weights"]) if config["weights"].exists() else None
        classifier = Classifier(phylotree, reference=reference, weights_path=weights_path)

        return classifier, config, reference

    @pytest.fixture
    def generate_23andme(self):
        """Return function to generate 23andMe file from VCF sample."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

        from generate_23andme_from_vcf import vcf_to_23andme, write_23andme_file, load_rsid_lookup

        rsid_lookup = {}
        if RSID_LOOKUP_PATH.exists():
            rsid_lookup = load_rsid_lookup(str(RSID_LOOKUP_PATH))

        def generate(sample_id: str) -> str:
            """Generate 23andMe file and return path."""
            lines = vcf_to_23andme(VCF_PATH, sample_id, rsid_lookup)

            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
                f.write('\n'.join(lines) + '\n')
                return f.name

        return generate

    @pytest.mark.integration
    def test_single_sample_round_trip(self, vcf_classifier, generate_23andme):
        """Test round-trip accuracy for a single sample."""
        from evehap.adapters.vcf import VCFAdapter
        from evehap.adapters.microarray import MicroarrayAdapter

        classifier, config, reference = vcf_classifier
        sample_id = "HG00096"

        # Classify from VCF
        vcf_adapter = VCFAdapter(
            reference_path=str(config["reference"]),
            sample_id=sample_id,
        )
        vcf_profile = vcf_adapter.extract_profile(VCF_PATH)
        vcf_result = classifier.classify(vcf_profile)

        # Generate and classify from 23andMe
        ma_path = generate_23andme(sample_id)
        try:
            ma_adapter = MicroarrayAdapter()
            ma_profile = ma_adapter.extract_profile(ma_path)
            ma_profile.assumes_full_coverage = True
            ma_result = classifier.classify(ma_profile)

            assert vcf_result.haplogroup == ma_result.haplogroup, \
                f"Round-trip mismatch for {sample_id}: VCF={vcf_result.haplogroup}, 23andMe={ma_result.haplogroup}"
        finally:
            os.unlink(ma_path)

    @pytest.mark.integration
    def test_multiple_samples_round_trip(self, vcf_classifier, generate_23andme):
        """Test round-trip accuracy for multiple samples."""
        from evehap.adapters.vcf import VCFAdapter
        from evehap.adapters.microarray import MicroarrayAdapter

        classifier, config, reference = vcf_classifier
        samples = ["HG00096", "HG00097", "HG00099"]

        for sample_id in samples:
            # Classify from VCF
            vcf_adapter = VCFAdapter(
                reference_path=str(config["reference"]),
                sample_id=sample_id,
            )
            vcf_profile = vcf_adapter.extract_profile(VCF_PATH)
            vcf_result = classifier.classify(vcf_profile)

            # Generate and classify from 23andMe
            ma_path = generate_23andme(sample_id)
            try:
                ma_adapter = MicroarrayAdapter()
                ma_profile = ma_adapter.extract_profile(ma_path)
                ma_profile.assumes_full_coverage = True
                ma_result = classifier.classify(ma_profile)

                assert vcf_result.haplogroup == ma_result.haplogroup, \
                    f"Round-trip mismatch for {sample_id}: VCF={vcf_result.haplogroup}, 23andMe={ma_result.haplogroup}"
            finally:
                os.unlink(ma_path)


class TestMicroarrayParser:
    """Test microarray file parsing."""

    @pytest.fixture
    def sample_23andme_file(self, tmp_path):
        """Create a sample 23andMe format file."""
        content = """# This data file generated by 23andMe
# rsid	chromosome	position	genotype
rs3087742	MT	73	G
rs2853499	MT	263	G
rs2853503	MT	750	G
rs2853506	MT	1438	G
rs3021086	MT	4769	G
rs2001031	MT	8860	G
rs2853508	MT	15326	G
rs3937033	MT	16519	C
"""
        file_path = tmp_path / "test.23andme.txt"
        file_path.write_text(content)
        return str(file_path)

    def test_parse_23andme_format(self, sample_23andme_file):
        """Test parsing 23andMe format file."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        from evehap.adapters.microarray import MicroarrayAdapter

        adapter = MicroarrayAdapter()
        profile = adapter.extract_profile(sample_23andme_file)

        assert len(profile.observations) == 8
        assert 73 in profile.observations
        assert profile.observations[73].major_allele == "G"
        assert 16519 in profile.observations
        assert profile.observations[16519].major_allele == "C"

    def test_detect_format(self, sample_23andme_file):
        """Test format detection for 23andMe files."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        from evehap.adapters.microarray import MicroarrayAdapter

        adapter = MicroarrayAdapter()
        detected = adapter._detect_format(sample_23andme_file)

        assert detected == "23andme"


class TestReal23andMeFiles:
    """Test classification with real 23andMe files.

    Uses actual 23andMe data from /home/a/mtdna_poc/data/23andme/
    to verify classification matches Haplogrep3.
    """

    @pytest.fixture
    def classifier_rcrs(self):
        """Create rCRS classifier."""
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        import pysam
        from evehap.core.phylotree import Phylotree
        from evehap.core.classifier import Classifier
        from evehap.cli import TREE_TYPES

        config = TREE_TYPES["rcrs"]

        with pysam.FastaFile(str(config["reference"])) as fasta:
            reference = fasta.fetch(fasta.references[0])

        phylotree = Phylotree.load(str(config["tree"]), None)
        weights_path = str(config["weights"]) if config["weights"].exists() else None
        classifier = Classifier(phylotree, reference=reference, weights_path=weights_path)

        return classifier, reference

    @pytest.mark.integration
    @pytest.mark.skipif(
        not all(f.exists() for f in REAL_23ANDME_FILES),
        reason="Real 23andMe test files not available"
    )
    def test_real_23andme_files_classify(self, classifier_rcrs):
        """Test that real 23andMe files can be classified."""
        from evehap.adapters.microarray import MicroarrayAdapter

        classifier, reference = classifier_rcrs
        adapter = MicroarrayAdapter()

        for file_path in REAL_23ANDME_FILES:
            profile = adapter.extract_profile(str(file_path))
            profile.assumes_full_coverage = True

            result = classifier.classify(profile)

            # Should produce a valid haplogroup
            assert result.haplogroup is not None
            assert len(result.haplogroup) > 0
            assert result.haplogroup != "unknown"

            # Should have reasonable confidence
            assert result.confidence > 0

            # Should extract observations
            assert len(profile.observations) > 1000  # 23andMe chips have 2000+ mtDNA positions

    @pytest.mark.integration
    @pytest.mark.skipif(
        not all(f.exists() for f in REAL_23ANDME_FILES),
        reason="Real 23andMe test files not available"
    )
    def test_real_23andme_matches_haplogrep3(self, classifier_rcrs):
        """Test that real 23andMe classification matches Haplogrep3."""
        import subprocess
        import tempfile
        import pysam
        from evehap.adapters.microarray import MicroarrayAdapter
        from evehap.cli import TREE_TYPES

        classifier, reference = classifier_rcrs
        adapter = MicroarrayAdapter()
        config = TREE_TYPES["rcrs"]

        with pysam.FastaFile(str(config["reference"])) as fasta:
            ref_seq = fasta.fetch(fasta.references[0])

        for file_path in REAL_23ANDME_FILES:
            profile = adapter.extract_profile(str(file_path))
            profile.assumes_full_coverage = True

            # Classify with evehap
            result = classifier.classify(profile)

            # Extract variants from rCRS for Haplogrep3
            variants = []
            for pos, obs in sorted(profile.observations.items()):
                if pos <= len(ref_seq):
                    ref_allele = ref_seq[pos - 1]
                    if obs.major_allele and obs.major_allele != ref_allele:
                        variants.append(f"{pos}{obs.major_allele}")

            # Run Haplogrep3
            hsd_line = f"sample\t1-16569\t?\t" + "\t".join(variants)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.hsd', delete=False) as f:
                f.write(hsd_line + "\n")
                hsd_path = f.name

            result_path = hsd_path.replace('.hsd', '.txt')
            try:
                subprocess.run([
                    "/home/a/anaconda3/envs/haplogrep3/bin/haplogrep3",
                    "classify",
                    "--tree", "phylotree-fu-rcrs@1.2",
                    "--input", hsd_path,
                    "--output", result_path,
                ], capture_output=True, check=True)

                with open(result_path) as f:
                    lines = f.readlines()
                    if len(lines) > 1:
                        hp3_hg = lines[1].strip().replace('"', '').split('\t')[1]

                        # Should match Haplogrep3
                        assert result.haplogroup == hp3_hg, \
                            f"Mismatch for {file_path.name}: evehap={result.haplogroup}, haplogrep3={hp3_hg}"
            finally:
                os.unlink(hsd_path)
                if os.path.exists(result_path):
                    os.unlink(result_path)


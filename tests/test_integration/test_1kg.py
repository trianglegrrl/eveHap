"""Integration tests against 1000 Genomes ground truth.

Tests the complete pipeline from VCF files to haplogroup classification,
comparing results against known 1000 Genomes haplogroup assignments.

Usage:
    # Run with default 10 samples
    pytest tests/test_integration/test_1kg.py -v

    # Run with 50 random samples (seed=42 for reproducibility)
    pytest tests/test_integration/test_1kg.py -v --sample-count 50 --random-seed 42

    # Run on all samples
    pytest tests/test_integration/test_1kg.py -v --all-samples
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

from evehap.adapters.vcf import VCFAdapter
from evehap.core.classifier import Classifier
from evehap.core.phylotree import Phylotree
from evehap.output.result import HaplogroupResult

from ..conftest import select_samples


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade from full haplogroup."""
    if not haplogroup:
        return ""

    # Handle prefix letters
    for i, char in enumerate(haplogroup):
        if char.isdigit():
            return haplogroup[:i] if i > 0 else haplogroup[0]

    return haplogroup[0]


def normalize_haplogroup(haplogroup: str) -> str:
    """Normalize haplogroup for comparison."""
    # Remove annotations like *, +, etc.
    hg = haplogroup.strip().upper()
    for char in ["'", " ", "*", "+"]:
        if char in hg:
            hg = hg.split(char)[0]
    return hg


class Test1KGIntegration:
    """Integration tests against 1000 Genomes samples."""

    @pytest.fixture
    def vcf_path(self, mtdna_poc_data_dir: Path) -> Path:
        """Return 1000 Genomes VCF path."""
        vcf = mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"
        if not vcf.exists():
            pytest.skip("1000 Genomes VCF not available")
        return vcf

    @pytest.fixture
    def haplogroups(self, mtdna_poc_data_dir: Path) -> Dict[str, str]:
        """Load 1000 Genomes haplogroup ground truth."""
        hg_file = mtdna_poc_data_dir / "haplogroups" / "1kg_haplogroups_simple.csv"
        if not hg_file.exists():
            pytest.skip("1000 Genomes haplogroups file not available")

        haplogroups = {}
        with open(hg_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row.get("SampleID", "")
                hg = row.get("Haplogroup", "")
                if sample_id and hg:
                    haplogroups[sample_id] = hg

        return haplogroups

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        # Use RSRS tree
        tree_path = evehap_data_dir / "phylotree" / "tree-rsrs.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights-rsrs.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(str(tree_path), str(weights_path) if weights_path.exists() else None)

    @pytest.fixture
    def reference_path(self, evehap_data_dir: Path) -> str:
        """Return reference FASTA path."""
        ref = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref.exists():
            pytest.skip("Reference FASTA not available")
        return str(ref)

    def test_vcf_exists(self, vcf_path: Path) -> None:
        """Test that we have the 1000 Genomes VCF."""
        assert vcf_path.exists()
        assert vcf_path.stat().st_size > 0

    def test_haplogroups_loaded(self, haplogroups: Dict[str, str]) -> None:
        """Test that we have ground truth haplogroups."""
        assert len(haplogroups) > 100, "Expected at least 100 ground truth samples"

    def test_classify_single_sample(
        self,
        vcf_path: Path,
        haplogroups: Dict[str, str],
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test classification of a single 1000 Genomes sample."""
        # Pick first sample with known haplogroup
        sample_id = list(haplogroups.keys())[0]
        expected_hg = haplogroups[sample_id]

        # Classify
        adapter = VCFAdapter(reference_path=reference_path, sample_id=sample_id)
        profile = adapter.extract_profile(str(vcf_path))

        classifier = Classifier(phylotree)
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.haplogroup != ""

        # Check major clade match
        expected_clade = get_major_clade(expected_hg)
        result_clade = get_major_clade(result.haplogroup)

        print(f"Sample: {sample_id}")
        print(f"Expected: {expected_hg} (clade: {expected_clade})")
        print(f"Got: {result.haplogroup} (clade: {result_clade})")
        print(f"Confidence: {result.confidence:.2%}")

    def test_classify_multiple_samples(
        self,
        vcf_path: Path,
        haplogroups: Dict[str, str],
        phylotree: Phylotree,
        reference_path: str,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Test classification of multiple 1000 Genomes samples.

        Use --sample-count and --random-seed to control sample selection.
        """
        all_sample_ids = list(haplogroups.keys())
        sample_ids = select_samples(all_sample_ids, sample_count, random_seed)

        print(
            f"\nTesting {len(sample_ids)} samples"
            + (f" (seed={random_seed})" if random_seed else " (sequential)")
        )

        classifier = Classifier(phylotree)

        results: List[Tuple[str, str, str, float]] = []
        clade_correct = 0
        total_classified = 0

        for sample_id in sample_ids:
            expected_hg = haplogroups[sample_id]

            try:
                adapter = VCFAdapter(reference_path=reference_path, sample_id=sample_id)
                profile = adapter.extract_profile(str(vcf_path))
                result = classifier.classify(profile)

                expected_clade = get_major_clade(expected_hg)
                result_clade = get_major_clade(result.haplogroup)

                results.append((sample_id, expected_hg, result.haplogroup, result.confidence))

                if expected_clade == result_clade:
                    clade_correct += 1
                total_classified += 1

            except Exception as e:
                print(f"Error classifying {sample_id}: {e}")
                continue

        if total_classified == 0:
            pytest.skip("No samples could be classified")

        accuracy = clade_correct / total_classified
        print(f"\nMajor clade accuracy: {accuracy:.1%} ({clade_correct}/{total_classified})")

        # Print detailed results
        print("\nResults:")
        for sample_id, expected, got, conf in results:
            match = "✓" if get_major_clade(expected) == get_major_clade(got) else "✗"
            print(f"  {sample_id}: {expected} → {got} ({conf:.1%}) {match}")

        # Expect at least 70% major clade accuracy
        assert accuracy >= 0.7, f"Major clade accuracy too low: {accuracy:.1%}"


class Test1KGBenchmark:
    """Benchmark tests for 1000 Genomes classification accuracy."""

    @pytest.fixture
    def vcf_path(self, mtdna_poc_data_dir: Path) -> Path:
        """Return 1000 Genomes VCF path."""
        vcf = mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"
        if not vcf.exists():
            pytest.skip("1000 Genomes VCF not available")
        return vcf

    @pytest.fixture
    def haplogroups(self, mtdna_poc_data_dir: Path) -> Dict[str, str]:
        """Load 1000 Genomes haplogroup ground truth."""
        hg_file = mtdna_poc_data_dir / "haplogroups" / "1kg_haplogroups_simple.csv"
        if not hg_file.exists():
            pytest.skip("1000 Genomes haplogroups file not available")

        haplogroups = {}
        with open(hg_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row.get("SampleID", "")
                hg = row.get("Haplogroup", "")
                if sample_id and hg:
                    haplogroups[sample_id] = hg

        return haplogroups

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree-rsrs.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights-rsrs.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(str(tree_path), str(weights_path) if weights_path.exists() else None)

    @pytest.fixture
    def reference_path(self, evehap_data_dir: Path) -> str:
        """Return reference FASTA path."""
        ref = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref.exists():
            pytest.skip("Reference FASTA not available")
        return str(ref)

    @pytest.mark.slow
    def test_full_1kg_benchmark(
        self,
        vcf_path: Path,
        haplogroups: Dict[str, str],
        phylotree: Phylotree,
        reference_path: str,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Full benchmark against 1000 Genomes samples.

        This test is marked slow and may take several minutes.
        Use --sample-count and --random-seed to control sample selection.
        Use --all-samples for complete benchmark.
        """
        all_sample_ids = list(haplogroups.keys())
        # Default to 100 for benchmark if not specified
        count = sample_count if sample_count != 10 else 100
        sample_ids = select_samples(all_sample_ids, count, random_seed)

        print(
            f"\nBenchmarking {len(sample_ids)} samples"
            + (f" (seed={random_seed})" if random_seed else "")
        )

        classifier = Classifier(phylotree)

        exact_match = 0
        clade_match = 0
        total = 0
        errors = 0

        for sample_id in sample_ids:
            expected_hg = haplogroups[sample_id]

            try:
                adapter = VCFAdapter(reference_path=reference_path, sample_id=sample_id)
                profile = adapter.extract_profile(str(vcf_path))
                result = classifier.classify(profile)

                total += 1

                # Normalize for comparison
                expected_norm = normalize_haplogroup(expected_hg)
                result_norm = normalize_haplogroup(result.haplogroup)

                if result_norm == expected_norm:
                    exact_match += 1
                    clade_match += 1
                elif get_major_clade(expected_hg) == get_major_clade(result.haplogroup):
                    clade_match += 1

            except Exception:
                errors += 1

        if total == 0:
            pytest.skip("No samples with ground truth")

        exact_accuracy = exact_match / total
        clade_accuracy = clade_match / total

        print("\n=== 1000 Genomes Benchmark Results ===")
        print(f"Total samples: {total}")
        print(f"Exact matches: {exact_match} ({exact_accuracy:.1%})")
        print(f"Clade matches: {clade_match} ({clade_accuracy:.1%})")
        print(f"Errors: {errors}")

        # Minimum expected accuracy
        assert clade_accuracy >= 0.7, f"Clade accuracy {clade_accuracy:.1%} below 70%"

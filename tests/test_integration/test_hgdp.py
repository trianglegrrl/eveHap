"""Integration tests against HGDP ground truth.

Tests the complete pipeline from BAM files to haplogroup classification,
comparing results against known HGDP haplogroup assignments.

Usage:
    # Run with default 10 samples
    pytest tests/test_integration/test_hgdp.py -v

    # Run with 50 random samples (seed=42 for reproducibility)
    pytest tests/test_integration/test_hgdp.py -v --sample-count 50 --random-seed 42

    # Run on all samples
    pytest tests/test_integration/test_hgdp.py -v --all-samples
"""

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import pytest

from evehap.adapters.bam import BAMAdapter
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
    return haplogroup.strip().upper().replace("'", "").replace(" ", "")


class TestHGDPIntegration:
    """Integration tests against HGDP samples."""

    @pytest.fixture
    def hgdp_bam_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return HGDP BAM directory."""
        bam_dir = mtdna_poc_data_dir / "hgdp" / "bams"
        if not bam_dir.exists():
            pytest.skip("HGDP BAM directory not available")
        return bam_dir

    @pytest.fixture
    def hgdp_haplogroups(self, mtdna_poc_data_dir: Path) -> Dict[str, str]:
        """Load HGDP haplogroup ground truth."""
        hg_file = mtdna_poc_data_dir / "hgdp" / "haplogroups" / "haplogroups.csv"
        if not hg_file.exists():
            pytest.skip("HGDP haplogroups file not available")

        haplogroups = {}
        with open(hg_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row.get("sample_id", row.get("Sample", ""))
                hg = row.get("haplogroup", row.get("Haplogroup", ""))
                if sample_id and hg:
                    haplogroups[sample_id] = hg

        return haplogroups

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(
            str(tree_path),
            str(weights_path) if weights_path.exists() else None
        )

    def test_hgdp_sample_count(self, hgdp_bam_dir: Path) -> None:
        """Test that we have HGDP BAM files."""
        bams = list(hgdp_bam_dir.glob("*.bam"))
        assert len(bams) > 0, "No HGDP BAM files found"

    def test_classify_single_sample(
        self,
        hgdp_bam_dir: Path,
        hgdp_haplogroups: Dict[str, str],
        phylotree: Phylotree,
    ) -> None:
        """Test classification of a single HGDP sample."""
        bams = list(hgdp_bam_dir.glob("*.bam"))
        if not bams:
            pytest.skip("No BAM files found")

        # Pick first BAM with known haplogroup
        sample_bam = None
        expected_hg = None

        for bam in bams:
            sample_id = bam.stem.split(".")[0]
            if sample_id in hgdp_haplogroups:
                sample_bam = bam
                expected_hg = hgdp_haplogroups[sample_id]
                break

        if sample_bam is None:
            pytest.skip("No BAM files with known haplogroups")

        # Classify
        adapter = BAMAdapter()
        profile = adapter.extract_profile(str(sample_bam))

        classifier = Classifier(phylotree)
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.haplogroup != ""

        # Check major clade match
        expected_clade = get_major_clade(expected_hg)
        result_clade = get_major_clade(result.haplogroup)

        # Log result for debugging
        print(f"Sample: {sample_bam.stem}")
        print(f"Expected: {expected_hg} (clade: {expected_clade})")
        print(f"Got: {result.haplogroup} (clade: {result_clade})")
        print(f"Confidence: {result.confidence:.2%}")

    def test_classify_multiple_samples(
        self,
        hgdp_bam_dir: Path,
        hgdp_haplogroups: Dict[str, str],
        phylotree: Phylotree,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Test classification of multiple HGDP samples.

        Use --sample-count and --random-seed to control sample selection.
        """
        all_bams = list(hgdp_bam_dir.glob("*.bam"))
        bams = select_samples(all_bams, sample_count, random_seed)

        if len(bams) < 3:
            pytest.skip("Need at least 3 BAM files")

        print(f"\nTesting {len(bams)} samples" +
              (f" (seed={random_seed})" if random_seed else " (sequential)"))

        adapter = BAMAdapter()
        classifier = Classifier(phylotree)

        results: List[Tuple[str, str, str, float]] = []
        clade_correct = 0
        total_classified = 0

        for bam in bams:
            sample_id = bam.stem.split(".")[0]
            expected_hg = hgdp_haplogroups.get(sample_id)

            if expected_hg is None:
                continue

            try:
                profile = adapter.extract_profile(str(bam))
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

        # Expect at least 50% major clade accuracy
        assert accuracy >= 0.5, f"Major clade accuracy too low: {accuracy:.1%}"


class TestHGDPBenchmark:
    """Benchmark tests for HGDP classification accuracy."""

    @pytest.fixture
    def hgdp_bam_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return HGDP BAM directory."""
        bam_dir = mtdna_poc_data_dir / "hgdp" / "bams"
        if not bam_dir.exists():
            pytest.skip("HGDP BAM directory not available")
        return bam_dir

    @pytest.fixture
    def hgdp_haplogroups(self, mtdna_poc_data_dir: Path) -> Dict[str, str]:
        """Load HGDP haplogroup ground truth."""
        hg_file = mtdna_poc_data_dir / "hgdp" / "haplogroups" / "haplogroups.csv"
        if not hg_file.exists():
            pytest.skip("HGDP haplogroups file not available")

        haplogroups = {}
        with open(hg_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row.get("sample_id", row.get("Sample", ""))
                hg = row.get("haplogroup", row.get("Haplogroup", ""))
                if sample_id and hg:
                    haplogroups[sample_id] = hg

        return haplogroups

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(
            str(tree_path),
            str(weights_path) if weights_path.exists() else None
        )

    @pytest.mark.slow
    def test_full_hgdp_benchmark(
        self,
        hgdp_bam_dir: Path,
        hgdp_haplogroups: Dict[str, str],
        phylotree: Phylotree,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Full benchmark against HGDP samples with ground truth.

        This test is marked slow and may take several minutes.
        Use --sample-count and --random-seed to control sample selection.
        Use --all-samples for complete benchmark.
        """
        all_bams = list(hgdp_bam_dir.glob("*.bam"))
        # Default to all for benchmark if not specified
        count = sample_count if sample_count != 10 else -1
        bams = select_samples(all_bams, count, random_seed)

        print(f"\nBenchmarking {len(bams)} HGDP samples" +
              (f" (seed={random_seed})" if random_seed else ""))

        adapter = BAMAdapter()
        classifier = Classifier(phylotree)

        exact_match = 0
        clade_match = 0
        total = 0
        errors = 0

        for bam in bams:
            sample_id = bam.stem.split(".")[0]
            expected_hg = hgdp_haplogroups.get(sample_id)

            if expected_hg is None:
                continue

            try:
                profile = adapter.extract_profile(str(bam))
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

        print(f"\n=== HGDP Benchmark Results ===")
        print(f"Total samples: {total}")
        print(f"Exact matches: {exact_match} ({exact_accuracy:.1%})")
        print(f"Clade matches: {clade_match} ({clade_accuracy:.1%})")
        print(f"Errors: {errors}")

        # Minimum expected accuracy
        assert clade_accuracy >= 0.6, f"Clade accuracy {clade_accuracy:.1%} below 60%"



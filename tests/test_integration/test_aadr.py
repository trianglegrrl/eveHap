"""Integration tests against AADR ancient DNA samples.

Tests the complete pipeline from BAM files to haplogroup classification
with ancient DNA damage filtering.

Usage:
    # Run with default 10 samples
    pytest tests/test_integration/test_aadr.py -v

    # Run with 20 random samples (seed=42 for reproducibility)
    pytest tests/test_integration/test_aadr.py -v --sample-count 20 --random-seed 42

    # Run on all samples
    pytest tests/test_integration/test_aadr.py -v --all-samples
"""

from pathlib import Path
from typing import List, Tuple

import pytest

from evehap.adapters.bam import BAMAdapter
from evehap.core.classifier import Classifier
from evehap.core.damage import DamageFilter
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


class TestAADRIntegration:
    """Integration tests against AADR ancient DNA samples."""

    @pytest.fixture
    def aadr_bam_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return AADR BAM directory."""
        bam_dir = mtdna_poc_data_dir / "aadr_bams"
        if not bam_dir.exists():
            pytest.skip("AADR BAM directory not available")
        return bam_dir

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        # Use RSRS tree for ancient DNA
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

    def test_aadr_bams_exist(self, aadr_bam_dir: Path) -> None:
        """Test that we have AADR BAM files."""
        bams = list(aadr_bam_dir.glob("*.bam"))
        assert len(bams) > 0, "No AADR BAM files found"
        print(f"Found {len(bams)} AADR BAM files")

    def test_classify_single_sample(
        self,
        aadr_bam_dir: Path,
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test classification of a single AADR sample."""
        bams = list(aadr_bam_dir.glob("*.bam"))
        if not bams:
            pytest.skip("No BAM files found")

        # Pick first non-empty BAM
        sample_bam = None
        for bam in bams:
            if bam.stat().st_size > 1000:  # Skip very small/empty files
                sample_bam = bam
                break

        if sample_bam is None:
            pytest.skip("No valid BAM files found")

        # Classify with damage filter
        adapter = BAMAdapter(reference_path=reference_path)
        profile = adapter.extract_profile(str(sample_bam))

        damage_filter = DamageFilter()
        classifier = Classifier(phylotree, method="traversal", damage_filter=damage_filter)
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.haplogroup != ""

        print(f"Sample: {sample_bam.stem}")
        print(f"Haplogroup: {result.haplogroup}")
        print(f"Confidence: {result.confidence:.2%}")
        print(f"Coverage: {result.coverage:.2%}")

    def test_classify_multiple_samples(
        self,
        aadr_bam_dir: Path,
        phylotree: Phylotree,
        reference_path: str,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Test classification of multiple AADR samples.

        Use --sample-count and --random-seed to control sample selection.
        """
        all_bams = [b for b in aadr_bam_dir.glob("*.bam") if b.stat().st_size > 1000]
        bams = select_samples(all_bams, sample_count, random_seed)

        if len(bams) < 3:
            pytest.skip("Need at least 3 valid BAM files")

        print(
            f"\nTesting {len(bams)} samples"
            + (f" (seed={random_seed})" if random_seed else " (sequential)")
        )

        adapter = BAMAdapter(reference_path=reference_path)
        damage_filter = DamageFilter()
        classifier = Classifier(phylotree, method="traversal", damage_filter=damage_filter)

        results: List[Tuple[str, str, float, float]] = []
        errors = 0

        for bam in bams:
            sample_id = bam.stem.split(".")[0]

            try:
                profile = adapter.extract_profile(str(bam))
                result = classifier.classify(profile)
                coverage = getattr(result, "coverage_fraction", 0.0) or 0.0
                results.append((sample_id, result.haplogroup, result.confidence, coverage))

            except Exception as e:
                print(f"Error classifying {sample_id}: {e}")
                errors += 1
                continue

        if not results:
            pytest.skip("No samples could be classified")

        print("\n=== AADR Classification Results ===")
        print(f"Total classified: {len(results)}")
        print(f"Errors: {errors}")

        # Print detailed results
        print("\nResults:")
        for sample_id, hg, conf, cov in results:
            print(f"  {sample_id}: {hg} (conf={conf:.1%}, cov={cov:.1%})")

        # Haplogroup distribution
        clades = [get_major_clade(r[1]) for r in results]
        from collections import Counter

        clade_counts = Counter(clades)
        print("\nClade distribution:")
        for clade, count in clade_counts.most_common(10):
            print(f"  {clade}: {count}")

        # At least 50% should classify successfully
        success_rate = len(results) / (len(results) + errors)
        assert success_rate >= 0.5, f"Success rate too low: {success_rate:.1%}"

    def test_damage_detection(
        self,
        aadr_bam_dir: Path,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Test ancient DNA damage detection on AADR samples."""
        all_bams = [b for b in aadr_bam_dir.glob("*.bam") if b.stat().st_size > 10000]
        bams = select_samples(all_bams, min(sample_count, 5), random_seed)

        if len(bams) < 2:
            pytest.skip("Need at least 2 BAM files for damage analysis")

        print(f"\nAnalyzing damage in {len(bams)} samples")

        damage_filter = DamageFilter()
        damage_results = []

        for bam in bams:
            sample_id = bam.stem.split(".")[0]
            try:
                stats = damage_filter.estimate_damage(str(bam), max_reads=1000)
                damage_results.append(
                    (sample_id, stats.ct_5prime_rate, stats.ga_3prime_rate, stats.is_ancient)
                )
            except Exception as e:
                print(f"Error analyzing {sample_id}: {e}")
                continue

        if not damage_results:
            pytest.skip("No samples could be analyzed for damage")

        print("\nDamage analysis results:")
        ancient_count = 0
        for sample_id, ct_rate, ga_rate, is_ancient in damage_results:
            status = "ANCIENT" if is_ancient else "modern"
            print(f"  {sample_id}: C→T={ct_rate:.1%}, G→A={ga_rate:.1%} [{status}]")
            if is_ancient:
                ancient_count += 1

        print(f"\nAncient DNA detected: {ancient_count}/{len(damage_results)}")


class TestAADRBenchmark:
    """Benchmark tests for AADR classification."""

    @pytest.fixture
    def aadr_bam_dir(self, mtdna_poc_data_dir: Path) -> Path:
        """Return AADR BAM directory."""
        bam_dir = mtdna_poc_data_dir / "aadr_bams"
        if not bam_dir.exists():
            pytest.skip("AADR BAM directory not available")
        return bam_dir

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
    def test_full_aadr_benchmark(
        self,
        aadr_bam_dir: Path,
        phylotree: Phylotree,
        reference_path: str,
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Full benchmark against AADR samples.

        This test is marked slow and may take several minutes.
        Use --sample-count and --random-seed to control sample selection.
        """
        all_bams = [b for b in aadr_bam_dir.glob("*.bam") if b.stat().st_size > 1000]
        # Default to 50 for benchmark if not specified
        count = sample_count if sample_count != 10 else 50
        bams = select_samples(all_bams, count, random_seed)

        print(
            f"\nBenchmarking {len(bams)} AADR samples"
            + (f" (seed={random_seed})" if random_seed else "")
        )

        adapter = BAMAdapter(reference_path=reference_path)
        damage_filter = DamageFilter()
        classifier = Classifier(phylotree, method="traversal", damage_filter=damage_filter)

        results = []
        errors = 0

        for bam in bams:
            sample_id = bam.stem.split(".")[0]
            try:
                profile = adapter.extract_profile(str(bam))
                result = classifier.classify(profile)
                coverage = getattr(result, "coverage_fraction", 0.0) or 0.0
                results.append((sample_id, result.haplogroup, result.confidence, coverage))
            except Exception:
                errors += 1

        if not results:
            pytest.skip("No samples could be classified")

        # Calculate statistics
        avg_confidence = sum(r[2] for r in results) / len(results)
        avg_coverage = sum(r[3] for r in results) / len(results)
        clades = [get_major_clade(r[1]) for r in results]
        from collections import Counter

        clade_counts = Counter(clades)

        print("\n=== AADR Benchmark Results ===")
        print(f"Total classified: {len(results)}")
        print(f"Errors: {errors}")
        print(f"Average confidence: {avg_confidence:.1%}")
        print(f"Average coverage: {avg_coverage:.1%}")
        print("\nTop clades:")
        for clade, count in clade_counts.most_common(10):
            print(f"  {clade}: {count} ({count/len(results):.1%})")

        # Success rate should be reasonable
        success_rate = len(results) / (len(results) + errors)
        assert success_rate >= 0.5, f"Success rate too low: {success_rate:.1%}"

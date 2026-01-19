#!/usr/bin/env python3
"""Run integration tests manually (bypasses pytest plugin conflicts).

Usage:
    python tests/run_integration_tests.py                    # Default: 10 samples each
    python tests/run_integration_tests.py --sample-count 50  # 50 samples each
    python tests/run_integration_tests.py --random-seed 42   # Reproducible random selection
    python tests/run_integration_tests.py --all-samples      # All samples (slow)
    python tests/run_integration_tests.py --dataset 1kg      # Only 1000 Genomes
    python tests/run_integration_tests.py --dataset hgdp     # Only HGDP
    python tests/run_integration_tests.py --dataset aadr     # Only AADR
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from evehap.adapters.bam import BAMAdapter
from evehap.adapters.vcf import VCFAdapter
from evehap.core.classifier import Classifier
from evehap.core.damage import DamageFilter
from evehap.core.phylotree import Phylotree

# Paths
EVEHAP_ROOT = Path(__file__).parent.parent
DATA_DIR = EVEHAP_ROOT.parent / "data"


def select_samples(items: List, count: int, seed: Optional[int] = None) -> List:
    """Select a subset of samples."""
    import random

    if count == -1 or count >= len(items):
        return items

    if seed is not None:
        rng = random.Random(seed)
        return rng.sample(items, min(count, len(items)))
    else:
        return items[:count]


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade from full haplogroup."""
    if not haplogroup:
        return ""
    for i, char in enumerate(haplogroup):
        if char.isdigit():
            return haplogroup[:i] if i > 0 else haplogroup[0]
    return haplogroup[0]


def load_1kg_haplogroups() -> Dict[str, str]:
    """Load 1000 Genomes ground truth."""
    hg_file = DATA_DIR / "haplogroups" / "1kg_haplogroups_simple.csv"
    if not hg_file.exists():
        return {}

    haplogroups = {}
    with open(hg_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            haplogroups[row["SampleID"]] = row["Haplogroup"]
    return haplogroups


def load_hgdp_haplogroups() -> Dict[str, str]:
    """Load HGDP ground truth."""
    hg_file = DATA_DIR / "hgdp" / "haplogroups" / "haplogroups.csv"
    if not hg_file.exists():
        return {}

    haplogroups = {}
    with open(hg_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row.get("sample_id", row.get("Sample", ""))
            hg = row.get("haplogroup", row.get("Haplogroup", ""))
            if sample_id and hg:
                haplogroups[sample_id] = hg
    return haplogroups


def test_1kg(
    phylotree: Phylotree,
    reference_path: str,
    sample_count: int,
    random_seed: Optional[int],
) -> Tuple[int, int, int]:
    """Test 1000 Genomes VCF classification."""
    print("\n" + "=" * 60)
    print("1000 GENOMES VCF INTEGRATION TEST")
    print("=" * 60)

    vcf_path = DATA_DIR / "raw" / "1kg_chrMT.vcf.gz"
    if not vcf_path.exists():
        print("SKIP: 1000 Genomes VCF not found")
        return (0, 0, 0)

    haplogroups = load_1kg_haplogroups()
    if not haplogroups:
        print("SKIP: Ground truth not found")
        return (0, 0, 0)

    sample_ids = select_samples(list(haplogroups.keys()), sample_count, random_seed)
    print(
        f"Testing {len(sample_ids)} samples"
        + (f" (seed={random_seed})" if random_seed else " (sequential)")
    )

    classifier = Classifier(phylotree)

    correct = 0
    errors = 0
    results = []

    for sid in sample_ids:
        try:
            adapter = VCFAdapter(reference_path=reference_path, sample_id=sid)
            profile = adapter.extract_profile(str(vcf_path))
            result = classifier.classify(profile)

            expected = haplogroups[sid]
            e_clade = get_major_clade(expected)
            r_clade = get_major_clade(result.haplogroup)

            if e_clade == r_clade:
                correct += 1
                status = "✓"
            else:
                status = "✗"

            results.append((sid, expected, result.haplogroup, result.confidence, status))
        except Exception as e:
            errors += 1
            print(f"  ERROR {sid}: {e}")

    # Print results
    print("\nResults:")
    for sid, expected, got, conf, status in results:
        print(f"  {sid}: {expected} → {got} (conf={conf:.1%}) {status}")

    total = len(results)
    accuracy = correct / total if total > 0 else 0
    print(f"\nClade accuracy: {correct}/{total} ({accuracy:.1%})")
    print(f"Errors: {errors}")

    return (correct, total, errors)


def test_hgdp(
    phylotree: Phylotree,
    reference_path: str,
    sample_count: int,
    random_seed: Optional[int],
) -> Tuple[int, int, int]:
    """Test HGDP BAM classification."""
    print("\n" + "=" * 60)
    print("HGDP BAM INTEGRATION TEST")
    print("=" * 60)

    bam_dir = DATA_DIR / "hgdp" / "bams"
    if not bam_dir.exists():
        print("SKIP: HGDP BAM directory not found")
        return (0, 0, 0)

    haplogroups = load_hgdp_haplogroups()
    if not haplogroups:
        print("SKIP: Ground truth not found")
        return (0, 0, 0)

    bams = list(bam_dir.glob("*.bam"))
    bams = select_samples(bams, sample_count, random_seed)
    print(
        f"Testing {len(bams)} samples"
        + (f" (seed={random_seed})" if random_seed else " (sequential)")
    )

    adapter = BAMAdapter(reference_path=reference_path)
    classifier = Classifier(phylotree)

    correct = 0
    errors = 0
    results = []

    for bam in bams:
        sample_id = bam.stem.split(".")[0]
        expected = haplogroups.get(sample_id)

        if expected is None:
            continue

        try:
            profile = adapter.extract_profile(str(bam))
            result = classifier.classify(profile)

            e_clade = get_major_clade(expected)
            r_clade = get_major_clade(result.haplogroup)

            if e_clade == r_clade:
                correct += 1
                status = "✓"
            else:
                status = "✗"

            results.append((sample_id, expected, result.haplogroup, result.confidence, status))
        except Exception as e:
            errors += 1
            print(f"  ERROR {sample_id}: {e}")

    # Print results
    print("\nResults:")
    for sid, expected, got, conf, status in results:
        print(f"  {sid}: {expected} → {got} (conf={conf:.1%}) {status}")

    total = len(results)
    accuracy = correct / total if total > 0 else 0
    print(f"\nClade accuracy: {correct}/{total} ({accuracy:.1%})")
    print(f"Errors: {errors}")

    return (correct, total, errors)


def test_aadr(
    phylotree: Phylotree,
    reference_path: str,
    sample_count: int,
    random_seed: Optional[int],
) -> Tuple[int, int, int]:
    """Test AADR ancient DNA BAM classification."""
    print("\n" + "=" * 60)
    print("AADR ANCIENT DNA INTEGRATION TEST")
    print("=" * 60)

    bam_dir = DATA_DIR / "aadr_bams"
    if not bam_dir.exists():
        print("SKIP: AADR BAM directory not found")
        return (0, 0, 0)

    # Filter to non-empty BAMs
    bams = [b for b in bam_dir.glob("*.bam") if b.stat().st_size > 1000]
    bams = select_samples(bams, sample_count, random_seed)
    print(
        f"Testing {len(bams)} samples"
        + (f" (seed={random_seed})" if random_seed else " (sequential)")
    )

    adapter = BAMAdapter(reference_path=reference_path)
    damage_filter = DamageFilter()
    classifier = Classifier(phylotree, method="traversal", damage_filter=damage_filter)

    errors = 0
    results = []

    for bam in bams:
        sample_id = bam.stem.split(".")[0]

        try:
            profile = adapter.extract_profile(str(bam))
            result = classifier.classify(profile)
            coverage = getattr(result, "coverage_fraction", 0.0) or 0.0
            results.append((sample_id, result.haplogroup, result.confidence, coverage))
        except Exception as e:
            errors += 1
            print(f"  ERROR {sample_id}: {e}")

    # Print results
    print("\nResults:")
    for sid, hg, conf, cov in results:
        print(f"  {sid}: {hg} (conf={conf:.1%}, cov={cov:.1%})")

    # Haplogroup distribution
    clades = [get_major_clade(r[1]) for r in results]
    clade_counts = Counter(clades)
    print("\nClade distribution:")
    for clade, count in clade_counts.most_common(10):
        print(f"  {clade}: {count}")

    print(f"\nClassified: {len(results)}")
    print(f"Errors: {errors}")

    # No ground truth for AADR, just return success count
    return (len(results), len(results), errors)


def main():
    parser = argparse.ArgumentParser(description="Run eveHap integration tests")
    parser.add_argument(
        "--sample-count", type=int, default=10, help="Number of samples to test (default: 10)"
    )
    parser.add_argument(
        "--random-seed", type=int, default=None, help="Random seed for sample selection"
    )
    parser.add_argument("--all-samples", action="store_true", help="Test all samples (slow)")
    parser.add_argument(
        "--dataset", choices=["all", "1kg", "hgdp", "aadr"], default="all", help="Dataset to test"
    )
    args = parser.parse_args()

    sample_count = -1 if args.all_samples else args.sample_count

    # Load resources
    tree_path = EVEHAP_ROOT / "data" / "phylotree" / "tree-rsrs.xml"
    ref_path = str(EVEHAP_ROOT / "data" / "reference" / "rCRS.fasta")

    if not tree_path.exists():
        print(f"ERROR: Phylotree not found: {tree_path}")
        sys.exit(1)

    print("Loading phylotree...")
    phylotree = Phylotree.load(str(tree_path), None)
    print(f"Loaded {len(phylotree.nodes)} haplogroups")

    # Run tests
    total_correct = 0
    total_samples = 0
    total_errors = 0

    if args.dataset in ["all", "1kg"]:
        c, t, e = test_1kg(phylotree, ref_path, sample_count, args.random_seed)
        total_correct += c
        total_samples += t
        total_errors += e

    if args.dataset in ["all", "hgdp"]:
        c, t, e = test_hgdp(phylotree, ref_path, sample_count, args.random_seed)
        total_correct += c
        total_samples += t
        total_errors += e

    if args.dataset in ["all", "aadr"]:
        c, t, e = test_aadr(phylotree, ref_path, sample_count, args.random_seed)
        total_correct += c
        total_samples += t
        total_errors += e

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total classified: {total_samples}")
    print(f"Total errors: {total_errors}")
    if total_samples > 0:
        accuracy = total_correct / total_samples
        print(f"Overall clade accuracy: {total_correct}/{total_samples} ({accuracy:.1%})")

    # Exit code
    if total_errors > total_samples * 0.5:
        print("\nFAILED: Too many errors")
        sys.exit(1)
    else:
        print("\nPASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()

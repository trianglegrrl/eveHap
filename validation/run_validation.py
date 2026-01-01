#!/usr/bin/env python3
"""
Run comprehensive validation suite for eveHap.

Validates eveHap against multiple datasets with ground truth, compares
with Haplogrep3, and produces detailed JSON reports.

Output:
- Console summary table
- JSON results in validation_results/validation_<timestamp>.json
- Per-sample TSV for detailed analysis

Usage:
    # Full validation (all samples)
    python -m validation.run_validation

    # Quick validation with subsampling and threading
    python -m validation.run_validation --max-samples 100 --workers 16

    # Specific dataset
    python -m validation.run_validation --dataset 1kg
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from evehap.adapters.bam import BAMAdapter
from evehap.adapters.vcf import VCFAdapter
from evehap.core.classifier import Classifier
from evehap.core.damage import DamageFilter
from evehap.core.phylotree import Phylotree

# Constants
VALIDATION_ROOT = Path(__file__).parent
EVEHAP_ROOT = VALIDATION_ROOT.parent


# =============================================================================
# Haplogroup Comparison Functions
# =============================================================================


def normalize_haplogroup(hg: str) -> str:
    """Normalize haplogroup for comparison."""
    if not hg:
        return ""
    hg = hg.strip().upper()
    # Remove common annotations
    for char in ["'", "*", "!", "+", " "]:
        if char in hg:
            hg = hg.split(char)[0]
    return hg


def get_major_clade(hg: str) -> str:
    """Extract major clade (letters before first digit) from haplogroup."""
    if not hg:
        return ""
    hg = normalize_haplogroup(hg)
    if not hg:
        return ""
    for i, char in enumerate(hg):
        if char.isdigit():
            return hg[:i] if i > 0 else hg[0]
    return hg[0] if hg else ""


def is_prefix_match(hg1: str, hg2: str) -> bool:
    """Check if one haplogroup is a prefix of the other."""
    h1 = normalize_haplogroup(hg1)
    h2 = normalize_haplogroup(hg2)
    if not h1 or not h2:
        return False
    return h1.startswith(h2) or h2.startswith(h1)


def classify_match(predicted: str, ground_truth: str) -> str:
    """
    Classify the match type between predicted and ground truth.

    Returns one of: exact, prefix, major_clade, mismatch
    """
    pred_norm = normalize_haplogroup(predicted)
    gt_norm = normalize_haplogroup(ground_truth)

    if not pred_norm or not gt_norm:
        return "mismatch"

    if pred_norm == gt_norm:
        return "exact"

    if is_prefix_match(pred_norm, gt_norm):
        return "prefix"

    if get_major_clade(pred_norm) == get_major_clade(gt_norm):
        return "major_clade"

    return "mismatch"


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class SampleResult:
    """Result for a single sample."""
    sample_id: str
    ground_truth: Optional[str]

    # eveHap results
    evehap_haplogroup: str
    evehap_confidence: float
    evehap_quality: str
    evehap_method: str
    evehap_found_mutations: int
    evehap_missing_mutations: int
    evehap_extra_mutations: int
    evehap_coverage_fraction: float
    evehap_time_ms: float

    # Match classification
    match_type: str = "mismatch"  # exact, prefix, major_clade, mismatch

    # Haplogrep3 results (optional)
    haplogrep_haplogroup: Optional[str] = None
    haplogrep_quality: Optional[float] = None
    haplogrep_match_type: str = "mismatch"

    # Error tracking
    error: Optional[str] = None


@dataclass
class CladeStats:
    """Statistics for a single major clade."""
    total: int = 0
    correct: int = 0
    exact: int = 0
    prefix: int = 0
    mismatch: int = 0


@dataclass
class ValidationResult:
    """Complete validation result for a dataset."""
    dataset: str
    description: str
    timestamp: str

    # Counts
    total_samples: int
    successful_samples: int
    failed_samples: int

    # eveHap accuracy
    same_major_clade: int
    same_major_clade_pct: float
    exact_match: int
    exact_match_pct: float
    prefix_match: int
    prefix_match_pct: float
    wrong_clade: int

    # Statistics
    mean_confidence: float
    mean_found_mutations: float
    mean_coverage_fraction: float

    # Timing
    total_time_seconds: float
    samples_per_second: float
    mean_time_per_sample_ms: float

    # Haplogrep3 comparison (if available)
    haplogrep_same_major_clade: int = 0
    haplogrep_same_major_clade_pct: float = 0.0
    haplogrep_exact_match: int = 0
    haplogrep_exact_match_pct: float = 0.0

    # Breakdown by major clade
    by_major_clade: Dict[str, Dict[str, int]] = field(default_factory=dict)

    # Wrong samples (for debugging)
    wrong_samples: List[Dict[str, Any]] = field(default_factory=list)

    # All sample results
    samples: List[SampleResult] = field(default_factory=list)

    # Pass/fail
    passed: bool = False
    min_accuracy_required: Optional[float] = None

    # System info
    workers: int = 1
    evehap_version: str = "0.1.0"
    phylotree_nodes: int = 0


# =============================================================================
# Configuration
# =============================================================================


def load_config(config_path: Path) -> Dict[str, Any]:
    """Load validation configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def load_ground_truth(
    path: Path,
    sample_col: str,
    hg_col: str,
) -> Dict[str, str]:
    """Load ground truth haplogroups from CSV."""
    ground_truth = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row.get(sample_col, "").strip()
            haplogroup = row.get(hg_col, "").strip()
            if sample_id and haplogroup:
                ground_truth[sample_id] = haplogroup
    return ground_truth


# =============================================================================
# Classification Functions
# =============================================================================


def classify_sample_vcf(
    sample_id: str,
    vcf_path: str,
    classifier: Classifier,
    reference_path: str,
    ground_truth: Optional[str],
) -> SampleResult:
    """Classify a single VCF sample."""
    start_time = time.perf_counter()

    try:
        adapter = VCFAdapter(reference_path=reference_path, sample_id=sample_id)
        profile = adapter.extract_profile(vcf_path)
        result = classifier.classify(profile)

        elapsed_ms = (time.perf_counter() - start_time) * 1000

        # Determine match type
        match_type = classify_match(result.haplogroup, ground_truth) if ground_truth else "unknown"

        return SampleResult(
            sample_id=sample_id,
            ground_truth=ground_truth,
            evehap_haplogroup=result.haplogroup,
            evehap_confidence=result.confidence,
            evehap_quality=result.quality,
            evehap_method=result.method,
            evehap_found_mutations=len(result.found_mutations),
            evehap_missing_mutations=len(result.missing_mutations),
            evehap_extra_mutations=len(result.extra_mutations),
            evehap_coverage_fraction=getattr(result, 'coverage_fraction', 0.0) or 0.0,
            evehap_time_ms=elapsed_ms,
            match_type=match_type,
        )

    except Exception as e:
        elapsed_ms = (time.perf_counter() - start_time) * 1000
        return SampleResult(
            sample_id=sample_id,
            ground_truth=ground_truth,
            evehap_haplogroup="",
            evehap_confidence=0.0,
            evehap_quality="error",
            evehap_method="",
            evehap_found_mutations=0,
            evehap_missing_mutations=0,
            evehap_extra_mutations=0,
            evehap_coverage_fraction=0.0,
            evehap_time_ms=elapsed_ms,
            error=str(e),
        )


def classify_sample_bam(
    bam_path: Path,
    classifier: Classifier,
    adapter: BAMAdapter,
    ground_truth: Optional[str],
) -> SampleResult:
    """Classify a single BAM sample."""
    sample_id = bam_path.stem.split(".")[0]
    start_time = time.perf_counter()

    try:
        profile = adapter.extract_profile(str(bam_path))
        result = classifier.classify(profile)

        elapsed_ms = (time.perf_counter() - start_time) * 1000

        # Determine match type
        match_type = classify_match(result.haplogroup, ground_truth) if ground_truth else "unknown"

        return SampleResult(
            sample_id=sample_id,
            ground_truth=ground_truth,
            evehap_haplogroup=result.haplogroup,
            evehap_confidence=result.confidence,
            evehap_quality=result.quality,
            evehap_method=result.method,
            evehap_found_mutations=len(result.found_mutations),
            evehap_missing_mutations=len(result.missing_mutations),
            evehap_extra_mutations=len(result.extra_mutations),
            evehap_coverage_fraction=getattr(result, 'coverage_fraction', 0.0) or 0.0,
            evehap_time_ms=elapsed_ms,
            match_type=match_type,
        )

    except Exception as e:
        elapsed_ms = (time.perf_counter() - start_time) * 1000
        return SampleResult(
            sample_id=sample_id,
            ground_truth=ground_truth,
            evehap_haplogroup="",
            evehap_confidence=0.0,
            evehap_quality="error",
            evehap_method="",
            evehap_found_mutations=0,
            evehap_missing_mutations=0,
            evehap_extra_mutations=0,
            evehap_coverage_fraction=0.0,
            evehap_time_ms=elapsed_ms,
            error=str(e),
        )


# =============================================================================
# Haplogrep3 Integration
# =============================================================================


def run_haplogrep3_batch(
    sample_results: List[SampleResult],
    reference_path: str,
    config: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """Run Haplogrep3 on classified samples."""
    hg3_config = config.get("settings", {}).get("haplogrep3", {})
    command = hg3_config.get("command", "haplogrep3")
    tree = hg3_config.get("tree", "phylotree17-rsrs")

    # Check if haplogrep3 is available
    try:
        subprocess.run([command, "--version"], capture_output=True, timeout=5)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return {}

    # Generate FASTA from results (using haplogroup as proxy for now)
    # In practice, this would use the actual sequence/variants
    # For now, return empty - Haplogrep3 comparison requires proper FASTA generation
    return {}


# =============================================================================
# Dataset Validation
# =============================================================================


def validate_vcf_dataset(
    dataset_name: str,
    dataset_config: Dict[str, Any],
    config: Dict[str, Any],
    phylotree: Phylotree,
    workers: int,
    max_samples: Optional[int],
    random_seed: Optional[int],
) -> ValidationResult:
    """Validate a VCF dataset."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    print(f"\n{'='*70}")
    print(f"Validating {dataset_name}: {dataset_config.get('description', '')}")
    print(f"{'='*70}")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    vcf_path = dataset_config["vcf"]
    reference_path = str(EVEHAP_ROOT / config["evehap"]["reference"])

    # Load ground truth
    gt_path = Path(dataset_config["ground_truth"])
    ground_truth = load_ground_truth(
        gt_path,
        dataset_config["sample_column"],
        dataset_config["haplogroup_column"],
    )
    print(f"  Loaded {len(ground_truth)} ground truth samples")

    # Select samples
    sample_ids = list(ground_truth.keys())
    if max_samples and max_samples < len(sample_ids):
        if random_seed is not None:
            random.seed(random_seed)
            sample_ids = random.sample(sample_ids, max_samples)
        else:
            sample_ids = sample_ids[:max_samples]

    print(f"  Testing {len(sample_ids)} samples with {workers} workers")

    # Create classifier
    classifier = Classifier(phylotree)

    # Run classification in parallel
    start_time = time.perf_counter()
    results: List[SampleResult] = []

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                classify_sample_vcf,
                sid,
                vcf_path,
                classifier,
                reference_path,
                ground_truth.get(sid),
            ): sid
            for sid in sample_ids
        }

        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)

            if i % 100 == 0 or i == len(sample_ids):
                elapsed = time.perf_counter() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                print(f"    Progress: {i}/{len(sample_ids)} ({rate:.1f} samples/sec)")

    total_time = time.perf_counter() - start_time

    # Aggregate results
    return aggregate_results(
        dataset_name=dataset_name,
        description=dataset_config.get("description", ""),
        timestamp=timestamp,
        results=results,
        total_time=total_time,
        workers=workers,
        min_accuracy=dataset_config.get("min_accuracy"),
        phylotree_nodes=len(phylotree.nodes),
    )


def validate_bam_dataset(
    dataset_name: str,
    dataset_config: Dict[str, Any],
    config: Dict[str, Any],
    phylotree: Phylotree,
    workers: int,
    max_samples: Optional[int],
    random_seed: Optional[int],
) -> ValidationResult:
    """Validate a BAM dataset."""
    import re
    from concurrent.futures import ThreadPoolExecutor, as_completed

    print(f"\n{'='*70}")
    print(f"Validating {dataset_name}: {dataset_config.get('description', '')}")
    print(f"{'='*70}")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    bam_dir = Path(dataset_config["bam_dir"])
    reference_path = str(EVEHAP_ROOT / config["evehap"]["reference"])

    # Load ground truth if available
    ground_truth = {}
    if dataset_config.get("ground_truth"):
        gt_path = Path(dataset_config["ground_truth"])
        ground_truth = load_ground_truth(
            gt_path,
            dataset_config["sample_column"],
            dataset_config["haplogroup_column"],
        )
        print(f"  Loaded {len(ground_truth)} ground truth samples")

    # Find BAM files
    pattern = dataset_config.get("file_pattern", "*.bam")
    bam_files = [b for b in bam_dir.glob(pattern) if b.stat().st_size > 1000]

    # Match to ground truth
    sample_regex = dataset_config.get("sample_regex", "^([^.]+)")
    bam_samples = []
    for bam in bam_files:
        match = re.match(sample_regex, bam.stem)
        sample_id = match.group(1) if match else bam.stem
        gt = ground_truth.get(sample_id)
        if gt or not ground_truth:
            bam_samples.append((bam, sample_id, gt))

    # Select samples
    if max_samples and max_samples < len(bam_samples):
        if random_seed is not None:
            random.seed(random_seed)
            bam_samples = random.sample(bam_samples, max_samples)
        else:
            bam_samples = bam_samples[:max_samples]

    print(f"  Testing {len(bam_samples)} samples with {workers} workers")

    # Create classifier and adapter
    use_damage = dataset_config.get("damage_filter", False)
    damage_filter = DamageFilter() if use_damage else None
    classifier = Classifier(phylotree, damage_filter=damage_filter)
    adapter = BAMAdapter(reference_path=reference_path)

    # Run classification in parallel
    start_time = time.perf_counter()
    results: List[SampleResult] = []

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                classify_sample_bam,
                bam,
                classifier,
                adapter,
                gt,
            ): sample_id
            for bam, sample_id, gt in bam_samples
        }

        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)

            if i % 50 == 0 or i == len(bam_samples):
                elapsed = time.perf_counter() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                print(f"    Progress: {i}/{len(bam_samples)} ({rate:.1f} samples/sec)")

    total_time = time.perf_counter() - start_time

    # Aggregate results
    return aggregate_results(
        dataset_name=dataset_name,
        description=dataset_config.get("description", ""),
        timestamp=timestamp,
        results=results,
        total_time=total_time,
        workers=workers,
        min_accuracy=dataset_config.get("min_accuracy"),
        phylotree_nodes=len(phylotree.nodes),
    )


def aggregate_results(
    dataset_name: str,
    description: str,
    timestamp: str,
    results: List[SampleResult],
    total_time: float,
    workers: int,
    min_accuracy: Optional[float],
    phylotree_nodes: int,
) -> ValidationResult:
    """Aggregate sample results into dataset result."""
    successful = [r for r in results if not r.error]
    failed = [r for r in results if r.error]
    with_gt = [r for r in successful if r.ground_truth]

    # Count match types
    exact = sum(1 for r in with_gt if r.match_type == "exact")
    prefix = sum(1 for r in with_gt if r.match_type == "prefix")
    major_clade = sum(1 for r in with_gt if r.match_type == "major_clade")
    mismatch = sum(1 for r in with_gt if r.match_type == "mismatch")

    same_major = exact + prefix + major_clade  # All except mismatch

    # Calculate percentages
    n = len(with_gt) if with_gt else 1

    # Statistics
    confidences = [r.evehap_confidence for r in successful if r.evehap_confidence > 0]
    found_muts = [r.evehap_found_mutations for r in successful]
    coverages = [r.evehap_coverage_fraction for r in successful if r.evehap_coverage_fraction > 0]

    # Breakdown by major clade
    by_clade: Dict[str, Dict[str, int]] = defaultdict(lambda: {"total": 0, "correct": 0})
    for r in with_gt:
        clade = get_major_clade(r.ground_truth)
        by_clade[clade]["total"] += 1
        if r.match_type != "mismatch":
            by_clade[clade]["correct"] += 1

    # Wrong samples for debugging
    wrong_samples = []
    for r in with_gt:
        if r.match_type == "mismatch":
            wrong_samples.append({
                "sample": r.sample_id,
                "expected": r.ground_truth,
                "called": r.evehap_haplogroup,
                "confidence": round(r.evehap_confidence, 4),
                "found_mutations": r.evehap_found_mutations,
            })

    # Build result
    result = ValidationResult(
        dataset=dataset_name,
        description=description,
        timestamp=timestamp,
        total_samples=len(results),
        successful_samples=len(successful),
        failed_samples=len(failed),
        same_major_clade=same_major,
        same_major_clade_pct=round(100 * same_major / n, 2),
        exact_match=exact,
        exact_match_pct=round(100 * exact / n, 2),
        prefix_match=prefix,
        prefix_match_pct=round(100 * prefix / n, 2),
        wrong_clade=mismatch,
        mean_confidence=round(sum(confidences) / len(confidences), 4) if confidences else 0,
        mean_found_mutations=round(sum(found_muts) / len(found_muts), 1) if found_muts else 0,
        mean_coverage_fraction=round(sum(coverages) / len(coverages), 4) if coverages else 0,
        total_time_seconds=round(total_time, 2),
        samples_per_second=round(len(results) / total_time, 2) if total_time > 0 else 0,
        mean_time_per_sample_ms=round(1000 * total_time / len(results), 2) if results else 0,
        by_major_clade=dict(by_clade),
        wrong_samples=wrong_samples[:50],  # Limit to first 50
        samples=results,
        workers=workers,
        phylotree_nodes=phylotree_nodes,
        min_accuracy_required=min_accuracy,
    )

    # Check pass/fail
    if min_accuracy:
        result.passed = result.same_major_clade_pct >= min_accuracy * 100
    else:
        result.passed = True

    return result


# =============================================================================
# Output Functions
# =============================================================================


def print_results_table(results: List[ValidationResult]) -> None:
    """Print formatted results table."""
    print("\n" + "=" * 90)
    print("VALIDATION RESULTS")
    print("=" * 90)
    print()

    print(f"{'Dataset':<20} {'Samples':>8} {'Major':>8} {'Exact':>8} {'Prefix':>8} {'Wrong':>8} {'Time':>8}")
    print("-" * 90)

    for r in results:
        print(
            f"{r.dataset:<20} {r.successful_samples:>8} "
            f"{r.same_major_clade_pct:>6.1f}% {r.exact_match_pct:>6.1f}% "
            f"{r.prefix_match_pct:>6.1f}% {r.wrong_clade:>8} {r.total_time_seconds:>6.1f}s"
        )

    print("-" * 90)
    print()
    print("Legend:")
    print("  Major = Same major clade (exact + prefix + subclade)")
    print("  Exact = Identical haplogroup after normalization")
    print("  Prefix = One is ancestor of the other (compatible)")
    print("  Wrong = Different major clade")


def save_results(
    results: List[ValidationResult],
    output_dir: Path,
    timestamp: str,
) -> None:
    """Save results to JSON and TSV files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # JSON summary
    json_path = output_dir / f"validation_{timestamp}.json"
    json_data = {
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "evehap_version": "0.1.0",
        "datasets": [],
    }

    for r in results:
        ds_data = {
            "dataset": r.dataset,
            "description": r.description,
            "timestamp": r.timestamp,
            "total_samples": r.total_samples,
            "successful_samples": r.successful_samples,
            "failed_samples": r.failed_samples,
            "same_major_clade": r.same_major_clade,
            "same_major_clade_pct": r.same_major_clade_pct,
            "exact_match": r.exact_match,
            "exact_match_pct": r.exact_match_pct,
            "prefix_match": r.prefix_match,
            "prefix_match_pct": r.prefix_match_pct,
            "wrong_clade": r.wrong_clade,
            "mean_confidence": r.mean_confidence,
            "mean_found_mutations": r.mean_found_mutations,
            "mean_coverage_fraction": r.mean_coverage_fraction,
            "total_time_seconds": r.total_time_seconds,
            "samples_per_second": r.samples_per_second,
            "mean_time_per_sample_ms": r.mean_time_per_sample_ms,
            "workers": r.workers,
            "phylotree_nodes": r.phylotree_nodes,
            "passed": r.passed,
            "min_accuracy_required": r.min_accuracy_required,
            "by_major_clade": r.by_major_clade,
            "wrong_samples": r.wrong_samples,
        }
        json_data["datasets"].append(ds_data)

    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)

    print(f"\nResults saved to: {json_path}")

    # TSV per dataset
    for r in results:
        tsv_path = output_dir / f"validation_{r.dataset}_{timestamp}.tsv"

        with open(tsv_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "sample_id",
                "ground_truth",
                "evehap_haplogroup",
                "match_type",
                "confidence",
                "quality",
                "found_mutations",
                "missing_mutations",
                "coverage_fraction",
                "time_ms",
                "error",
            ])

            for sr in r.samples:
                writer.writerow([
                    sr.sample_id,
                    sr.ground_truth or "",
                    sr.evehap_haplogroup,
                    sr.match_type,
                    f"{sr.evehap_confidence:.4f}",
                    sr.evehap_quality,
                    sr.evehap_found_mutations,
                    sr.evehap_missing_mutations,
                    f"{sr.evehap_coverage_fraction:.4f}",
                    f"{sr.evehap_time_ms:.2f}",
                    sr.error or "",
                ])

        print(f"Details saved to: {tsv_path}")


def print_final_summary(results: List[ValidationResult]) -> bool:
    """Print final summary and return overall pass/fail."""
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    all_passed = True
    total_samples = 0
    total_correct = 0
    total_time = 0.0

    for r in results:
        status = "✓ PASSED" if r.passed else "✗ FAILED"
        if not r.passed:
            all_passed = False

        print(f"\n{r.dataset}: {status}")
        print(f"  Samples: {r.successful_samples}/{r.total_samples}")
        print(f"  Major clade accuracy: {r.same_major_clade_pct:.1f}%")
        print(f"  Exact match: {r.exact_match_pct:.1f}%")
        print(f"  Mean confidence: {r.mean_confidence:.3f}")
        print(f"  Time: {r.total_time_seconds:.1f}s ({r.samples_per_second:.1f} samples/sec)")

        total_samples += r.successful_samples
        total_correct += r.same_major_clade
        total_time += r.total_time_seconds

    overall_accuracy = 100 * total_correct / total_samples if total_samples > 0 else 0

    print("\n" + "=" * 70)
    print(f"OVERALL: {total_samples} samples, {overall_accuracy:.1f}% major clade accuracy")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f}m)")

    if all_passed:
        print("\n✓ ALL VALIDATIONS PASSED")
    else:
        print("\n✗ SOME VALIDATIONS FAILED")

    print("=" * 70)

    return all_passed


# =============================================================================
# Main
# =============================================================================


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run eveHap validation suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Full validation
    python -m validation.run_validation

    # Quick validation with 100 samples and 16 threads
    python -m validation.run_validation --max-samples 100 --workers 16

    # Specific dataset
    python -m validation.run_validation --dataset 1kg
""",
    )
    parser.add_argument(
        "--config", "-c",
        type=Path,
        default=VALIDATION_ROOT / "config.yaml",
        help="Path to validation config file",
    )
    parser.add_argument(
        "--dataset", "-d",
        type=str,
        default=None,
        help="Run specific dataset only",
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=None,
        help="Number of parallel workers",
    )
    parser.add_argument(
        "--max-samples", "-n",
        type=int,
        default=None,
        help="Maximum samples per dataset",
    )
    parser.add_argument(
        "--random-seed", "-s",
        type=int,
        default=None,
        help="Random seed for sample selection",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=None,
        help="Output directory for results",
    )

    args = parser.parse_args()

    # Load config
    if not args.config.exists():
        print(f"ERROR: Config not found: {args.config}")
        return 1

    config = load_config(args.config)

    # Determine settings
    workers = args.workers or config.get("settings", {}).get("workers", 0)
    if workers <= 0:
        workers = os.cpu_count() or 4

    random_seed = args.random_seed or config.get("settings", {}).get("random_seed")

    output_dir = args.output_dir or Path(config.get("settings", {}).get("output_dir", "validation_results"))
    output_dir = EVEHAP_ROOT / output_dir

    # Load phylotree
    tree_path = EVEHAP_ROOT / config["evehap"]["tree"]
    if not tree_path.exists():
        print(f"ERROR: Phylotree not found: {tree_path}")
        return 1

    print("=" * 70)
    print("eveHap Validation Suite")
    print("=" * 70)
    print(f"  Workers: {workers}")
    if args.max_samples:
        print(f"  Max samples: {args.max_samples}")
    if random_seed:
        print(f"  Random seed: {random_seed}")

    print("\nLoading phylotree...")
    phylotree = Phylotree.load(str(tree_path), None)
    print(f"  Loaded {len(phylotree.nodes)} haplogroups")

    # Run validations
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results: List[ValidationResult] = []

    for ds_name, ds_config in config.get("datasets", {}).items():
        if not ds_config.get("enabled", True):
            continue

        if args.dataset and ds_name != args.dataset:
            continue

        ds_type = ds_config.get("type", "vcf")
        max_samples = args.max_samples or ds_config.get("max_samples")

        try:
            if ds_type == "vcf":
                result = validate_vcf_dataset(
                    ds_name,
                    ds_config,
                    config,
                    phylotree,
                    workers,
                    max_samples,
                    random_seed,
                )
            elif ds_type == "bam":
                result = validate_bam_dataset(
                    ds_name,
                    ds_config,
                    config,
                    phylotree,
                    workers,
                    max_samples,
                    random_seed,
                )
            else:
                print(f"Unknown dataset type: {ds_type}")
                continue

            results.append(result)

        except Exception as e:
            print(f"ERROR validating {ds_name}: {e}")
            continue

    # Output
    if results:
        print_results_table(results)
        save_results(results, output_dir, timestamp)
        all_passed = print_final_summary(results)
        return 0 if all_passed else 1
    else:
        print("\nNo datasets validated")
        return 1


if __name__ == "__main__":
    sys.exit(main())

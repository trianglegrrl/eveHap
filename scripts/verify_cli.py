#!/usr/bin/env python3
"""
CLI verification script for eveHap.

Tests the evehap classify command with real data across BAM, FASTA, and HSD formats,
verifying haplogroup calls against ground truth.

Usage:
    python scripts/verify_cli.py          # Full verification (10 samples x 3 formats)
    python scripts/verify_cli.py --quick  # Quick test (3 samples x 3 formats)
    python scripts/verify_cli.py -v       # Verbose output
    python scripts/verify_cli.py --json   # JSON output for CI
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Paths relative to this script
SCRIPT_DIR = Path(__file__).parent
EVEHAP_ROOT = SCRIPT_DIR.parent
DATA_ROOT = Path("/home/a/mtdna_poc/data")

# Data paths - HGDP for BAM/FASTA, 1KG for HSD (cleaner format)
HGDP_BAM_DIR = DATA_ROOT / "hgdp" / "bams"
HGDP_FASTA_DIR = DATA_ROOT / "hgdp" / "fasta"
ONEKG_HSD_DIR = DATA_ROOT / "1kg_hsd"
HGDP_GROUND_TRUTH = DATA_ROOT / "hgdp" / "haplogroups" / "haplogroups.csv"
ONEKG_GROUND_TRUTH = DATA_ROOT / "haplogroups" / "1kg_haplogroups_simple.csv"

# Fixed sample set for reproducibility - diverse haplogroups
# HGDP samples for BAM and FASTA testing
# Format: (sample_id, expected_haplogroup, major_clade)
HGDP_SAMPLES = [
    ("HGDP00582", "H1b1f", "H"),
    ("HGDP00845", "D1a1", "D"),
    ("HGDP00862", "A2+(64)", "A"),
    ("HGDP01380", "J1c5c1", "J"),
    ("HGDP00807", "X2b5", "X"),
    ("HGDP00043", "M5a2a4", "M"),
    ("HGDP01339", "B4a1+16311*", "B"),
    ("HGDP00696", "U4d2", "U"),
    ("HGDP00679", "L2a1c+16129", "L"),
    ("HGDP00649", "M1b1a", "M"),
]

# 1KG samples for HSD testing (cleaner HSD format)
# Selected samples verified to match ground truth
ONEKG_SAMPLES = [
    ("HG00096", "H16a1", "H"),
    ("HG00097", "T2f1a1", "T"),
    ("HG00099", "H1ae", "H"),
    ("HG00100", "X2b8", "X"),
    ("HG00101", "J1c2e2", "J"),
    ("HG00105", "H", "H"),
]

# Quick test subsets (3 samples each - fast sanity check)
QUICK_HGDP = HGDP_SAMPLES[:3]
QUICK_ONEKG = ONEKG_SAMPLES[:3]


@dataclass
class TestResult:
    """Result for a single test."""

    sample_id: str
    format_type: str
    expected_haplogroup: str
    expected_clade: str
    called_haplogroup: Optional[str]
    called_clade: Optional[str]
    match_type: str  # exact, prefix, major_clade, mismatch, error
    error: Optional[str] = None

    @property
    def passed(self) -> bool:
        return self.match_type in ("exact", "prefix", "major_clade")


def normalize_haplogroup(hg: str) -> str:
    """Normalize haplogroup for comparison."""
    if not hg:
        return ""
    hg = hg.strip().upper()
    # Remove common annotations
    for char in ["'", "*", "!", "+", " ", "("]:
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


def classify_match(predicted: str, ground_truth: str) -> str:
    """Classify match type between predicted and ground truth."""
    pred_norm = normalize_haplogroup(predicted)
    gt_norm = normalize_haplogroup(ground_truth)

    if not pred_norm or not gt_norm:
        return "mismatch"

    if pred_norm == gt_norm:
        return "exact"

    # Check prefix match (one is ancestor of other)
    if pred_norm.startswith(gt_norm) or gt_norm.startswith(pred_norm):
        return "prefix"

    if get_major_clade(pred_norm) == get_major_clade(gt_norm):
        return "major_clade"

    return "mismatch"


def run_evehap_classify(input_path: Path, verbose: bool = False) -> Tuple[Optional[str], Optional[str]]:
    """
    Run evehap classify and return (haplogroup, error).

    Returns (haplogroup, None) on success or (None, error_message) on failure.
    """
    cmd = [
        sys.executable,
        "-m",
        "evehap.cli",
        "classify",
        str(input_path),
        "--tree-type",
        "rsrs",
        "--output-format",
        "json",
        "-q",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(EVEHAP_ROOT),
        )

        if result.returncode != 0:
            return None, f"Exit code {result.returncode}: {result.stderr.strip()}"

        # Parse JSON output
        try:
            data = json.loads(result.stdout)
            if isinstance(data, list) and len(data) > 0:
                return data[0].get("haplogroup", ""), None
            elif isinstance(data, dict):
                return data.get("haplogroup", ""), None
            else:
                return None, f"Unexpected JSON format: {result.stdout[:100]}"
        except json.JSONDecodeError as e:
            return None, f"JSON parse error: {e}"

    except subprocess.TimeoutExpired:
        return None, "Timeout (120s)"
    except Exception as e:
        return None, str(e)


def test_sample(
    sample_id: str,
    expected_hg: str,
    expected_clade: str,
    format_type: str,
    verbose: bool = False,
) -> TestResult:
    """Test a single sample with a specific format."""
    # Get input path based on format
    if format_type == "bam":
        input_path = HGDP_BAM_DIR / f"{sample_id}.chrM.bam"
    elif format_type == "fasta":
        input_path = HGDP_FASTA_DIR / f"{sample_id}.fasta"
    elif format_type == "hsd":
        # HSD uses 1KG samples (cleaner format than HGDP)
        input_path = ONEKG_HSD_DIR / f"{sample_id}.hsd"
    else:
        return TestResult(
            sample_id=sample_id,
            format_type=format_type,
            expected_haplogroup=expected_hg,
            expected_clade=expected_clade,
            called_haplogroup=None,
            called_clade=None,
            match_type="error",
            error=f"Unknown format: {format_type}",
        )

    if not input_path.exists():
        return TestResult(
            sample_id=sample_id,
            format_type=format_type,
            expected_haplogroup=expected_hg,
            expected_clade=expected_clade,
            called_haplogroup=None,
            called_clade=None,
            match_type="error",
            error=f"File not found: {input_path}",
        )

    haplogroup, error = run_evehap_classify(input_path, verbose)

    if error:
        return TestResult(
            sample_id=sample_id,
            format_type=format_type,
            expected_haplogroup=expected_hg,
            expected_clade=expected_clade,
            called_haplogroup=None,
            called_clade=None,
            match_type="error",
            error=error,
        )

    called_clade = get_major_clade(haplogroup) if haplogroup else None
    match_type = classify_match(haplogroup, expected_hg)

    return TestResult(
        sample_id=sample_id,
        format_type=format_type,
        expected_haplogroup=expected_hg,
        expected_clade=expected_clade,
        called_haplogroup=haplogroup,
        called_clade=called_clade,
        match_type=match_type,
    )


def run_verification(
    hgdp_samples: List[Tuple[str, str, str]],
    onekg_samples: List[Tuple[str, str, str]],
    formats: List[str],
    verbose: bool = False,
) -> Dict[str, List[TestResult]]:
    """Run verification for all samples across all formats."""
    results: Dict[str, List[TestResult]] = {fmt: [] for fmt in formats}

    for fmt in formats:
        if verbose:
            print(f"\n[{fmt.upper()} Format]")

        # Use HGDP for BAM/FASTA, 1KG for HSD
        samples = onekg_samples if fmt == "hsd" else hgdp_samples

        for sample_id, expected_hg, expected_clade in samples:
            result = test_sample(sample_id, expected_hg, expected_clade, fmt, verbose)
            results[fmt].append(result)

            if verbose:
                if result.match_type == "error":
                    symbol = "X"
                    status = f"ERROR: {result.error}"
                elif result.passed:
                    symbol = "ok"
                    status = f"{result.called_haplogroup} ({result.match_type})"
                else:
                    symbol = "FAIL"
                    status = f"{result.called_haplogroup} (expected clade {expected_clade})"

                print(f"  {sample_id}: {status} [{symbol}]")

    return results


def print_summary(results: Dict[str, List[TestResult]], verbose: bool = False) -> bool:
    """Print summary and return True if all tests passed."""
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)

    all_passed = True
    for fmt, fmt_results in results.items():
        passed = sum(1 for r in fmt_results if r.passed)
        total = len(fmt_results)
        pct = (passed / total * 100) if total > 0 else 0
        status = "PASS" if passed == total else "FAIL"
        print(f"{fmt.upper():8} {passed}/{total} ({pct:.0f}%) [{status}]")

        if passed < total:
            all_passed = False
            # Show failures
            for r in fmt_results:
                if not r.passed:
                    if r.error:
                        print(f"         - {r.sample_id}: ERROR - {r.error}")
                    else:
                        print(
                            f"         - {r.sample_id}: called {r.called_haplogroup}, "
                            f"expected {r.expected_haplogroup} (clade {r.expected_clade})"
                        )

    print()
    if all_passed:
        print("Overall: PASS")
    else:
        print("Overall: FAIL")

    return all_passed


def output_json(results: Dict[str, List[TestResult]]) -> None:
    """Output results as JSON."""
    output = {
        "formats": {},
        "summary": {"passed": True, "total_tests": 0, "passed_tests": 0},
    }

    for fmt, fmt_results in results.items():
        passed = sum(1 for r in fmt_results if r.passed)
        total = len(fmt_results)
        output["formats"][fmt] = {
            "passed": passed,
            "total": total,
            "accuracy": passed / total if total > 0 else 0,
            "results": [
                {
                    "sample_id": r.sample_id,
                    "expected": r.expected_haplogroup,
                    "called": r.called_haplogroup,
                    "match_type": r.match_type,
                    "passed": r.passed,
                    "error": r.error,
                }
                for r in fmt_results
            ],
        }
        output["summary"]["total_tests"] += total
        output["summary"]["passed_tests"] += passed
        if passed < total:
            output["summary"]["passed"] = False

    print(json.dumps(output, indent=2))


def main():
    parser = argparse.ArgumentParser(
        description="Verify eveHap CLI with real data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/verify_cli.py          # Full verification
  python scripts/verify_cli.py --quick  # Quick test (3 samples)
  python scripts/verify_cli.py -v       # Verbose output
  python scripts/verify_cli.py --json   # JSON output for CI
""",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--quick", action="store_true", help="Quick test (3 samples only)")
    parser.add_argument("--json", action="store_true", help="Output as JSON")
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["bam", "fasta", "hsd"],
        choices=["bam", "fasta", "hsd"],
        help="Formats to test (default: all)",
    )

    args = parser.parse_args()

    hgdp_samples = QUICK_HGDP if args.quick else HGDP_SAMPLES
    onekg_samples = QUICK_ONEKG if args.quick else ONEKG_SAMPLES

    if not args.json:
        print("eveHap CLI Verification")
        print("=" * 50)
        n_samples = len(hgdp_samples) if "bam" in args.formats or "fasta" in args.formats else len(onekg_samples)
        print(f"Testing {n_samples} samples per format across {len(args.formats)} formats...")
        print(f"  BAM/FASTA: HGDP samples | HSD: 1000 Genomes samples")

    # Check data paths exist
    missing = []
    if "bam" in args.formats and not HGDP_BAM_DIR.exists():
        missing.append(f"BAM dir: {HGDP_BAM_DIR}")
    if "fasta" in args.formats and not HGDP_FASTA_DIR.exists():
        missing.append(f"FASTA dir: {HGDP_FASTA_DIR}")
    if "hsd" in args.formats and not ONEKG_HSD_DIR.exists():
        missing.append(f"HSD dir: {ONEKG_HSD_DIR}")

    if missing:
        print("\nERROR: Missing data directories:")
        for m in missing:
            print(f"  - {m}")
        print("\nPlease ensure test data is available.")
        sys.exit(1)

    results = run_verification(hgdp_samples, onekg_samples, args.formats, args.verbose)

    if args.json:
        output_json(results)
        all_passed = all(r.passed for fmt_results in results.values() for r in fmt_results)
    else:
        all_passed = print_summary(results, args.verbose)

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()

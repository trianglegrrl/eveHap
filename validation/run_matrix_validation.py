#!/usr/bin/env python3
"""
Run validation matrix: Input Format × Tree Type × Haplogrep3 Comparison

Tests all combinations of:
  - Input formats: VCF, BAM, FASTA, 23andMe (microarray)
  - Tree types: rcrs, rsrs

Compares results with Haplogrep3 and produces a detailed JSON report.

Usage:
    python -m validation.run_matrix_validation
    python -m validation.run_matrix_validation --max-samples 50 --workers 8
"""

import argparse
import json
import random
import subprocess
import sys
import tempfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

# Suppress pysam warnings about duplicate header entries
warnings.filterwarnings("ignore", category=UserWarning)
pysam.set_verbosity(0)

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import csv

from evehap.adapters.bam import BAMAdapter
from evehap.adapters.fasta import FASTAAdapter
from evehap.adapters.hsd import HSDAdapter
from evehap.adapters.microarray import MicroarrayAdapter
from evehap.adapters.vcf import VCFAdapter
from evehap.cli import TREE_TYPES
from evehap.core.classifier import Classifier
from evehap.core.phylotree import Phylotree

# Constants
VALIDATION_ROOT = Path(__file__).parent
EVEHAP_ROOT = VALIDATION_ROOT.parent

# Haplogrep3 configuration
HAPLOGREP3_CMD = "/home/a/anaconda3/envs/haplogrep3/bin/haplogrep3"
HAPLOGREP3_TREES = {
    "rcrs": "phylotree-fu-rcrs@1.2",
    "rsrs": "phylotree-rsrs@17.0",  # Fixed: use 17.0 (available version)
    "mitoleaf": "phylotree-fu-rcrs@1.2",  # Compare mitoleaf eveHap vs rCRS Haplogrep3
}

# Mitoleaf tree configuration (not in TREE_TYPES, loaded separately)
MITOLEAF_TREE_PATH = EVEHAP_ROOT / "data" / "phylotree" / "mitoleaf" / "tree.json"
MITOLEAF_REFERENCE_PATH = EVEHAP_ROOT / "data" / "reference" / "rCRS.fasta"

# Ground truth file paths
AADR_GROUND_TRUTH = Path("/home/a/mtdna_poc/data/aadr_ground_truth.csv")
KG_GROUND_TRUTH = Path("/home/a/mtdna_poc/data/haplogroups/1kg_haplogroups_simple.csv")
HGDP_GROUND_TRUTH = Path("/home/a/mtdna_poc/data/hgdp/haplogroups/haplogroups.csv")


def load_ground_truth() -> Dict[str, str]:
    """Load ground truth haplogroups from all sources.

    Returns dict mapping sample_id -> haplogroup
    """
    ground_truth = {}

    # Load AADR ground truth (ERR accession -> haplogroup)
    if AADR_GROUND_TRUTH.exists():
        with open(AADR_GROUND_TRUTH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("match_status") == "matched" and row.get("haplogroup"):
                    ground_truth[row["err_accession"]] = row["haplogroup"]

    # Load 1KG ground truth (sample ID -> haplogroup)
    # NOTE: This is from Haplogrep3, so not independent ground truth
    if KG_GROUND_TRUTH.exists():
        with open(KG_GROUND_TRUTH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("Haplogroup"):
                    ground_truth[row["SampleID"]] = row["Haplogroup"]

    # Load HGDP ground truth (sample ID -> haplogroup)
    # NOTE: This is from Haplogrep3, so not independent ground truth
    if HGDP_GROUND_TRUTH.exists():
        with open(HGDP_GROUND_TRUTH) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get("haplogroup"):
                    ground_truth[row["sample_id"]] = row["haplogroup"]

    return ground_truth


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class MatrixCell:
    """Results for one cell in the validation matrix."""
    format: str  # vcf, bam, fasta, 23andme
    tree_type: str  # rcrs, rsrs

    samples_tested: int = 0
    successful: int = 0
    failed: int = 0

    # Accuracy vs Haplogrep3 (all samples)
    evehap_hp3_exact: int = 0
    evehap_hp3_major: int = 0
    evehap_hp3_mismatch: int = 0
    evehap_hp3_exact_pct: float = 0.0
    evehap_hp3_major_pct: float = 0.0

    # Accuracy vs Ground Truth (only for samples with ground truth)
    ground_truth_total: int = 0
    evehap_gt_exact: int = 0
    evehap_gt_major: int = 0
    evehap_gt_mismatch: int = 0
    evehap_gt_exact_pct: float = 0.0
    evehap_gt_major_pct: float = 0.0
    hp3_gt_exact: int = 0
    hp3_gt_major: int = 0
    hp3_gt_mismatch: int = 0
    hp3_gt_exact_pct: float = 0.0
    hp3_gt_major_pct: float = 0.0

    # High-confidence breakdown (confidence > 0.8)
    high_conf_total: int = 0
    high_conf_exact: int = 0
    high_conf_exact_pct: float = 0.0
    high_conf_major: int = 0
    high_conf_major_pct: float = 0.0

    # Low-confidence breakdown (confidence <= 0.8)
    low_conf_total: int = 0
    low_conf_exact: int = 0
    low_conf_exact_pct: float = 0.0

    # Coverage stats (for BAMs)
    mean_coverage_pct: float = 0.0  # Average % of mtDNA covered
    mean_depth: float = 0.0         # Average read depth at covered positions
    min_coverage_pct: float = 0.0
    max_coverage_pct: float = 0.0

    # Timing
    total_time_seconds: float = 0.0
    mean_time_per_sample_ms: float = 0.0

    # Sample details
    samples: List[Dict[str, Any]] = field(default_factory=list)
    mismatches: List[Dict[str, Any]] = field(default_factory=list)
    ground_truth_disagreements: List[Dict[str, Any]] = field(default_factory=list)


@dataclass
class MatrixResult:
    """Complete matrix validation result."""
    timestamp: str
    evehap_version: str

    # Configuration
    max_samples: Optional[int]
    random_seed: Optional[int]
    workers: int

    # Total counts
    total_cells: int
    total_samples: int
    total_time_seconds: float

    # Matrix cells
    cells: List[MatrixCell] = field(default_factory=list)

    # Summary by dimension
    by_format: Dict[str, Dict[str, float]] = field(default_factory=dict)
    by_tree_type: Dict[str, Dict[str, float]] = field(default_factory=dict)


# =============================================================================
# Haplogrep3 Functions
# =============================================================================


def run_haplogrep3(variants: List[str], tree_type: str) -> Optional[str]:
    """Run Haplogrep3 on a list of variants.

    Args:
        variants: List of variant strings like "73G", "263G"
        tree_type: 'rcrs' or 'rsrs'

    Returns:
        Haplogroup called by Haplogrep3, or None on error
    """
    if not Path(HAPLOGREP3_CMD).exists():
        return None

    tree = HAPLOGREP3_TREES.get(tree_type)
    if not tree:
        return None

    with tempfile.NamedTemporaryFile(mode='w', suffix='.hsd', delete=False) as f:
        hsd_line = f"sample\t1-16569\t?\t" + "\t".join(variants)
        f.write(hsd_line + "\n")
        hsd_path = f.name

    result_path = hsd_path.replace('.hsd', '.txt')

    try:
        subprocess.run([
            HAPLOGREP3_CMD,
            "classify",
            "--tree", tree,
            "--input", hsd_path,
            "--output", result_path,
        ], capture_output=True, text=True, timeout=30)

        if Path(result_path).exists():
            with open(result_path) as f:
                lines = f.readlines()
                if len(lines) > 1:
                    parts = lines[1].strip().replace('"', '').split('\t')
                    return parts[1] if len(parts) > 1 else None
    except Exception:
        pass
    finally:
        Path(hsd_path).unlink(missing_ok=True)
        Path(result_path).unlink(missing_ok=True)

    return None


def get_major_clade(hg: str) -> str:
    """Extract major clade (first letter or two) from haplogroup."""
    if not hg:
        return ""
    hg = hg.strip().upper()
    for i, char in enumerate(hg):
        if char.isdigit():
            return hg[:i] if i > 0 else hg[0]
    return hg[0] if hg else ""


def compare_haplogroups(evehap: str, hp3: str) -> str:
    """Compare evehap and haplogrep3 calls.

    Returns: 'exact', 'major', or 'mismatch'
    """
    if not evehap or not hp3:
        return "mismatch"

    e = evehap.strip().upper()
    h = hp3.strip().upper()

    if e == h:
        return "exact"

    if get_major_clade(e) == get_major_clade(h):
        return "major"

    return "mismatch"


# =============================================================================
# Classification Functions
# =============================================================================


def classify_vcf_sample(
    sample_id: str,
    vcf_path: str,
    classifier: Classifier,
    reference_path: str,
    tree_type: str,
) -> Tuple[str, float, List[str]]:
    """Classify VCF sample and return (haplogroup, confidence, variants)."""
    adapter = VCFAdapter(reference_path=reference_path, sample_id=sample_id)
    profile = adapter.extract_profile(vcf_path)
    result = classifier.classify(profile)

    # Extract variants for Haplogrep3
    variants = [
        f"{pos}{obs.major_allele}"
        for pos, obs in sorted(profile.observations.items())
        if obs.is_variant and obs.major_allele != "-"
    ]

    return result.haplogroup, result.confidence, variants


def classify_bam_sample(
    bam_path: str,
    classifier: Classifier,
    reference_path: str,
    tree_type: str,
) -> Tuple[str, float, List[str], Dict[str, float]]:
    """Classify BAM sample and return (haplogroup, confidence, variants, coverage_info)."""
    adapter = BAMAdapter(reference_path=reference_path)
    profile = adapter.extract_profile(bam_path)
    result = classifier.classify(profile)

    # Extract variants for Haplogrep3
    variants = [
        f"{pos}{obs.major_allele}"
        for pos, obs in sorted(profile.observations.items())
        if obs.is_variant and obs.major_allele != "-"
    ]

    # Calculate coverage info
    coverage_pct = 100 * len(profile.observations) / 16569
    depths = [obs.depth for obs in profile.observations.values() if obs.depth is not None]
    mean_depth = sum(depths) / len(depths) if depths else 0.0

    coverage_info = {
        "coverage_pct": round(coverage_pct, 2),
        "mean_depth": round(mean_depth, 2),
        "positions_covered": len(profile.observations),
    }

    return result.haplogroup, result.confidence, variants, coverage_info


def classify_23andme_sample(
    file_path: str,
    classifier: Classifier,
    tree_type: str,
) -> Tuple[str, float, List[str]]:
    """Classify 23andMe sample and return (haplogroup, confidence, variants)."""
    adapter = MicroarrayAdapter()
    profile = adapter.extract_profile(file_path)
    profile.assumes_full_coverage = True
    result = classifier.classify(profile)

    # Extract variants for Haplogrep3
    variants = [
        f"{pos}{obs.major_allele}"
        for pos, obs in sorted(profile.observations.items())
        if obs.major_allele and obs.major_allele not in ("-", "N")
    ]

    return result.haplogroup, result.confidence, variants


def classify_hsd_sample(
    file_path: str,
    sample_id: str,
    classifier: Classifier,
    tree_type: str,
) -> Tuple[str, float, List[str]]:
    """Classify HSD sample and return (haplogroup, confidence, variants)."""
    adapter = HSDAdapter()
    profile = adapter.extract_profile(file_path, sample_id=sample_id)
    profile.assumes_full_coverage = True
    result = classifier.classify(profile)

    # Extract variants for Haplogrep3
    variants = [
        f"{pos}{obs.major_allele}"
        for pos, obs in sorted(profile.observations.items())
        if obs.major_allele and obs.major_allele not in ("-", "N")
    ]

    return result.haplogroup, result.confidence, variants


def classify_fasta_sample(
    file_path: str,
    classifier: Classifier,
    reference_path: str,
    tree_type: str,
) -> Tuple[str, float, List[str]]:
    """Classify FASTA sample and return (haplogroup, confidence, variants)."""
    adapter = FASTAAdapter(reference_path=reference_path)
    profile = adapter.extract_profile(file_path)
    result = classifier.classify(profile)

    # Extract variants for Haplogrep3
    variants = [
        f"{pos}{obs.major_allele}"
        for pos, obs in sorted(profile.observations.items())
        if obs.is_variant and obs.major_allele != "-"
    ]

    return result.haplogroup, result.confidence, variants


# =============================================================================
# Matrix Validation
# =============================================================================


def run_matrix_cell(
    format_type: str,
    tree_type: str,
    config: Dict[str, Any],
    max_samples: Optional[int],
    random_seed: Optional[int],
    workers: int = 4,
    ground_truth: Optional[Dict[str, str]] = None,
) -> MatrixCell:
    """Run validation for one cell in the matrix with parallel processing."""
    cell = MatrixCell(format=format_type, tree_type=tree_type)
    ground_truth = ground_truth or {}

    # Get tree config (mitoleaf is handled separately)
    if tree_type == "mitoleaf":
        tree_path = str(MITOLEAF_TREE_PATH)
        weights_path = None  # Mitoleaf doesn't use weights
        reference_path = str(MITOLEAF_REFERENCE_PATH)
    else:
        tree_config = TREE_TYPES[tree_type]
        tree_path = str(tree_config["tree"])
        weights_path = str(tree_config["weights"]) if tree_config["weights"].exists() else None
        reference_path = str(tree_config["reference"])

    # Load reference sequence
    with pysam.FastaFile(reference_path) as fasta:
        reference = fasta.fetch(fasta.references[0])

    # Load phylotree and classifier
    phylotree = Phylotree.load(tree_path, None)
    classifier = Classifier(
        phylotree,
        reference=reference,
        weights_path=weights_path,
    )

    # Get samples based on format and tree type
    samples = get_samples_for_format(format_type, tree_type, config, max_samples, random_seed)

    if not samples:
        print(f"    No samples found for {format_type}", flush=True)
        return cell

    # HGDP BAMs are ultra-high depth (~11,000x) and I/O-bound
    # Sequential processing avoids disk contention
    effective_workers = 1 if format_type == "hgdp" else workers

    if format_type == "hgdp":
        print(f"    Testing {len(samples)} samples sequentially (I/O-bound)...", flush=True)
    else:
        print(f"    Testing {len(samples)} samples with {effective_workers} workers...", flush=True)
    start_time = time.perf_counter()

    def process_sample(sample_info: Dict[str, Any]) -> Dict[str, Any]:
        """Process a single sample (for parallel execution)."""
        coverage_info = None  # Only set for BAM/HGDP samples
        try:
            if format_type in ("vcf", "hgdp_vcf", "aadr_vcf"):
                vcf_path = sample_info["vcf_path"]
                # For single-sample VCFs (aadr_vcf), use None to take first sample
                # For multi-sample VCFs (vcf, hgdp_vcf), use the sample_id
                if format_type == "aadr_vcf":
                    sample_id = None  # Single-sample VCF, use first sample
                else:
                    sample_id = sample_info["sample_id"]
                hg, conf, variants = classify_vcf_sample(
                    sample_id, vcf_path, classifier, reference_path, tree_type
                )
                # For aadr_vcf, use filename as sample_id for reporting
                if format_type == "aadr_vcf":
                    sample_id = sample_info["sample_id"]
            elif format_type in ("bam", "hgdp"):
                bam_path = sample_info["bam_path"]
                hg, conf, variants, coverage_info = classify_bam_sample(
                    bam_path, classifier, reference_path, tree_type
                )
                sample_id = Path(bam_path).stem.split(".")[0]
            elif format_type == "23andme":
                file_path = sample_info["file_path"]
                hg, conf, variants = classify_23andme_sample(
                    file_path, classifier, tree_type
                )
                sample_id = Path(file_path).stem
            elif format_type in ("hsd", "hgdp_hsd", "1kg_hsd"):
                file_path = sample_info["file_path"]
                sample_id = sample_info.get("sample_id", Path(file_path).stem)
                hg, conf, variants = classify_hsd_sample(
                    file_path, sample_id, classifier, tree_type
                )
            elif format_type in ("fasta", "hgdp_fasta", "1kg_fasta"):
                file_path = sample_info["file_path"]
                hg, conf, variants = classify_fasta_sample(
                    file_path, classifier, reference_path, tree_type
                )
                sample_id = Path(file_path).stem
            else:
                return {"skip": True}

            # Run Haplogrep3
            hp3_hg = run_haplogrep3(variants, tree_type)
            match_type = compare_haplogroups(hg, hp3_hg) if hp3_hg else "unknown"

            result = {
                "sample_id": sample_id,
                "evehap": hg,
                "evehap_conf": round(conf, 4),
                "haplogrep3": hp3_hg or "N/A",
                "match": match_type,
                "match_type": match_type,
                "variants": variants,
            }
            if coverage_info:
                result["coverage"] = coverage_info
            return result
        except Exception as e:
            return {
                "sample_id": str(sample_info),
                "error": str(e),
            }

    # Process samples in parallel (or sequentially for I/O-bound HGDP)
    with ThreadPoolExecutor(max_workers=effective_workers) as executor:
        futures = {executor.submit(process_sample, s): s for s in samples}

        completed = 0
        for future in as_completed(futures):
            completed += 1
            # Print progress every 10 samples or at milestones
            if completed % 10 == 0 or completed == len(samples):
                elapsed = time.perf_counter() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                remaining = (len(samples) - completed) / rate if rate > 0 else 0
                exact_pct = 100 * cell.evehap_hp3_exact / cell.successful if cell.successful > 0 else 0
                print(f"      {completed}/{len(samples)} ({rate:.1f}/s, ~{remaining:.0f}s left) exact:{exact_pct:.0f}%", flush=True)

            result = future.result()

            if result.get("skip"):
                continue
            if "error" in result:
                cell.failed += 1
                cell.samples.append(result)
                continue

            # Ground truth comparison
            sample_id = result["sample_id"]
            gt_hg = ground_truth.get(sample_id, "")
            result["ground_truth"] = gt_hg

            if gt_hg:
                evehap_gt_match = compare_haplogroups(result["evehap"], gt_hg)
                hp3_gt_match = compare_haplogroups(result["haplogrep3"], gt_hg) if result["haplogrep3"] != "N/A" else "mismatch"
                result["evehap_gt_match"] = evehap_gt_match
                result["hp3_gt_match"] = hp3_gt_match

                cell.ground_truth_total += 1

                # Track eveHap vs ground truth
                if evehap_gt_match == "exact":
                    cell.evehap_gt_exact += 1
                elif evehap_gt_match == "major":
                    cell.evehap_gt_major += 1
                else:
                    cell.evehap_gt_mismatch += 1

                # Track Haplogrep3 vs ground truth
                if hp3_gt_match == "exact":
                    cell.hp3_gt_exact += 1
                elif hp3_gt_match == "major":
                    cell.hp3_gt_major += 1
                else:
                    cell.hp3_gt_mismatch += 1

                # Track disagreements with ground truth
                if evehap_gt_match == "mismatch" or hp3_gt_match == "mismatch":
                    cell.ground_truth_disagreements.append(result)

            cell.samples.append(result)
            cell.successful += 1

            match_type = result["match"]
            conf = result.get("evehap_conf", 0)
            is_high_conf = conf > 0.8

            # Track by confidence level
            if is_high_conf:
                cell.high_conf_total += 1
                if match_type == "exact":
                    cell.high_conf_exact += 1
                elif match_type == "major":
                    cell.high_conf_major += 1
            else:
                cell.low_conf_total += 1
                if match_type == "exact":
                    cell.low_conf_exact += 1

            # Overall tracking
            if match_type == "exact":
                cell.evehap_hp3_exact += 1
            elif match_type == "major":
                cell.evehap_hp3_major += 1
            elif match_type == "mismatch":
                cell.evehap_hp3_mismatch += 1
                cell.mismatches.append(result)

    cell.total_time_seconds = time.perf_counter() - start_time
    cell.samples_tested = cell.successful + cell.failed

    # Calculate percentages
    if cell.successful > 0:
        cell.evehap_hp3_exact_pct = round(100 * cell.evehap_hp3_exact / cell.successful, 2)
        cell.evehap_hp3_major_pct = round(100 * (cell.evehap_hp3_exact + cell.evehap_hp3_major) / cell.successful, 2)
        cell.mean_time_per_sample_ms = round(1000 * cell.total_time_seconds / cell.successful, 2)

    # Ground truth percentages
    if cell.ground_truth_total > 0:
        cell.evehap_gt_exact_pct = round(100 * cell.evehap_gt_exact / cell.ground_truth_total, 2)
        cell.evehap_gt_major_pct = round(100 * (cell.evehap_gt_exact + cell.evehap_gt_major) / cell.ground_truth_total, 2)
        cell.hp3_gt_exact_pct = round(100 * cell.hp3_gt_exact / cell.ground_truth_total, 2)
        cell.hp3_gt_major_pct = round(100 * (cell.hp3_gt_exact + cell.hp3_gt_major) / cell.ground_truth_total, 2)

    # High/low confidence percentages
    if cell.high_conf_total > 0:
        cell.high_conf_exact_pct = round(100 * cell.high_conf_exact / cell.high_conf_total, 2)
        cell.high_conf_major_pct = round(100 * (cell.high_conf_exact + cell.high_conf_major) / cell.high_conf_total, 2)
    if cell.low_conf_total > 0:
        cell.low_conf_exact_pct = round(100 * cell.low_conf_exact / cell.low_conf_total, 2)

    # Coverage stats (for BAM/HGDP formats)
    coverage_samples = [s for s in cell.samples if "coverage" in s]
    if coverage_samples:
        coverages = [s["coverage"]["coverage_pct"] for s in coverage_samples]
        depths = [s["coverage"]["mean_depth"] for s in coverage_samples]
        cell.mean_coverage_pct = round(sum(coverages) / len(coverages), 2)
        cell.mean_depth = round(sum(depths) / len(depths), 2)
        cell.min_coverage_pct = round(min(coverages), 2)
        cell.max_coverage_pct = round(max(coverages), 2)

    return cell


def get_samples_for_format(
    format_type: str,
    tree_type: str,
    config: Dict[str, Any],
    max_samples: Optional[int],
    random_seed: Optional[int],
) -> List[Dict[str, Any]]:
    """Get sample list for a format type, respecting tree type for directory selection."""
    samples = []
    effective_max = max_samples

    if format_type == "vcf":
        # For VCF, check if there's a tree-specific directory
        if tree_type == "rsrs":
            vcf_path = config.get("vcf_path_rsrs") or config.get("vcf_path")
        else:
            vcf_path = config.get("vcf_path")

        if vcf_path and Path(vcf_path).exists():
            try:
                if Path(vcf_path).is_file():
                    # Multi-sample VCF file
                    vcf = pysam.VariantFile(vcf_path)
                    sample_ids = list(vcf.header.samples)
                    vcf.close()
                    samples = [{"sample_id": sid, "vcf_path": vcf_path} for sid in sample_ids]
                else:
                    # Directory of single-sample VCFs
                    vcf_files = list(Path(vcf_path).glob("*.vcf"))
                    samples = [{"sample_id": f.stem, "vcf_path": str(f)} for f in vcf_files]
            except Exception:
                pass

    elif format_type == "bam":
        # Use tree-specific BAM directory for RSRS
        if tree_type == "rsrs":
            bam_dir = config.get("bam_dir_rsrs") or config.get("bam_dir")
        else:
            bam_dir = config.get("bam_dir")

        if bam_dir and Path(bam_dir).exists():
            bam_files = list(Path(bam_dir).glob("*.bam"))
            samples = [{"bam_path": str(b)} for b in bam_files]

    elif format_type == "hgdp":
        # HGDP BAMs (slow, ~30-60s each) - use separate sample limit
        # HGDP is only available in rCRS
        if tree_type != "rcrs":
            return []
        hgdp_dir = config.get("hgdp_dir")
        if hgdp_dir and Path(hgdp_dir).exists():
            bam_files = list(Path(hgdp_dir).glob("*.bam"))
            samples = [{"bam_path": str(b)} for b in bam_files]
        effective_max = config.get("hgdp_samples", 5)

    elif format_type == "hgdp_vcf":
        # HGDP merged VCF (fast - multi-sample VCF)
        # HGDP is only available in rCRS
        if tree_type != "rcrs":
            return []
        hgdp_vcf = config.get("hgdp_vcf")
        if hgdp_vcf and Path(hgdp_vcf).exists():
            try:
                vcf = pysam.VariantFile(hgdp_vcf)
                sample_ids = list(vcf.header.samples)
                vcf.close()
                samples = [{"sample_id": sid, "vcf_path": hgdp_vcf} for sid in sample_ids]
            except Exception:
                pass

    elif format_type == "hgdp_hsd":
        # HGDP HSD files (fast)
        # HGDP is only available in rCRS
        if tree_type != "rcrs":
            return []
        hgdp_hsd_dir = config.get("hgdp_hsd_dir")
        if hgdp_hsd_dir and Path(hgdp_hsd_dir).exists():
            hsd_files = [f for f in Path(hgdp_hsd_dir).glob("*.hsd") if f.name != "all_samples.hsd"]
            samples = [{"file_path": str(f), "sample_id": f.stem} for f in hsd_files]

    elif format_type == "hgdp_fasta":
        # HGDP FASTA files (fast)
        # HGDP is only available in rCRS
        if tree_type != "rcrs":
            return []
        hgdp_fasta_dir = config.get("hgdp_fasta_dir")
        if hgdp_fasta_dir and Path(hgdp_fasta_dir).exists():
            fasta_files = [f for f in Path(hgdp_fasta_dir).glob("*.fasta") if f.name != "all_samples.fasta"]
            samples = [{"file_path": str(f)} for f in fasta_files]

    elif format_type == "aadr_vcf":
        # AADR individual VCF files
        # Use tree-specific directory (rCRS, RSRS, or mitoleaf uses rCRS VCFs)
        if tree_type == "rsrs":
            aadr_vcf_dir = config.get("aadr_vcf_dir_rsrs")
        else:
            # Both rcrs and mitoleaf use rCRS-aligned VCFs
            aadr_vcf_dir = config.get("aadr_vcf_dir")

        if aadr_vcf_dir and Path(aadr_vcf_dir).exists():
            # Get compressed VCFs (.vcf.gz) preferentially
            vcf_files = list(Path(aadr_vcf_dir).glob("*.vcf.gz"))
            if not vcf_files:
                vcf_files = list(Path(aadr_vcf_dir).glob("*.vcf"))
            samples = [{"sample_id": f.stem.replace(".vcf", ""), "vcf_path": str(f)} for f in vcf_files]

    elif format_type == "23andme":
        # Microarray is only tested with rCRS
        if tree_type != "rcrs":
            return []
        ma_dir = config.get("microarray_dir")
        if ma_dir and Path(ma_dir).exists():
            ma_files = list(Path(ma_dir).glob("*.txt"))
            samples = [{"file_path": str(f)} for f in ma_files]

    elif format_type == "hsd":
        # HSD format - use tree-specific directory
        if tree_type == "rsrs":
            hsd_dir = config.get("hsd_dir_rsrs") or config.get("hsd_dir")
        else:
            hsd_dir = config.get("hsd_dir")

        if hsd_dir and Path(hsd_dir).exists():
            hsd_files = [f for f in Path(hsd_dir).glob("*.hsd") if f.name != "all_samples.hsd"]
            samples = [{"file_path": str(f), "sample_id": f.stem} for f in hsd_files]

    elif format_type == "fasta":
        # FASTA format - use tree-specific directory
        if tree_type == "rsrs":
            fasta_dir = config.get("fasta_dir_rsrs") or config.get("fasta_dir")
        else:
            fasta_dir = config.get("fasta_dir")

        if fasta_dir and Path(fasta_dir).exists():
            fasta_files = [f for f in Path(fasta_dir).glob("*.fasta") if f.name != "all_samples.fasta"]
            if not fasta_files:
                fasta_files = [f for f in Path(fasta_dir).glob("*.fa") if f.name != "all_samples.fa"]
            samples = [{"file_path": str(f)} for f in fasta_files]

    elif format_type == "1kg_hsd":
        # 1KG HSD files (rCRS only, generated from 1kg VCF)
        if tree_type != "rcrs":
            return []
        kg_hsd_dir = config.get("kg_hsd_dir")
        if kg_hsd_dir and Path(kg_hsd_dir).exists():
            hsd_files = [f for f in Path(kg_hsd_dir).glob("*.hsd") if f.name != "all_samples.hsd"]
            samples = [{"file_path": str(f), "sample_id": f.stem} for f in hsd_files]

    elif format_type == "1kg_fasta":
        # 1KG FASTA files (rCRS only, generated from 1kg VCF)
        if tree_type != "rcrs":
            return []
        kg_fasta_dir = config.get("kg_fasta_dir")
        if kg_fasta_dir and Path(kg_fasta_dir).exists():
            fasta_files = [f for f in Path(kg_fasta_dir).glob("*.fasta") if f.name != "all_samples.fasta"]
            samples = [{"file_path": str(f)} for f in fasta_files]

    # Apply sampling
    if effective_max and len(samples) > effective_max:
        if random_seed is not None:
            random.seed(random_seed)
        samples = random.sample(samples, effective_max)

    return samples


def run_matrix_validation(
    config: Dict[str, Any],
    formats: List[str],
    tree_types: List[str],
    max_samples: Optional[int],
    random_seed: Optional[int],
    workers: int,
) -> MatrixResult:
    """Run the full validation matrix."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    result = MatrixResult(
        timestamp=timestamp,
        evehap_version="0.1.0",
        max_samples=max_samples,
        random_seed=random_seed,
        workers=workers,
        total_cells=len(formats) * len(tree_types),
        total_samples=0,
        total_time_seconds=0.0,
    )

    # Load ground truth
    print("Loading ground truth...", flush=True)
    ground_truth = load_ground_truth()
    print(f"  Loaded {len(ground_truth)} sample -> haplogroup mappings", flush=True)

    print("\n" + "=" * 70, flush=True)
    print("VALIDATION MATRIX", flush=True)
    print("=" * 70, flush=True)
    print(f"Formats: {formats}", flush=True)
    print(f"Tree types: {tree_types}", flush=True)
    if max_samples:
        print(f"Max samples: {max_samples}", flush=True)
    print(flush=True)

    start_time = time.perf_counter()

    for tree_type in tree_types:
        for format_type in formats:
            print(f"\n[{format_type} × {tree_type}]", flush=True)

            cell = run_matrix_cell(
                format_type,
                tree_type,
                config,
                max_samples,
                random_seed,
                workers=workers,
                ground_truth=ground_truth,
            )

            result.cells.append(cell)
            result.total_samples += cell.samples_tested

            if cell.successful > 0:
                print(f"    ✓ vs HP3 exact: {cell.evehap_hp3_exact_pct:.1f}% | Major clade: {cell.evehap_hp3_major_pct:.1f}%", flush=True)
                if cell.ground_truth_total > 0:
                    print(f"    ✓ vs Ground Truth ({cell.ground_truth_total} samples):", flush=True)
                    print(f"        eveHap:     {cell.evehap_gt_exact_pct:.1f}% exact, {cell.evehap_gt_major_pct:.1f}% major", flush=True)
                    print(f"        Haplogrep3: {cell.hp3_gt_exact_pct:.1f}% exact, {cell.hp3_gt_major_pct:.1f}% major", flush=True)
                if cell.high_conf_total > 0:
                    print(f"    ✓ High-confidence (>0.8): {cell.high_conf_exact_pct:.1f}% exact ({cell.high_conf_total} samples)", flush=True)
                if cell.low_conf_total > 0:
                    print(f"    ⚠ Low-confidence (≤0.8): {cell.low_conf_exact_pct:.1f}% exact ({cell.low_conf_total} samples)", flush=True)
                # Coverage stats for BAM formats
                if cell.mean_coverage_pct > 0:
                    print(f"    📊 Coverage: {cell.mean_coverage_pct:.1f}% avg (range: {cell.min_coverage_pct:.1f}-{cell.max_coverage_pct:.1f}%), depth: {cell.mean_depth:.1f}x", flush=True)
                print(f"    ⏱ Time: {cell.total_time_seconds:.1f}s ({cell.successful} samples)", flush=True)

    result.total_time_seconds = time.perf_counter() - start_time

    # Aggregate by dimension
    for fmt in formats:
        cells = [c for c in result.cells if c.format == fmt]
        total = sum(c.successful for c in cells)
        exact = sum(c.evehap_hp3_exact for c in cells)
        major = sum(c.evehap_hp3_exact + c.evehap_hp3_major for c in cells)
        result.by_format[fmt] = {
            "total_samples": total,
            "exact_pct": round(100 * exact / total, 2) if total > 0 else 0,
            "major_pct": round(100 * major / total, 2) if total > 0 else 0,
        }

    for tt in tree_types:
        cells = [c for c in result.cells if c.tree_type == tt]
        total = sum(c.successful for c in cells)
        exact = sum(c.evehap_hp3_exact for c in cells)
        major = sum(c.evehap_hp3_exact + c.evehap_hp3_major for c in cells)
        result.by_tree_type[tt] = {
            "total_samples": total,
            "exact_pct": round(100 * exact / total, 2) if total > 0 else 0,
            "major_pct": round(100 * major / total, 2) if total > 0 else 0,
        }

    return result


def print_matrix_table(result: MatrixResult) -> None:
    """Print matrix results as a table."""
    print("\n" + "=" * 70)
    print("MATRIX RESULTS")
    print("=" * 70)

    # Build table
    formats = sorted(set(c.format for c in result.cells))
    tree_types = sorted(set(c.tree_type for c in result.cells))

    # Header
    header = f"{'Format':<12}"
    for tt in tree_types:
        header += f" | {tt:>12}"
    print(header)
    print("-" * len(header))

    # Rows
    for fmt in formats:
        row = f"{fmt:<12}"
        for tt in tree_types:
            cell = next((c for c in result.cells if c.format == fmt and c.tree_type == tt), None)
            if cell and cell.successful > 0:
                row += f" | {cell.evehap_hp3_exact_pct:>6.1f}% ex"
            else:
                row += " | " + "-".center(12)
        print(row)

    print("-" * len(header))

    # Summary by confidence level
    print("\n" + "-" * 70)
    print("CONFIDENCE BREAKDOWN")
    print("-" * 70)

    for cell in result.cells:
        if cell.successful == 0:
            continue
        print(f"\n[{cell.format} × {cell.tree_type}]")
        print(f"  High confidence (>0.8): {cell.high_conf_total}/{cell.successful} samples")
        if cell.high_conf_total > 0:
            print(f"    Exact match: {cell.high_conf_exact_pct:.1f}%")
            print(f"    Major clade: {cell.high_conf_major_pct:.1f}%")
        print(f"  Low confidence (≤0.8):  {cell.low_conf_total}/{cell.successful} samples")
        if cell.low_conf_total > 0:
            print(f"    Exact match: {cell.low_conf_exact_pct:.1f}%")

    # Summary
    print("\n" + "-" * 70)
    print("By Format:")
    for fmt, stats in result.by_format.items():
        print(f"  {fmt}: {stats['exact_pct']:.1f}% exact, {stats['major_pct']:.1f}% major clade")

    print("\nBy Tree Type:")
    for tt, stats in result.by_tree_type.items():
        print(f"  {tt}: {stats['exact_pct']:.1f}% exact, {stats['major_pct']:.1f}% major clade")

    # Note about methodology
    print("\n" + "=" * 70)
    print("NOTE ON METHODOLOGY")
    print("=" * 70)
    print("""
When comparing eveHap vs Haplogrep3 on damaged ancient DNA samples:

- Haplogrep3 defaults to H2a2a1 (rCRS reference) when no haplogroup-defining
  mutations are found, returning ~0.5 quality score
- eveHap uses Kulczynski scoring to find the best-matching haplogroup even
  with sparse/damaged data, returning low confidence (<0.8) for uncertain calls

For high-confidence samples (>0.8), both tools should agree closely.
For low-confidence samples, disagreement is expected and neither tool's
call should be considered definitive.

Recommended reporting metrics:
- High-confidence exact match rate (primary accuracy metric)
- Overall major clade match rate (secondary metric)
""")


def save_matrix_results(result: MatrixResult, output_dir: Path) -> None:
    """Save matrix results to JSON."""
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    json_path = output_dir / f"matrix_validation_{timestamp}.json"

    data = {
        "timestamp": result.timestamp,
        "evehap_version": result.evehap_version,
        "config": {
            "max_samples": result.max_samples,
            "random_seed": result.random_seed,
            "workers": result.workers,
        },
        "summary": {
            "total_cells": result.total_cells,
            "total_samples": result.total_samples,
            "total_time_seconds": round(result.total_time_seconds, 2),
        },
        "by_format": result.by_format,
        "by_tree_type": result.by_tree_type,
        "cells": [],
    }

    for cell in result.cells:
        cell_data = {
            "format": cell.format,
            "tree_type": cell.tree_type,
            "samples_tested": cell.samples_tested,
            "successful": cell.successful,
            "failed": cell.failed,
            # Overall accuracy
            "evehap_hp3_exact": cell.evehap_hp3_exact,
            "evehap_hp3_exact_pct": cell.evehap_hp3_exact_pct,
            "evehap_hp3_major": cell.evehap_hp3_major,
            "evehap_hp3_major_pct": cell.evehap_hp3_major_pct,
            "evehap_hp3_mismatch": cell.evehap_hp3_mismatch,
            # High-confidence accuracy (primary metric)
            "high_confidence": {
                "total": cell.high_conf_total,
                "exact": cell.high_conf_exact,
                "exact_pct": cell.high_conf_exact_pct,
                "major": cell.high_conf_major,
                "major_pct": cell.high_conf_major_pct,
            },
            # Low-confidence samples (uncertain calls)
            "low_confidence": {
                "total": cell.low_conf_total,
                "exact": cell.low_conf_exact,
                "exact_pct": cell.low_conf_exact_pct,
            },
            "total_time_seconds": round(cell.total_time_seconds, 2),
            "mean_time_per_sample_ms": cell.mean_time_per_sample_ms,
            "mismatches": cell.mismatches[:10],  # Limit to 10
        }
        data["cells"].append(cell_data)

    # Add methodology note
    data["methodology_note"] = {
        "summary": "High-confidence (>0.8) accuracy is the primary metric for comparing eveHap vs Haplogrep3.",
        "hp3_behavior": "Haplogrep3 defaults to H2a2a1 (rCRS reference) with ~0.5 quality when no haplogroup-defining mutations are found in damaged/sparse samples.",
        "evehap_behavior": "eveHap uses Kulczynski scoring to find best-matching haplogroup even with sparse data, reporting low confidence (<0.8) for uncertain calls.",
        "recommendation": "For ancient DNA with heavy damage, low-confidence calls from either tool should not be considered definitive.",
    }

    with open(json_path, "w") as f:
        json.dump(data, f, indent=2)

    print(f"\nResults saved to: {json_path}")


def save_tsv_results(result: MatrixResult, output_dir: Path) -> Path:
    """Save all sample calls to TSV for supplementary data."""
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tsv_path = output_dir / f"all_calls_{timestamp}.tsv"

    with open(tsv_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "format", "tree_type", "evehap_haplogroup", "evehap_confidence",
            "haplogrep3_haplogroup", "ground_truth", "evehap_hp3_match", "evehap_gt_match",
            "hp3_gt_match", "variants"
        ])

        for cell in result.cells:
            for sample in cell.samples:
                sample_id = sample.get("sample_id", "")
                evehap = sample.get("evehap", "")
                conf = sample.get("evehap_conf", 0)
                hp3 = sample.get("haplogrep3", "")
                gt = sample.get("ground_truth", "")
                evehap_hp3_match = sample.get("match_type", "")
                evehap_gt_match = sample.get("evehap_gt_match", "")
                hp3_gt_match = sample.get("hp3_gt_match", "")
                variants = " ".join(sample.get("variants", [])[:20])  # Limit variants

                writer.writerow([
                    sample_id, cell.format, cell.tree_type, evehap, f"{conf:.4f}",
                    hp3, gt, evehap_hp3_match, evehap_gt_match, hp3_gt_match, variants
                ])

    print(f"TSV saved to: {tsv_path}")
    return tsv_path


def save_disagreement_files(result: MatrixResult, output_dir: Path) -> Path:
    """Save HSD files and call summaries for disagreements with ground truth."""
    disagreements_dir = output_dir / f"disagreements_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    disagreements_dir.mkdir(parents=True, exist_ok=True)

    summary_path = disagreements_dir / "summary.tsv"

    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "sample_id", "format", "tree_type", "evehap", "evehap_conf",
            "haplogrep3", "ground_truth", "evehap_correct", "hp3_correct", "hsd_file"
        ])

        for cell in result.cells:
            for disagreement in cell.ground_truth_disagreements:
                sample_id = disagreement.get("sample_id", "")
                evehap = disagreement.get("evehap", "")
                conf = disagreement.get("evehap_conf", 0)
                hp3 = disagreement.get("haplogrep3", "")
                gt = disagreement.get("ground_truth", "")
                variants = disagreement.get("variants", [])
                evehap_correct = disagreement.get("evehap_gt_match") in ["exact", "major"]
                hp3_correct = disagreement.get("hp3_gt_match") in ["exact", "major"]

                # Generate HSD file for this sample
                hsd_filename = f"{sample_id}.hsd"
                hsd_path = disagreements_dir / hsd_filename

                with open(hsd_path, "w") as hsd_f:
                    var_str = "\t".join(variants) if variants else ""
                    hsd_f.write(f"{sample_id}\t1-16569\t?\t{var_str}\n")

                # Write call summary file
                summary_file = disagreements_dir / f"{sample_id}_calls.txt"
                with open(summary_file, "w") as sf:
                    sf.write(f"Sample: {sample_id}\n")
                    sf.write(f"Format: {cell.format}\n")
                    sf.write(f"Tree: {cell.tree_type}\n")
                    sf.write(f"{'='*40}\n")
                    sf.write(f"Ground Truth: {gt}\n")
                    sf.write(f"eveHap Call:  {evehap} (confidence: {conf:.4f})\n")
                    sf.write(f"Haplogrep3:   {hp3}\n")
                    sf.write(f"{'='*40}\n")
                    sf.write(f"eveHap correct: {'Yes' if evehap_correct else 'No'}\n")
                    sf.write(f"Haplogrep3 correct: {'Yes' if hp3_correct else 'No'}\n")
                    sf.write(f"\nVariants ({len(variants)}):\n")
                    sf.write("  " + " ".join(variants[:50]) + "\n")
                    if len(variants) > 50:
                        sf.write(f"  ... and {len(variants) - 50} more\n")

                writer.writerow([
                    sample_id, cell.format, cell.tree_type, evehap, f"{conf:.4f}",
                    hp3, gt, evehap_correct, hp3_correct, hsd_filename
                ])

    count = sum(len(c.ground_truth_disagreements) for c in result.cells)
    print(f"Disagreement files saved to: {disagreements_dir} ({count} samples)")
    return disagreements_dir


# =============================================================================
# Main
# =============================================================================


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run eveHap validation matrix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Quick validation with VCFs only
    python -m validation.run_matrix_validation --formats vcf --max-samples 50

    # Full validation with all formats
    python -m validation.run_matrix_validation --formats vcf bam hsd fasta --workers 14

    # RSRS-only validation with AADR data
    python -m validation.run_matrix_validation --formats bam hsd --tree-types rsrs \\
        --bam-dir-rsrs /home/a/mtdna_poc/data/aadr_bams_rsrs
        """,
    )
    # rCRS data sources
    parser.add_argument(
        "--vcf", "-v",
        type=str,
        default="/home/a/mtdna_poc/data/raw/1kg_chrMT.vcf.gz",
        help="Path to VCF file for rCRS (can be multi-sample VCF or directory)",
    )
    parser.add_argument(
        "--bam-dir", "-b",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams",
        help="Directory containing rCRS-aligned BAM files",
    )
    parser.add_argument(
        "--hsd-dir",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_hsd",
        help="Directory containing rCRS HSD files",
    )
    parser.add_argument(
        "--fasta-dir",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_fasta",
        help="Directory containing rCRS FASTA files",
    )

    # RSRS data sources (separate directories for RSRS-aligned data)
    parser.add_argument(
        "--vcf-rsrs",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_rsrs_vcf",
        help="Path to VCF file/dir for RSRS-aligned data",
    )
    parser.add_argument(
        "--bam-dir-rsrs",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_rsrs",
        help="Directory containing RSRS-aligned BAM files",
    )
    parser.add_argument(
        "--hsd-dir-rsrs",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_rsrs_hsd",
        help="Directory containing RSRS HSD files",
    )
    parser.add_argument(
        "--fasta-dir-rsrs",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_rsrs_fasta",
        help="Directory containing RSRS FASTA files",
    )

    # Other data sources
    parser.add_argument(
        "--hgdp-dir",
        type=str,
        default="/home/a/mtdna_poc/data/hgdp/bams",
        help="Directory containing HGDP BAM files (rCRS only)",
    )
    parser.add_argument(
        "--hgdp-vcf",
        type=str,
        default="/home/a/mtdna_poc/data/hgdp/vcfs/hgdp_chrM.vcf.gz",
        help="Path to HGDP merged VCF file (rCRS only)",
    )
    parser.add_argument(
        "--hgdp-hsd-dir",
        type=str,
        default="/home/a/mtdna_poc/data/hgdp/hsd",
        help="Directory containing HGDP HSD files (rCRS only)",
    )
    parser.add_argument(
        "--hgdp-fasta-dir",
        type=str,
        default="/home/a/mtdna_poc/data/hgdp/fasta",
        help="Directory containing HGDP FASTA files (rCRS only)",
    )
    parser.add_argument(
        "--hgdp-samples",
        type=int,
        default=5,
        help="Max HGDP BAM samples (these are slow, ~30-60s each)",
    )
    parser.add_argument(
        "--aadr-vcf-dir",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_vcf",
        help="Directory containing AADR VCF files (rCRS-aligned)",
    )
    parser.add_argument(
        "--aadr-vcf-dir-rsrs",
        type=str,
        default="/home/a/mtdna_poc/data/aadr_bams_rsrs_vcf",
        help="Directory containing AADR VCF files (RSRS-aligned)",
    )
    parser.add_argument(
        "--microarray-dir", "-m",
        type=str,
        default="/home/a/mtdna_poc/data/23andme",
        help="Directory containing 23andMe files (rCRS only)",
    )

    # Validation options
    parser.add_argument(
        "--formats", "-f",
        type=str,
        nargs="+",
        default=["vcf"],
        choices=["vcf", "bam", "hgdp", "hgdp_vcf", "hgdp_hsd", "hgdp_fasta", "aadr_vcf", "23andme", "hsd", "fasta", "1kg_hsd", "1kg_fasta"],
        help="Formats to test (hgdp=BAM, hgdp_vcf/hgdp_hsd/hgdp_fasta for HGDP, aadr_vcf for AADR VCFs)",
    )
    parser.add_argument(
        "--tree-types", "-t",
        type=str,
        nargs="+",
        default=["rcrs", "rsrs"],
        choices=["rcrs", "rsrs", "mitoleaf"],
        help="Tree types to test (mitoleaf uses 6400+ haplogroup tree)",
    )
    parser.add_argument(
        "--max-samples", "-n",
        type=int,
        default=50,
        help="Maximum samples per cell (overridden by --hgdp-samples for HGDP)",
    )
    parser.add_argument(
        "--random-seed", "-s",
        type=int,
        default=42,
        help="Random seed for sample selection",
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=4,
        help="Number of parallel workers",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=EVEHAP_ROOT / "validation_results",
        help="Output directory",
    )
    parser.add_argument(
        "--output-tsv",
        action="store_true",
        help="Output TSV file with all sample calls",
    )
    parser.add_argument(
        "--output-disagreements",
        action="store_true",
        help="Output HSD files and summaries for ground truth disagreements",
    )
    parser.add_argument(
        "--run-full-analysis",
        action="store_true",
        help="Run complete analysis pipeline: matrix validation + AADR reports + mismatch analysis + individual vs merged + supplementary material",
    )

    # 1KG HSD and FASTA directories
    parser.add_argument(
        "--1kg-hsd-dir",
        type=str,
        default="/home/a/mtdna_poc/data/1kg_hsd",
        dest="kg_hsd_dir",
        help="Directory containing 1KG HSD files",
    )
    parser.add_argument(
        "--1kg-fasta-dir",
        type=str,
        default="/home/a/mtdna_poc/data/1kg_fasta",
        dest="kg_fasta_dir",
        help="Directory containing 1KG FASTA files",
    )

    args = parser.parse_args()

    config = {
        # rCRS data
        "vcf_path": args.vcf,
        "bam_dir": args.bam_dir,
        "hsd_dir": args.hsd_dir,
        "fasta_dir": args.fasta_dir,
        # RSRS data
        "vcf_path_rsrs": args.vcf_rsrs,
        "bam_dir_rsrs": args.bam_dir_rsrs,
        "hsd_dir_rsrs": args.hsd_dir_rsrs,
        "fasta_dir_rsrs": args.fasta_dir_rsrs,
        # HGDP data (rCRS only)
        "hgdp_dir": args.hgdp_dir,
        "hgdp_vcf": args.hgdp_vcf,
        "hgdp_hsd_dir": args.hgdp_hsd_dir,
        "hgdp_fasta_dir": args.hgdp_fasta_dir,
        "hgdp_samples": args.hgdp_samples,
        # AADR VCFs
        "aadr_vcf_dir": args.aadr_vcf_dir,
        "aadr_vcf_dir_rsrs": args.aadr_vcf_dir_rsrs,
        # 1KG HSD and FASTA
        "kg_hsd_dir": args.kg_hsd_dir,
        "kg_fasta_dir": args.kg_fasta_dir,
        # Other
        "microarray_dir": args.microarray_dir,
    }

    result = run_matrix_validation(
        config=config,
        formats=args.formats,
        tree_types=args.tree_types,
        max_samples=args.max_samples,
        random_seed=args.random_seed,
        workers=args.workers,
    )

    print_matrix_table(result)
    save_matrix_results(result, args.output_dir)

    if args.output_tsv:
        save_tsv_results(result, args.output_dir)

    if args.output_disagreements:
        save_disagreement_files(result, args.output_dir)

    # Run full analysis pipeline if requested
    if args.run_full_analysis:
        print("\n" + "=" * 70)
        print("RUNNING FULL ANALYSIS PIPELINE")
        print("=" * 70)

        import subprocess
        from pathlib import Path

        project_root = Path(__file__).parent.parent.parent
        scripts_dir = project_root / "scripts"
        validation_results_dir = project_root / "evehap" / "validation_results"

        # Step 1: Generate comprehensive AADR report (individual)
        print("\n[1/5] Generating comprehensive AADR report (individual BAMs)...")
        subprocess.run([
            "python", str(scripts_dir / "generate_aadr_validation_report.py"),
            "--output", str(validation_results_dir / "aadr_comprehensive_calls.tsv"),
            "--generate-missing",
            "--workers", str(args.workers),
        ], cwd=project_root, check=False)

        # Step 2: Generate comprehensive AADR report (merged)
        print("\n[2/5] Generating comprehensive AADR report (merged BAMs)...")
        subprocess.run([
            "python", str(scripts_dir / "generate_aadr_validation_report.py"),
            "--output", str(validation_results_dir / "aadr_comprehensive_calls_merged.tsv"),
            "--use-merged",
            "--generate-missing",
            "--workers", str(args.workers),
        ], cwd=project_root, check=False)

        # Step 3: Run mismatch analysis
        print("\n[3/5] Running mismatch analysis...")
        subprocess.run([
            "python", str(scripts_dir / "analyze_bam_mismatches.py"),
            "--input", str(validation_results_dir / "aadr_comprehensive_calls.tsv"),
            "--output-dir", str(validation_results_dir / "mismatch_analysis"),
        ], cwd=project_root, check=False)

        # Step 4: Run individual vs merged analysis
        print("\n[4/5] Running individual vs merged analysis...")
        individual_metadata = project_root / "data" / "aadr_individual_metadata.tsv"
        if individual_metadata.exists():
            subprocess.run([
                "python", str(scripts_dir / "analyze_individual_vs_merged.py"),
                "--individual-calls", str(validation_results_dir / "aadr_comprehensive_calls.tsv"),
                "--merged-calls", str(validation_results_dir / "aadr_comprehensive_calls_merged.tsv"),
                "--individual-metadata", str(individual_metadata),
                "--output-dir", str(validation_results_dir / "individual_vs_merged_analysis"),
            ], cwd=project_root, check=False)
        else:
            print(f"  Skipping: {individual_metadata} not found (run merge_bams_by_individual.py first)")

        # Step 5: Generate supplementary material
        print("\n[5/5] Generating supplementary material...")
        docs_dir = project_root / "docs"
        docs_dir.mkdir(exist_ok=True)
        subprocess.run([
            "python", str(scripts_dir / "generate_validation_supplementary.py"),
            "--output", str(docs_dir / "supplementary_validation.md"),
            "--mismatch-dir", str(validation_results_dir / "mismatch_analysis"),
            "--individual-dir", str(validation_results_dir / "individual_vs_merged_analysis"),
        ], cwd=project_root, check=False)

        print("\n" + "=" * 70)
        print("FULL ANALYSIS PIPELINE COMPLETE")
        print("=" * 70)
        print(f"Results:")
        print(f"  - Matrix validation: {args.output_dir}/matrix_validation_*.json")
        print(f"  - AADR comprehensive calls: {validation_results_dir}/aadr_comprehensive_calls*.tsv")
        print(f"  - Mismatch analysis: {validation_results_dir}/mismatch_analysis/")
        print(f"  - Individual vs merged: {validation_results_dir}/individual_vs_merged_analysis/")
        print(f"  - Supplementary material: {docs_dir}/supplementary_validation.md")
        print("=" * 70)

    print("\n" + "=" * 70)
    print(f"Completed in {result.total_time_seconds:.1f}s")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())


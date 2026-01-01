"""Pytest configuration and fixtures for eveHap tests."""

import os
import random
from pathlib import Path
from typing import Generator, List

import pytest

# Project root directories
EVEHAP_ROOT = Path(__file__).parent.parent
MTDNA_POC_ROOT = EVEHAP_ROOT.parent


def pytest_addoption(parser: pytest.Parser) -> None:
    """Add custom command line options for integration tests."""
    parser.addoption(
        "--sample-count",
        action="store",
        default=None,
        type=int,
        help="Number of samples to test (default: 10 for fast, all for benchmark)",
    )
    parser.addoption(
        "--random-seed",
        action="store",
        default=None,
        type=int,
        help="Random seed for sample selection (default: None = sequential)",
    )
    parser.addoption(
        "--all-samples",
        action="store_true",
        default=False,
        help="Run on all available samples (overrides --sample-count)",
    )


@pytest.fixture
def sample_count(request: pytest.FixtureRequest) -> int:
    """Return the number of samples to test."""
    if request.config.getoption("--all-samples"):
        return -1  # Signal to use all samples
    count = request.config.getoption("--sample-count")
    return count if count is not None else 10  # Default to 10


@pytest.fixture
def random_seed(request: pytest.FixtureRequest) -> int:
    """Return the random seed for sample selection."""
    seed = request.config.getoption("--random-seed")
    return seed


def select_samples(items: List, count: int, seed: int = None) -> List:
    """Select a subset of samples, optionally randomized.

    Args:
        items: List of items to select from
        count: Number of items to select (-1 for all)
        seed: Random seed (None for sequential selection)

    Returns:
        Selected subset of items
    """
    if count == -1 or count >= len(items):
        return items

    if seed is not None:
        rng = random.Random(seed)
        return rng.sample(items, min(count, len(items)))
    else:
        return items[:count]


@pytest.fixture
def evehap_data_dir() -> Path:
    """Return the evehap data directory."""
    return EVEHAP_ROOT / "data"


@pytest.fixture
def mtdna_poc_data_dir() -> Path:
    """Return the parent mtdna_poc data directory."""
    return MTDNA_POC_ROOT / "data"


@pytest.fixture
def phylotree_path(evehap_data_dir: Path) -> Path:
    """Return path to phylotree data file."""
    return evehap_data_dir / "phylotree" / "phylotree-fu-rcrs.txt"


@pytest.fixture
def reference_path(evehap_data_dir: Path) -> Path:
    """Return path to rCRS reference FASTA."""
    return evehap_data_dir / "reference" / "rCRS.fasta"


@pytest.fixture
def hgdp_bam_dir(mtdna_poc_data_dir: Path) -> Path:
    """Return path to HGDP BAM directory."""
    return mtdna_poc_data_dir / "hgdp" / "bams"


@pytest.fixture
def hgdp_haplogroups_path(mtdna_poc_data_dir: Path) -> Path:
    """Return path to HGDP haplogroups CSV."""
    return mtdna_poc_data_dir / "hgdp" / "haplogroups" / "haplogroups.csv"


@pytest.fixture
def aadr_bam_dir(mtdna_poc_data_dir: Path) -> Path:
    """Return path to AADR BAM directory."""
    return mtdna_poc_data_dir / "aadr_bams"


@pytest.fixture
def onekg_vcf_path(mtdna_poc_data_dir: Path) -> Path:
    """Return path to 1000 Genomes mtDNA VCF."""
    return mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"


@pytest.fixture
def onekg_haplogroups_path(mtdna_poc_data_dir: Path) -> Path:
    """Return path to 1000 Genomes haplogroups CSV."""
    return mtdna_poc_data_dir / "haplogroups" / "1kg_haplogroups_simple.csv"


@pytest.fixture
def sample_bam_path(mtdna_poc_data_dir: Path) -> Path:
    """Return path to a sample BAM file for testing."""
    bam_dir = mtdna_poc_data_dir / "hgdp" / "bams"
    if bam_dir.exists():
        bams = list(bam_dir.glob("*.bam"))
        if bams:
            return bams[0]
    # Fallback to test fixture
    return EVEHAP_ROOT / "data" / "test" / "sample.bam"


@pytest.fixture
def temp_output_dir(tmp_path: Path) -> Path:
    """Return a temporary directory for test outputs."""
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)
    return output_dir


# eveHap - Technical Specification

## Version 0.1.0

**Author**: AI-assisted development
**Date**: December 2025
**Status**: Initial Specification

---

## 1. Overview

### 1.1 Purpose

eveHap is a Python-based mitochondrial DNA (mtDNA) haplogroup classification tool designed to:

1. Accept multiple input formats (BAM, VCF, FASTA, FASTQ, 23andMe, HSD)
2. Classify samples against the Phylotree reference phylogeny
3. Handle ancient DNA with damage-aware filtering
4. Provide uncertainty quantification for low-coverage samples
5. Output results in multiple formats (JSON, TSV, human-readable)

### 1.2 Design Philosophy

- **Unified Interface**: All input formats converge to a common `AlleleProfile` representation
- **Phylotree-First**: Classification uses the authoritative Phylotree reference
- **TDD Approach**: All components are test-driven with >90% coverage
- **Ancient DNA Ready**: First-class support for degraded/damaged samples
- **Python Native**: Pure Python with minimal compiled dependencies

### 1.3 Target Users

- Population geneticists working with modern and ancient DNA
- Forensic scientists analyzing mtDNA profiles
- Consumer genetics enthusiasts with 23andMe/Ancestry data
- Ancient DNA researchers needing damage-aware classification

---

## 2. Architecture

### 2.1 High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                           CLI / API                                  │
└─────────────────────────────────┬───────────────────────────────────┘
                                  │
┌─────────────────────────────────▼───────────────────────────────────┐
│                         Input Adapters                               │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐       │
│  │   BAM   │ │   VCF   │ │  FASTA  │ │  FASTQ  │ │ 23andMe │ ...   │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘       │
│       └───────────┴───────────┴───────────┴───────────┘             │
└─────────────────────────────────┬───────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        AlleleProfile                                 │
│   - Position → Allele observations                                   │
│   - Depth, quality, source metadata                                  │
└─────────────────────────────────┬───────────────────────────────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    ▼                           ▼
┌───────────────────────────────┐   ┌───────────────────────────────┐
│        Phylotree              │   │     Damage Filter             │
│   - Tree structure            │   │   - C→T at 5' ends            │
│   - Defining mutations        │   │   - G→A at 3' ends            │
│   - Position weights          │   │   - Configurable thresholds   │
└───────────────┬───────────────┘   └───────────────┬───────────────┘
                └─────────────┬─────────────────────┘
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│                        Classifier                                    │
│   - Kulczynski distance scoring                                      │
│   - Tree traversal algorithm                                         │
│   - Uncertainty quantification                                       │
└─────────────────────────────────┬───────────────────────────────────┘
                                  ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     HaplogroupResult                                 │
│   - Haplogroup assignment                                            │
│   - Confidence score                                                 │
│   - Supporting/conflicting mutations                                 │
│   - QC warnings                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### 2.2 Module Overview

| Module | Purpose | Key Classes |
|--------|---------|-------------|
| `evehap.core.profile` | Data structures for allele observations | `AlleleObservation`, `AlleleProfile` |
| `evehap.core.phylotree` | Phylotree parsing and navigation | `Phylotree`, `PhyloNode`, `Mutation` |
| `evehap.core.classifier` | Classification algorithms | `Classifier`, `KulczynskiScorer`, `TreeTraverser` |
| `evehap.core.damage` | Ancient DNA damage filtering | `DamageFilter`, `DamageStats` |
| `evehap.adapters.*` | Input format adapters | `BAMAdapter`, `VCFAdapter`, etc. |
| `evehap.output.*` | Result formatting | `HaplogroupResult`, `QCReport` |

---

## 3. Core Data Structures

### 3.1 AlleleObservation

Represents allele information at a single mtDNA position.

```python
@dataclass
class AlleleObservation:
    """Allele observation at a single mtDNA position."""

    position: int                    # 1-based position (1-16569)
    ref_allele: str                  # rCRS reference allele
    alleles: Dict[str, float]        # {allele: frequency} e.g., {'G': 0.98, 'A': 0.02}
    depth: Optional[int] = None      # Read depth (None for FASTA/microarray)
    quality: Optional[float] = None  # Mean base quality (None if unavailable)
    source: str = "unknown"          # Source format identifier

    @property
    def major_allele(self) -> str:
        """Return the most frequent allele."""

    @property
    def is_variant(self) -> bool:
        """True if differs from reference allele."""

    @property
    def is_heteroplasmic(self) -> bool:
        """True if multiple alleles present above threshold (default 10%)."""

    def allele_frequency(self, allele: str) -> float:
        """Get frequency of specific allele, 0.0 if absent."""
```

### 3.2 AlleleProfile

Complete mtDNA profile from any input source.

```python
class AlleleProfile:
    """Complete mtDNA allele profile for a sample."""

    def __init__(self, sample_id: str, source_format: str):
        self.sample_id: str
        self.source_format: str
        self.observations: Dict[int, AlleleObservation]  # position -> observation
        self.covered_positions: Set[int]
        self.metadata: Dict[str, Any]

    def add_observation(self, obs: AlleleObservation) -> None:
        """Add an allele observation at a position."""

    def get_allele(self, position: int) -> Optional[str]:
        """Get major allele at position, None if not covered."""

    def get_observation(self, position: int) -> Optional[AlleleObservation]:
        """Get full observation at position."""

    def is_derived(self, position: int, derived_allele: str) -> Optional[bool]:
        """Check if position has derived allele. None if not covered."""

    def coverage_fraction(self) -> float:
        """Fraction of mtDNA genome covered (0.0-1.0)."""

    def mean_depth(self) -> Optional[float]:
        """Mean read depth across covered positions."""

    def get_variants(self, reference: str) -> List[Tuple[int, str, str]]:
        """Return list of (position, ref, alt) variants vs reference."""
```

### 3.3 Phylotree Structures

```python
@dataclass
class Mutation:
    """A phylogenetically informative mutation."""

    position: int           # 1-based position
    ancestral: str          # Ancestral allele
    derived: str            # Derived allele
    weight: float           # Phylogenetic weight (rarer = higher)
    is_backmutation: bool   # True if reverts to ancestral state

    def matches(self, profile: AlleleProfile) -> Optional[bool]:
        """Check if profile has this mutation. None if position not covered."""


@dataclass
class PhyloNode:
    """A node in the phylotree (represents a haplogroup)."""

    haplogroup: str                      # e.g., "H2a2a1"
    parent: Optional['PhyloNode']        # None for root
    children: List['PhyloNode']          # Child haplogroups
    defining_mutations: List[Mutation]   # Mutations defining this branch

    def get_path_from_root(self) -> List['PhyloNode']:
        """Return path from root to this node."""

    def get_all_expected_mutations(self) -> List[Mutation]:
        """Return all mutations expected for this haplogroup (cumulative from root)."""

    def distance_to(self, other: 'PhyloNode') -> int:
        """Phylogenetic distance (number of branches) to another node."""


class Phylotree:
    """Complete phylotree structure."""

    def __init__(self, tree_file: str, weights_file: Optional[str] = None):
        self.root: PhyloNode
        self.nodes: Dict[str, PhyloNode]  # haplogroup -> node
        self.mutations: Dict[int, List[Mutation]]  # position -> mutations
        self.hotspots: Set[int]  # Hypermutable positions

    @classmethod
    def load(cls, tree_path: str) -> 'Phylotree':
        """Load phylotree from file."""

    def get_haplogroup(self, name: str) -> Optional[PhyloNode]:
        """Get node by haplogroup name."""

    def get_defining_positions(self) -> Set[int]:
        """Get all positions that define haplogroups."""

    def get_mutations_at(self, position: int) -> List[Mutation]:
        """Get all mutations at a position."""

    def is_hotspot(self, position: int) -> bool:
        """Check if position is a known hotspot."""
```

---

## 4. Input Adapters

### 4.1 Adapter Interface

All adapters implement a common interface:

```python
from abc import ABC, abstractmethod

class InputAdapter(ABC):
    """Abstract base class for input format adapters."""

    @abstractmethod
    def can_handle(self, file_path: str) -> bool:
        """Check if this adapter can handle the given file."""

    @abstractmethod
    def extract_profile(self, file_path: str, **kwargs) -> AlleleProfile:
        """Extract AlleleProfile from input file."""

    @property
    @abstractmethod
    def format_name(self) -> str:
        """Return the format name (e.g., 'BAM', 'VCF')."""
```

### 4.2 BAM/CRAM Adapter

```python
class BAMAdapter(InputAdapter):
    """Adapter for BAM/CRAM aligned read files."""

    def __init__(self,
                 reference_path: Optional[str] = None,
                 min_base_quality: int = 20,
                 min_mapping_quality: int = 20):
        """
        Args:
            reference_path: Path to reference FASTA (required for CRAM)
            min_base_quality: Minimum base quality to include
            min_mapping_quality: Minimum mapping quality to include
        """

    def extract_profile(self,
                       file_path: str,
                       region: str = "chrM",
                       damage_filter: Optional['DamageFilter'] = None) -> AlleleProfile:
        """
        Extract profile from BAM/CRAM.

        Uses pysam.count_coverage() for efficient allele counting.
        Optionally applies damage filtering for ancient DNA.
        """
```

### 4.3 VCF Adapter

```python
class VCFAdapter(InputAdapter):
    """Adapter for VCF variant call files."""

    def __init__(self,
                 reference_path: str,
                 sample_id: Optional[str] = None):
        """
        Args:
            reference_path: Path to reference FASTA
            sample_id: Sample to extract (None = first sample)
        """

    def extract_profile(self, file_path: str) -> AlleleProfile:
        """
        Extract profile from VCF.

        Parses GT field for genotypes, AD field for allele depths.
        Positions not in VCF are assumed to match reference.
        """
```

### 4.4 FASTA Adapter

```python
class FASTAAdapter(InputAdapter):
    """Adapter for FASTA consensus sequences."""

    def __init__(self, reference_path: str):
        """
        Args:
            reference_path: Path to rCRS reference FASTA
        """

    def extract_profile(self, file_path: str) -> AlleleProfile:
        """
        Extract profile from FASTA.

        Compares each position to rCRS reference.
        Handles IUPAC ambiguity codes as heteroplasmy.
        """
```

### 4.5 FASTQ Adapter

```python
class FASTQAdapter(InputAdapter):
    """Adapter for raw FASTQ read files."""

    def __init__(self,
                 reference_path: str,
                 aligner: str = "minimap2"):
        """
        Args:
            reference_path: Path to reference FASTA
            aligner: Alignment tool to use
        """

    def extract_profile(self, file_path: str) -> AlleleProfile:
        """
        Extract profile from FASTQ.

        Aligns reads to reference and extracts allele observations.
        """
```

### 4.6 Microarray Adapter (23andMe/AncestryDNA)

```python
class MicroarrayAdapter(InputAdapter):
    """Adapter for consumer genotyping data (23andMe, AncestryDNA)."""

    # Position mapping from rsID to mtDNA position
    RSID_TO_POSITION: Dict[str, int]

    def __init__(self, format_type: str = "auto"):
        """
        Args:
            format_type: '23andme', 'ancestry', or 'auto' for detection
        """

    def extract_profile(self, file_path: str) -> AlleleProfile:
        """
        Extract profile from raw genotype file.

        Parses format-specific columns and maps rsIDs to positions.
        Only mtDNA positions are extracted.
        """
```

### 4.7 HSD Adapter

```python
class HSDAdapter(InputAdapter):
    """Adapter for Haplogrep HSD format."""

    def extract_profile(self, file_path: str, sample_id: Optional[str] = None) -> AlleleProfile:
        """
        Extract profile from HSD file.

        Parses polymorphism notation (e.g., 73G, 315.1C, 522-523d).
        """
```

---

## 5. Classification Algorithms

### 5.1 Kulczynski Scoring (from Haplogrep3)

The Kulczynski measure balances two perspectives:

1. **Haplogroup Weight**: How much of the expected haplogroup signature is present?
2. **Sample Weight**: How much of the sample's variation is explained by the haplogroup?

```
KulczynskiScore = (HaplogroupWeight + SampleWeight) / 2

Where:
  HaplogroupWeight = sum(FoundMutationWeights) / sum(ExpectedMutationWeights)
  SampleWeight = sum(FoundMutationWeights) / sum(SampleMutationWeights)
```

```python
class KulczynskiScorer:
    """Calculate Kulczynski similarity scores."""

    def __init__(self, phylotree: Phylotree):
        self.phylotree = phylotree

    def score(self,
              profile: AlleleProfile,
              haplogroup: PhyloNode) -> float:
        """
        Calculate Kulczynski score for a haplogroup.

        Returns score between 0.0 (no match) and 1.0 (perfect match).
        """

    def score_all(self,
                  profile: AlleleProfile,
                  top_n: int = 10) -> List[Tuple[str, float]]:
        """
        Score all haplogroups and return top N.

        Returns list of (haplogroup, score) tuples sorted by score descending.
        """
```

### 5.2 Tree Traversal Algorithm (from pathPhynder)

For low-coverage ancient DNA, traverse the tree evaluating support at each branch:

```python
class TreeTraverser:
    """Tree traversal classification for sparse data."""

    def __init__(self,
                 phylotree: Phylotree,
                 max_conflicts: int = 3,
                 min_support_ratio: float = 0.7):
        """
        Args:
            phylotree: Phylotree to traverse
            max_conflicts: Max ancestral alleles before stopping
            min_support_ratio: Min derived/(derived+ancestral) to continue
        """

    def traverse(self, profile: AlleleProfile) -> TraversalResult:
        """
        Traverse tree to find best haplogroup placement.

        Starting from root, evaluate each branch's support.
        Stop when conflicts exceed threshold or no children have support.

        Returns:
            TraversalResult with path, final haplogroup, and confidence
        """

@dataclass
class TraversalResult:
    """Result of tree traversal classification."""

    path: List[PhyloNode]           # Path from root to final node
    haplogroup: str                 # Final haplogroup assignment
    confidence: float               # Confidence score (0.0-1.0)
    derived_count: int              # Total derived alleles found
    ancestral_count: int            # Total unexpected ancestral alleles
    uncovered_count: int            # Positions not covered in profile
    stopped_early: bool             # True if traversal stopped before leaf
```

### 5.3 Unified Classifier

```python
class Classifier:
    """Unified haplogroup classifier combining multiple strategies."""

    def __init__(self,
                 phylotree: Phylotree,
                 method: str = "auto",
                 damage_filter: Optional['DamageFilter'] = None):
        """
        Args:
            phylotree: Phylotree reference
            method: 'kulczynski', 'traversal', or 'auto'
            damage_filter: Optional damage filter for ancient DNA
        """

    def classify(self, profile: AlleleProfile) -> HaplogroupResult:
        """
        Classify sample and return haplogroup result.

        If method='auto', chooses strategy based on coverage:
        - High coverage (>80%): Kulczynski scoring
        - Low coverage (<80%): Tree traversal
        """
```

---

## 6. Ancient DNA Damage Filtering

### 6.1 Damage Patterns

Ancient DNA exhibits characteristic damage patterns:
- **C→T substitutions** at 5' ends of reads (deamination of cytosine)
- **G→A substitutions** at 3' ends of reads (deamination on opposite strand)

### 6.2 Damage Filter Implementation

```python
@dataclass
class DamageStats:
    """Statistics about damage patterns in a sample."""

    ct_5prime_rate: float    # C→T rate at 5' positions
    ga_3prime_rate: float    # G→A rate at 3' positions
    total_reads: int         # Total reads analyzed
    damaged_reads: int       # Reads with damage signature

    @property
    def is_ancient(self) -> bool:
        """True if damage patterns suggest ancient DNA."""


class DamageFilter:
    """Filter for ancient DNA damage patterns."""

    def __init__(self,
                 filter_ct_5prime: bool = True,
                 filter_ga_3prime: bool = True,
                 end_positions: int = 3,
                 min_damage_rate: float = 0.1):
        """
        Args:
            filter_ct_5prime: Filter C→T at 5' ends
            filter_ga_3prime: Filter G→A at 3' ends
            end_positions: Number of positions from read ends to filter
            min_damage_rate: Minimum damage rate to apply filtering
        """

    def should_filter(self,
                      position: int,
                      ref: str,
                      alt: str,
                      read_position: int,
                      read_length: int,
                      is_reverse: bool) -> bool:
        """
        Determine if this observation should be filtered.

        Args:
            position: mtDNA position
            ref: Reference allele
            alt: Observed allele
            read_position: Position within read (0-based)
            read_length: Total read length
            is_reverse: True if read is reverse strand

        Returns:
            True if observation matches damage pattern and should be filtered
        """

    def estimate_damage(self, bam_path: str) -> DamageStats:
        """
        Estimate damage rates from BAM file.

        Samples reads to estimate C→T and G→A rates at read ends.
        """
```

---

## 7. Output Structures

### 7.1 HaplogroupResult

```python
@dataclass
class HaplogroupResult:
    """Complete result of haplogroup classification."""

    sample_id: str
    haplogroup: str                          # Best haplogroup assignment
    confidence: float                        # Confidence score (0.0-1.0)
    quality: str                             # 'high', 'medium', 'low'

    # Mutation details
    expected_mutations: List[Mutation]       # Mutations expected for haplogroup
    found_mutations: List[Mutation]          # Mutations found in sample
    missing_mutations: List[Mutation]        # Expected but not found
    extra_mutations: List[Mutation]          # Found but not expected (private)

    # Classification details
    method: str                              # 'kulczynski' or 'traversal'
    alternative_haplogroups: List[Tuple[str, float]]  # Other candidates

    # QC information
    coverage_fraction: float                 # Fraction of mtDNA covered
    mean_depth: Optional[float]              # Mean read depth
    warnings: List[str]                      # QC warnings

    # Ancient DNA specific
    damage_filtered: bool                    # True if damage filter applied
    damage_stats: Optional[DamageStats]      # Damage statistics if available

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""

    def to_tsv_row(self) -> str:
        """Convert to TSV row string."""
```

### 7.2 QC Warnings

```python
QC_WARNINGS = {
    'LOW_COVERAGE': 'Coverage below 50% - reduced confidence',
    'LOW_DEPTH': 'Mean depth below 10x - may have errors',
    'HETEROPLASMY': 'Heteroplasmic positions detected',
    'MISSING_KEY_MUTATIONS': 'Key defining mutations not covered',
    'CONFLICTING_MUTATIONS': 'Multiple conflicting mutations found',
    'DAMAGE_DETECTED': 'Ancient DNA damage patterns detected',
    'HOTSPOT_VARIANT': 'Variant at known hotspot position',
}
```

---

## 8. CLI Interface

### 8.1 Command Structure

```bash
evehap classify INPUT [OPTIONS]
evehap validate INPUT --truth GROUND_TRUTH
evehap info INPUT
evehap version
```

### 8.2 Classify Command

```bash
evehap classify sample.bam \
    --format auto \
    --tree phylotree-fu-rcrs \
    --output results.json \
    --output-format json \
    --damage-filter \
    --method auto \
    --top-n 5
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | str | auto | Input format (auto, bam, vcf, fasta, fastq, 23andme, hsd) |
| `--tree` | str | phylotree-fu-rcrs | Phylotree version to use |
| `--output` | path | stdout | Output file path |
| `--output-format` | str | tsv | Output format (json, tsv, text) |
| `--damage-filter` | flag | False | Apply ancient DNA damage filtering |
| `--method` | str | auto | Classification method (auto, kulczynski, traversal) |
| `--top-n` | int | 5 | Number of alternative haplogroups to report |
| `--reference` | path | built-in | Reference FASTA (rCRS) |
| `--sample` | str | None | Sample ID to extract (for multi-sample files) |

### 8.3 Batch Processing

```bash
evehap classify samples/*.bam \
    --output results.tsv \
    --output-format tsv \
    --parallel 8
```

---

## 9. Test Requirements

### 9.1 Unit Tests

Each module must have corresponding unit tests with >90% coverage:

```
tests/
├── test_core/
│   ├── test_profile.py      # AlleleObservation, AlleleProfile
│   ├── test_phylotree.py    # Phylotree parsing and navigation
│   ├── test_classifier.py   # Scoring and traversal algorithms
│   └── test_damage.py       # Damage filtering
├── test_adapters/
│   ├── test_bam.py          # BAM adapter
│   ├── test_vcf.py          # VCF adapter
│   ├── test_fasta.py        # FASTA adapter
│   ├── test_fastq.py        # FASTQ adapter
│   ├── test_microarray.py   # 23andMe adapter
│   └── test_hsd.py          # HSD adapter
└── test_output/
    ├── test_result.py       # Result structures
    └── test_report.py       # Report generation
```

### 9.2 Integration Tests

Integration tests validate against ground truth datasets:

```python
# test_integration/test_hgdp.py
class TestHGDPIntegration:
    """Test against HGDP ground truth (827 samples)."""

    def test_major_clade_accuracy(self):
        """Major clade accuracy should be >= 95%."""

    def test_full_haplogroup_accuracy(self):
        """Full haplogroup accuracy should be >= 85%."""

    def test_bam_vs_fasta_consistency(self):
        """BAM and FASTA should give same results for same sample."""


# test_integration/test_aadr.py
class TestAADRIntegration:
    """Test against AADR ancient DNA (70 samples)."""

    def test_with_damage_filter(self):
        """Damage-filtered results should match ground truth."""

    def test_low_coverage_handling(self):
        """Low coverage samples should have appropriate uncertainty."""


# test_integration/test_1kg.py
class Test1KGIntegration:
    """Test against 1000 Genomes (2504 samples)."""

    def test_vcf_classification(self):
        """VCF classification should match Haplogrep3 results."""
```

### 9.3 Acceptance Criteria

| Metric | Target | Dataset |
|--------|--------|---------|
| Major clade accuracy | >= 95% | HGDP BAM |
| Full haplogroup accuracy | >= 85% | HGDP BAM |
| VCF vs Haplogrep3 agreement | >= 95% | 1000 Genomes |
| Ancient DNA major clade | >= 90% | AADR |
| Code coverage | >= 90% | All modules |
| All input formats working | 100% | Test fixtures |

---

## 10. Dependencies

### 10.1 Runtime Dependencies

```
pysam>=0.21.0          # BAM/CRAM handling
pyfaidx>=0.7.0         # FASTA indexing
numpy>=1.24.0          # Numerical operations
pandas>=2.0.0          # Data handling
click>=8.0.0           # CLI framework
```

### 10.2 Development Dependencies

```
pytest>=7.0.0          # Testing framework
pytest-cov>=4.0.0      # Coverage reporting
pytest-xdist>=3.0.0    # Parallel test execution
mypy>=1.0.0            # Type checking
black>=23.0.0          # Code formatting
ruff>=0.1.0            # Linting
```

---

## 11. File Formats

### 11.1 Phylotree Format

The phylotree format used by Haplogrep3:

```
# Haplogroup\tParent\tMutations
L0  root    8701G 9540C 10398G ...
L1  L0      3594T 7055G ...
L2  L1      2416C 8206A ...
```

### 11.2 Output TSV Format

```
sample_id   haplogroup  confidence  quality coverage    depth   method  warnings
HGDP00001   H2a2a1     0.98        high    1.0         15234   kulczynski
HGDP00002   L3e2b      0.85        medium  0.65        45      traversal   LOW_COVERAGE
```

### 11.3 Output JSON Format

```json
{
  "sample_id": "HGDP00001",
  "haplogroup": "H2a2a1",
  "confidence": 0.98,
  "quality": "high",
  "coverage_fraction": 1.0,
  "mean_depth": 15234,
  "method": "kulczynski",
  "expected_mutations": [...],
  "found_mutations": [...],
  "missing_mutations": [],
  "extra_mutations": [...],
  "alternative_haplogroups": [
    ["H2a2a", 0.95],
    ["H2a2", 0.92]
  ],
  "warnings": []
}
```

---

## 12. Error Handling

### 12.1 Exception Hierarchy

```python
class EveHapError(Exception):
    """Base exception for eveHap."""
    pass

class InputFormatError(EveHapError):
    """Error reading/parsing input file."""
    pass

class PhylotreeError(EveHapError):
    """Error with phylotree data."""
    pass

class ClassificationError(EveHapError):
    """Error during classification."""
    pass

class ConfigurationError(EveHapError):
    """Invalid configuration."""
    pass
```

### 12.2 Error Messages

All errors should include:
- Clear description of what went wrong
- File/position where error occurred (if applicable)
- Suggested remediation

---

## 13. Future Extensions

### 13.1 Planned Features

- Web API server mode
- Interactive visualization
- Contamination detection
- Population frequency annotation
- Heteroplasmy quantification

### 13.2 Plugin Architecture

Future versions may support plugins for:
- Custom phylotrees
- Additional input formats
- Custom scoring algorithms
- Specialized reports

---

## Appendix A: Reference Data

### A.1 rCRS Reference

The revised Cambridge Reference Sequence (rCRS, NC_012920.1) is the standard mtDNA reference. Length: 16,569 bp.

### A.2 Phylotree Versions

- `phylotree-fu-rcrs@1.2` - Current default, rCRS-based
- `phylotree-17` - Phylotree Build 17

### A.3 Test Data Locations

```
data/hgdp/bams/                 # 827 HGDP BAM files
data/hgdp/haplogroups/          # Ground truth haplogroups
data/aadr_bams/                 # 70 AADR ancient DNA BAMs
data/raw/1kg_chrMT.vcf.gz       # 1000 Genomes mtDNA VCF
data/haplogroups/               # Ground truth files
data/reference/rCRS.fasta       # rCRS reference
```


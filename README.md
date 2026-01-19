# eveHap

**Universal mtDNA Haplogroup Classifier**

eveHap is a Python tool for classifying mitochondrial DNA (mtDNA) haplogroups from various input formats. It supports modern high-coverage sequencing data as well as ancient DNA with damage filtering.

## Features

- **Pipeline-Friendly**: Explicit file paths required - no hidden auto-downloads
- **Offline-Ready**: Bundled `rsrs` and `rcrs` resources work without internet
- **Multiple Input Formats**: BAM, CRAM, VCF, FASTA, HSD, and consumer genotyping files (23andMe, AncestryDNA)
- **Dual Classification Strategy**:
  - **Kulczynski scoring** for high-coverage modern DNA
  - **Tree traversal** for low-coverage/ancient DNA
- **mitoLeaf Phylotree**: Optionally download the complete mitoLeaf phylotree (6400+ haplogroups)
- **Ancient DNA Support**: Damage pattern detection and filtering for C→T/G→A substitutions
- **Quality Metrics**: Confidence scores, coverage statistics, and QC warnings
- **Flexible Output**: TSV, JSON, or human-readable text

## Installation

```bash
cd evehap
pip install -e .
```

## Quick Start

```bash
# Classify using bundled resources (offline-ready, no download required)
evehap classify sample.bam --tree rsrs --reference rcrs

# Classify a VCF file with specific sample
evehap classify samples.vcf.gz --tree rsrs --reference rcrs --sample-id HG00096

# Classify multiple files with output to TSV
evehap classify *.bam --tree rsrs --reference rcrs -o results.tsv

# Classify ancient DNA with damage filtering
evehap classify ancient.bam --tree rsrs --reference rcrs --damage-filter --method traversal

# Download mitoLeaf tree (more haplogroups) and use it
evehap download --outdir ./evehap_data/
evehap classify sample.bam --tree ./evehap_data/mitoleaf_tree.json --reference ./evehap_data/rCRS.fasta

# Estimate ancient DNA damage rates
evehap damage ancient.bam

# Show version and bundled resources
evehap version
```

## Usage

### classify

Classify mtDNA haplogroup from input files.

```bash
evehap classify [OPTIONS] INPUT_FILES...

Options:
  --tree TEXT                     REQUIRED: 'rsrs', 'rcrs', or path to tree file
  --reference TEXT                REQUIRED: 'rcrs' or path to reference FASTA
  --format TEXT                   Input format (auto, bam, vcf, fasta, hsd, microarray)
  -o, --output PATH               Output file (default: stdout)
  --output-format [json|tsv|text] Output format
  --sample-id TEXT                Sample ID for multi-sample VCF/HSD files
  --damage-filter                 Apply ancient DNA damage filtering
  --method [auto|kulczynski|traversal]
                                  Classification method
  --top-n INTEGER                 Number of alternatives to report
  -q, --quiet                     Suppress progress output
```

**Tree Options:**

| Name | Description | Haplogroups | Requires Download |
|------|-------------|-------------|-------------------|
| `rsrs` | Bundled RSRS-based XML tree | ~5400 | No (bundled) |
| `rcrs` | Bundled rCRS-based XML tree | ~2400 | No (bundled) |
| `mitoleaf_tree.json` | mitoLeaf JSON tree | ~6400 | Yes |

```bash
# Use bundled RSRS tree (offline-ready)
evehap classify sample.bam --tree rsrs --reference rcrs

# Use bundled rCRS tree
evehap classify sample.bam --tree rcrs --reference rcrs

# Download and use mitoLeaf tree (most comprehensive)
evehap download --outdir ./evehap_data/
evehap classify sample.bam --tree ./evehap_data/mitoleaf_tree.json --reference ./evehap_data/rCRS.fasta
```

### download

Download phylotree and reference resources.

```bash
evehap download --outdir ./evehap_data/

Options:
  -o, --outdir PATH               REQUIRED: Directory to save downloaded files
  --category [all|reference|phylotree]
                                  Category to download (default: all)
  --resource TEXT                 Specific resource(s) to download
  --force                         Overwrite existing files
  --check-updates                 Check for updates without downloading
  --list                          List available resources
```

```bash
# Download all resources
evehap download --outdir ./evehap_data/

# List available resources
evehap download --outdir ./evehap_data/ --list

# Download only the mitoLeaf tree
evehap download --outdir ./evehap_data/ --resource mitoleaf-tree

# Check for updates
evehap download --outdir ./evehap_data/ --check-updates
```

### info

Show information about an input file.

```bash
evehap info sample.bam
```

### damage

Estimate ancient DNA damage rates from BAM file.

```bash
evehap damage ancient.bam
```

### download

Download reference sequences and phylotree resources.

```bash
evehap download [OPTIONS]

Options:
  --category [all|reference|phylotree]  Category of resources (default: all)
  --resource TEXT                       Specific resource(s) to download
  --force                               Overwrite existing files
  --check-updates                       Check for updates without downloading
  --list                                List available resources
```

**Available Resources:**

| Resource | Source | Description |
|----------|--------|-------------|
| **Reference Sequences** | | |
| `rsrs` | [phylotree.org](https://www.phylotree.org/resources/RSRS.fasta) | Reconstructed Sapiens Reference Sequence |
| `rcrs` | [phylotree.org](https://www.phylotree.org/resources/rCRS.fasta) | Revised Cambridge Reference Sequence |
| **mitoLeaf Phylotree Data** | | |
| `mitoleaf-tree` | [forensicgenomics.github.io](https://forensicgenomics.github.io/mitoLeaf/data/tree.json) | Complete phylotree in JSON format (~15MB) |
| `mitoleaf-motifs` | [forensicgenomics.github.io](https://forensicgenomics.github.io/mitoLeaf/data/hgmotifs.json) | Haplogroup defining mutations |
| `mitoleaf-representatives` | [forensicgenomics.github.io](https://forensicgenomics.github.io/mitoLeaf/data/mito_representatives.csv) | Representative sequences per haplogroup |

**Examples:**

```bash
# Download all resources
evehap download

# Download only reference sequences
evehap download --category reference

# Download only phylotree data
evehap download --category phylotree

# Download specific resource
evehap download --resource mitoleaf-tree

# Check for updates to installed resources
evehap download --check-updates

# Force re-download (update all)
evehap download --force

# List all available resources
evehap download --list
```

### version

Show version information and installed resources.

```bash
evehap version [OPTIONS]

Options:
  --check-updates    Check for resource updates
```

Output shows:
- Package and data directory locations
- Default tree and reference paths
- Bundled phylotrees with sizes
- Downloaded resources with sizes and dates

**Example:**

```bash
$ evehap version
eveHap v0.1.0
Package directory: /path/to/evehap
Data directory: /path/to/evehap/data
Default tree: tree-rsrs.xml
Default reference: rCRS.fasta

Bundled phylotrees:
  ✓ tree-rsrs.xml (1837.7 KB)
  ✓ tree.xml (2369.2 KB)

Downloaded resources:

  [reference]
    ✓ rsrs (16.7 KB, 2025-01-01)
    ✓ rcrs (16.5 KB, 2025-01-01)

  [phylotree]
    ✓ mitoleaf-tree (15293.6 KB, 2025-01-01)
    ...
```

## Python API

```python
from evehap.adapters.bam import BAMAdapter
from evehap.core.classifier import Classifier
from evehap.core.phylotree import Phylotree

# Load phylotree (mitoLeaf JSON or Haplogrep XML - auto-detected)
phylotree = Phylotree.load("data/phylotree/mitoleaf/tree.json")

# Extract profile from BAM
adapter = BAMAdapter(reference_path="data/reference/rCRS.fasta")
profile = adapter.extract_profile("sample.bam")

# Classify
classifier = Classifier(phylotree)
result = classifier.classify(profile)

print(f"Haplogroup: {result.haplogroup}")
print(f"Confidence: {result.confidence:.1%}")
print(f"Quality: {result.quality}")
```

## Supported Input Formats

| Format | Extensions | Description |
|--------|-----------|-------------|
| BAM/CRAM | `.bam`, `.cram` | Aligned sequencing reads |
| VCF | `.vcf`, `.vcf.gz` | Variant call format |
| FASTA | `.fasta`, `.fa`, `.fna` | Consensus sequences |
| HSD | `.hsd` | Haplogrep polymorphism format |
| Microarray | `.txt`, `.csv` | 23andMe, AncestryDNA raw data |

## Classification Methods

### Kulczynski Scoring (default for high-coverage)

Uses the Kulczynski similarity measure to score each haplogroup based on the proportion of expected mutations present in the sample. Best for samples with >80% mtDNA coverage.

### Tree Traversal (default for low-coverage)

Traverses the phylogenetic tree from root, evaluating support for each branch based on derived/ancestral allele counts. Best for ancient DNA or samples with sparse coverage.

## Reference Sequences

eveHap supports two mtDNA reference sequences:

### RSRS (Reconstructed Sapiens Reference Sequence)

The ancestral human mtDNA sequence reconstructed from phylogenetic analysis. The RSRS-based phylotree has the complete haplogroup structure starting from mtMRCA (mitochondrial Most Recent Common Ancestor).

- **Use for**: Phylotree structure, haplogroup classification
- **Source**: [phylotree.org/resources/RSRS.fasta](https://www.phylotree.org/resources/RSRS.fasta)

### rCRS (Revised Cambridge Reference Sequence)

The standard reference for modern mtDNA sequencing. Most BAM/VCF files are aligned to rCRS. It represents haplogroup H2a2a1 (a European lineage).

- **Use for**: Sequence alignment, variant calling
- **Source**: [phylotree.org/resources/rCRS.fasta](https://www.phylotree.org/resources/rCRS.fasta)

**Note**: eveHap automatically handles the translation between rCRS-aligned input files and the RSRS-based phylotree.

## Ancient DNA Damage Filtering

Ancient DNA exhibits characteristic damage patterns:
- C→T substitutions at 5' read ends
- G→A substitutions at 3' read ends

Use `--damage-filter` to automatically detect and filter potentially damaged bases.

## Output Formats

### TSV (default)

```
sample_id	haplogroup	confidence	quality	coverage	depth	method	warnings
HGDP00582	H2a2a1g	1.0000	high	1.0000	19004.1	kulczynski	HETEROPLASMY
```

### JSON

```json
{
  "sample_id": "HGDP00582",
  "haplogroup": "H2a2a1g",
  "confidence": 1.0,
  "quality": "high",
  ...
}
```

### Text

```
Sample: HGDP00582
  Haplogroup: H2a2a1g
  Confidence: 100.0%
  Quality: high
  Method: kulczynski
  Coverage: 100.0%
```

## Pipeline Integration

eveHap is designed to be pipeline-friendly for batch processing on HPC and cloud environments. All resource paths must be explicitly specified - no hidden auto-downloads.

### Using Bundled Resources (Offline-Ready)

```bash
# Bundled resources work without internet access
evehap classify sample.bam --tree rsrs --reference rcrs
```

### Using Downloaded Resources

```bash
# Download resources to a directory
evehap download --outdir ./evehap_data/

# Use downloaded files explicitly
evehap classify sample.bam \
    --tree ./evehap_data/mitoleaf_tree.json \
    --reference ./evehap_data/rCRS.fasta
```

### Nextflow

```bash
# Using bundled resources
nextflow run pipelines/nextflow/main.nf \
    --input 'samples/*.bam' \
    --tree rsrs \
    --reference rcrs \
    --outdir results/

# Using downloaded resources
nextflow run pipelines/nextflow/main.nf \
    --input 'samples/*.bam' \
    --tree ./evehap_data/mitoleaf_tree.json \
    --reference ./evehap_data/rCRS.fasta \
    --outdir results/

# SLURM cluster
nextflow run pipelines/nextflow/main.nf \
    --input 'samples/*.bam' \
    --tree rsrs \
    --reference rcrs \
    -profile slurm
```

### Snakemake

```bash
cd pipelines/snakemake

# Using bundled resources
snakemake --cores 4 --config input_dir=../../samples tree=rsrs reference=rcrs

# SLURM cluster
snakemake --profile profiles/slurm --config tree=rsrs reference=rcrs
```

See `pipelines/README.md` for complete documentation.

## Testing

```bash
# Run all tests
pytest tests/

# Run fast tests only (exclude slow benchmarks)
pytest tests/ -m "not slow"

# Run with verbose output
pytest tests/ -v
```

## Project Structure

```
evehap/
├── evehap/
│   ├── adapters/       # Input format adapters
│   │   ├── bam.py      # BAM/CRAM adapter
│   │   ├── vcf.py      # VCF adapter
│   │   ├── fasta.py    # FASTA adapter
│   │   ├── hsd.py      # HSD adapter
│   │   └── microarray.py  # 23andMe/Ancestry adapter
│   ├── core/           # Core classification logic
│   │   ├── classifier.py  # Classifier algorithms
│   │   ├── damage.py   # Ancient DNA damage filtering
│   │   ├── phylotree.py   # Phylotree parser
│   │   └── profile.py  # AlleleProfile data structures
│   ├── output/         # Output formatting
│   │   ├── report.py   # Report generation
│   │   └── result.py   # Result data structures
│   └── cli.py          # Command-line interface
├── data/
│   ├── phylotree/      # Phylotree data files
│   │   ├── tree-rsrs.xml  # RSRS-based tree (default, complete phylogeny)
│   │   ├── tree.xml       # rCRS-based tree
│   │   └── mitoleaf/      # mitoLeaf data (downloaded)
│   │       ├── tree.json          # Complete phylotree (JSON)
│   │       ├── hgmotifs.json      # Haplogroup motifs
│   │       └── mito_representatives.csv
│   └── reference/      # Reference sequences (downloaded)
│       ├── RSRS.fasta  # Reconstructed Sapiens Reference Sequence
│       └── rCRS.fasta  # Revised Cambridge Reference Sequence
└── tests/              # Test suite
```

## Requirements

- Python 3.8+
- pysam
- click
- numpy

## License

MIT License

## Citation

If you use eveHap in your research, please cite:

```
eveHap: Universal mtDNA Haplogroup Classifier
https://github.com/trianglegrrl/eveHap
```



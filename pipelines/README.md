# eveHap Pipeline Examples

Production-ready pipeline templates for batch mtDNA haplogroup classification.

## Features

- **Pipeline-friendly**: Explicit file paths required - no hidden auto-downloads
- **Offline-ready**: Use bundled `rsrs`/`rcrs` resources or pre-downloaded files
- **HPC-ready**: Profiles for SLURM, PBS, and other schedulers
- **Container support**: Docker and Singularity configurations

## Quick Start

### Using Bundled Resources (Offline-Ready)

```bash
# Bundled resources require NO internet access
evehap classify sample.bam --tree rsrs --reference rcrs
```

### Using Downloaded Resources

```bash
# 1. Download mitoLeaf (more complete tree, 6400+ haplogroups)
evehap download --outdir ./evehap_data/

# 2. Use downloaded files explicitly
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

# Using custom/downloaded resources
nextflow run pipelines/nextflow/main.nf \
    --input 'samples/*.bam' \
    --tree /path/to/mitoleaf_tree.json \
    --reference /path/to/rCRS.fasta \
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
snakemake --cores 4 --config \
    input_dir=../../samples \
    tree=rsrs \
    reference=rcrs

# Using custom resources
snakemake --cores 4 --config \
    input_dir=../../samples \
    tree=/path/to/mitoleaf_tree.json \
    reference=/path/to/rCRS.fasta

# SLURM cluster
snakemake --profile profiles/slurm --config \
    input_dir=../../samples \
    tree=rsrs \
    reference=rcrs
```

## Download Resources for Offline Use

   ```bash
# Download all resources to a directory
evehap download --outdir ./evehap_data/

# Download specific resources
evehap download --outdir ./evehap_data/ --resource mitoleaf-tree --resource rcrs

# Check for updates
evehap download --outdir ./evehap_data/ --check-updates

# Copy to HPC/offline environment
scp -r ./evehap_data/ cluster:~/evehap_data/
   ```

## Configuration Options

### Required Arguments

| Argument | Description | Examples |
|----------|-------------|----------|
| `--tree` | REQUIRED: Phylotree | `rsrs`, `rcrs`, `/path/to/tree.json` |
| `--reference` | REQUIRED: Reference FASTA | `rcrs`, `/path/to/reference.fasta` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--method` | Classification method | `auto` |
| `--damage-filter` | Enable ancient DNA damage filtering | false |
| `--sample-id` | Sample ID for multi-sample VCF | none |
| `-o, --output` | Output file | stdout |
| `--output-format` | Output format: json, tsv, text | tsv |

### Available Trees

| Name | Haplogroups | Type | Notes |
|------|-------------|------|-------|
| `rsrs` | ~5400 | Bundled | RSRS-rooted, offline-ready |
| `rcrs` | ~2400 | Bundled | rCRS-rooted, offline-ready |
| `mitoleaf_tree.json` | ~6400 | Downloaded | Most complete, requires download |

## Output Files

```
results/
├── {sample}.haplogroup.tsv   # Per-sample TSV results
├── {sample}.haplogroup.json  # Per-sample JSON results
├── all_haplogroups.tsv       # Merged results
├── summary.txt               # Distribution statistics
├── damage/                   # Damage analysis (BAM only)
│   └── {sample}.damage.txt
└── logs/                     # Execution logs
```

## HPC Integration

### SLURM

```bash
# Nextflow
nextflow run pipelines/nextflow/main.nf \
    --input 'samples/*.bam' \
    --tree rsrs \
    --reference rcrs \
    -profile slurm

# Snakemake
snakemake --profile profiles/slurm --config tree=rsrs reference=rcrs
```

### PBS/Torque

```bash
nextflow run pipelines/nextflow/main.nf --tree rsrs --reference rcrs -profile pbs
```

### Docker/Singularity

```bash
# Docker
nextflow run pipelines/nextflow/main.nf --tree rsrs --reference rcrs -profile docker

# Singularity (HPC-friendly)
nextflow run pipelines/nextflow/main.nf --tree rsrs --reference rcrs -profile singularity
```

## Example: Large-Scale Analysis

```bash
# 1. Download resources (on machine with internet)
evehap download --outdir /shared/evehap_data/

# 2. Run on 1000 samples with SLURM
nextflow run pipelines/nextflow/main.nf \
    --input '/data/samples/*.bam' \
    --outdir /scratch/results/ \
    --tree /shared/evehap_data/mitoleaf_tree.json \
    --reference /shared/evehap_data/rCRS.fasta \
    -profile slurm \
    -resume

# 3. Check progress
tail -f /scratch/results/reports/execution_report.html
```

## Troubleshooting

### "Error: --tree is required"

All runs require explicit `--tree` and `--reference` arguments:

```bash
# Use bundled resources
evehap classify sample.bam --tree rsrs --reference rcrs

# Or use downloaded files
evehap classify sample.bam --tree ./data/mitoleaf_tree.json --reference ./data/rCRS.fasta
```

### "Tree file not found"

Make sure the path exists or use bundled resource names (`rsrs`, `rcrs`):

```bash
# Check what's bundled
evehap version

# Download resources
evehap download --outdir ./evehap_data/ --list
evehap download --outdir ./evehap_data/
```

### Memory issues with large VCFs

Multi-sample VCFs are memory-intensive. Use `--sample-id` to extract one sample at a time:

```bash
evehap classify large.vcf.gz --tree rsrs --reference rcrs --sample-id SAMPLE001
```

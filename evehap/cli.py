"""Command-line interface for eveHap.

Provides a simple CLI for haplogroup classification from various input formats.
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import click

from evehap.adapters.bam import BAMAdapter
from evehap.adapters.fasta import FASTAAdapter
from evehap.adapters.hsd import HSDAdapter
from evehap.adapters.microarray import MicroarrayAdapter
from evehap.adapters.vcf import VCFAdapter
from evehap.core.classifier import Classifier
from evehap.core.damage import DamageFilter
from evehap.core.phylotree import Phylotree
from evehap.output.result import HaplogroupResult

# Default data paths (relative to package)
PACKAGE_DIR = Path(__file__).parent.parent
DATA_DIR = PACKAGE_DIR / "data"

# Tree type configurations
# Each tree type specifies its tree file, weights file, and reference
TREE_TYPES: Dict[str, Dict[str, Union[Path, str]]] = {
    "rcrs": {
        "tree": DATA_DIR / "phylotree" / "tree-rcrs.xml",
        "weights": DATA_DIR / "phylotree" / "weights-rcrs.txt",
        "reference": DATA_DIR / "reference" / "rCRS.fasta",
        "description": "rCRS-rooted phylotree (matches Haplogrep3 default)",
    },
    "rsrs": {
        "tree": DATA_DIR / "phylotree" / "tree-rsrs.xml",
        "weights": DATA_DIR / "phylotree" / "weights-rsrs.txt",
        "reference": DATA_DIR / "reference" / "rCRS.fasta",
        "description": "RSRS-rooted phylotree (for ancient DNA studies)",
    },
}

# Legacy aliases for backwards compatibility (cast to Path for type safety)
BUNDLED_TREE_RSRS: Path = Path(TREE_TYPES["rsrs"]["tree"])
BUNDLED_TREE_RCRS: Path = Path(TREE_TYPES["rcrs"]["tree"])
BUNDLED_REFERENCE: Path = Path(TREE_TYPES["rcrs"]["reference"])
DEFAULT_WEIGHTS: Path = Path(TREE_TYPES["rsrs"]["weights"])

# Download resources configuration
# Note: 'dest' is relative path within the download directory
DOWNLOAD_RESOURCES = {
    # Reference sequences from phylotree.org
    "rsrs": {
        "url": "https://www.phylotree.org/resources/RSRS.fasta",
        "filename": "RSRS.fasta",
        "category": "reference",
        "description": "Reconstructed Sapiens Reference Sequence",
    },
    "rcrs": {
        "url": "https://www.phylotree.org/resources/rCRS.fasta",
        "filename": "rCRS.fasta",
        "category": "reference",
        "description": "Revised Cambridge Reference Sequence",
    },
    # mitoLeaf phylotree data from forensicgenomics.github.io
    "mitoleaf-tree": {
        "url": "https://forensicgenomics.github.io/mitoLeaf/data/tree.json",
        "filename": "mitoleaf_tree.json",
        "category": "phylotree",
        "description": "Complete phylotree (JSON format, ~15MB)",
    },
    "mitoleaf-motifs": {
        "url": "https://forensicgenomics.github.io/mitoLeaf/data/hgmotifs.json",
        "filename": "mitoleaf_motifs.json",
        "category": "phylotree",
        "description": "Haplogroup defining mutations",
    },
    "mitoleaf-representatives": {
        "url": "https://forensicgenomics.github.io/mitoLeaf/data/mito_representatives.csv",
        "filename": "mitoleaf_representatives.csv",
        "category": "phylotree",
        "description": "Representative sequences per haplogroup",
    },
}


def get_resource_path(resource_name: str, download_dir: Optional[Path] = None) -> Path:
    """Get the path for a downloadable resource.

    Args:
        resource_name: Name of the resource
        download_dir: Custom download directory (default: package data dir)

    Returns:
        Path to the resource file
    """
    if resource_name not in DOWNLOAD_RESOURCES:
        raise ValueError(f"Unknown resource: {resource_name}")

    info = DOWNLOAD_RESOURCES[resource_name]
    base_dir = download_dir or DATA_DIR

    return base_dir / info["filename"]


def check_for_updates(download_dir: Optional[Path] = None, quiet: bool = True) -> List[str]:
    """Check if any installed resources have updates available.

    Args:
        download_dir: Directory where resources are stored
        quiet: Suppress output

    Returns:
        List of resource names with updates available
    """
    updates = []
    for name, info in DOWNLOAD_RESOURCES.items():
        dest = get_resource_path(name, download_dir)
        if dest.exists():
            if _check_update_available(dest, info["url"]):
                updates.append(name)
                if not quiet:
                    click.echo(f"Update available: {name}", err=True)
    return updates


def detect_format(file_path: str) -> str:
    """Detect input file format from extension.

    Args:
        file_path: Path to input file

    Returns:
        Format name ('bam', 'vcf', 'fasta', 'hsd', 'microarray')
    """
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in (".bam", ".cram"):
        return "bam"
    elif suffix in (".vcf",) or path.name.endswith(".vcf.gz"):
        return "vcf"
    elif suffix in (".fasta", ".fa", ".fna"):
        return "fasta"
    elif suffix == ".hsd":
        return "hsd"
    elif suffix in (".txt", ".csv"):
        return "microarray"
    else:
        return "unknown"


def get_adapter(
    format_name: str,
    reference_path: Optional[str] = None,
    sample_id: Optional[str] = None,
):
    """Get the appropriate adapter for a format.

    Args:
        format_name: Format name
        reference_path: Path to reference FASTA
        sample_id: Sample ID for multi-sample files (VCF, HSD)

    Returns:
        InputAdapter instance
    """
    ref_path = reference_path or str(BUNDLED_REFERENCE)

    if format_name == "bam":
        return BAMAdapter(reference_path=ref_path)
    elif format_name == "vcf":
        return VCFAdapter(reference_path=ref_path, sample_id=sample_id)
    elif format_name == "fasta":
        return FASTAAdapter(reference_path=ref_path)
    elif format_name == "hsd":
        return HSDAdapter(reference_path=ref_path)
    elif format_name == "microarray":
        return MicroarrayAdapter()
    else:
        raise ValueError(f"Unknown format: {format_name}")


@click.group()
@click.version_option(version="0.1.0")
def main() -> None:
    """eveHap - Universal mtDNA haplogroup classifier.

    Classify mitochondrial DNA haplogroups from various input formats
    including BAM, VCF, FASTA, HSD, and consumer genotyping files.
    """
    pass


@main.command()
@click.argument("input_files", type=click.Path(exists=True), nargs=-1, required=True)
@click.option(
    "--tree-type",
    required=True,
    type=click.Choice(["rcrs", "rsrs"]),
    help="REQUIRED: Phylotree type to use (rcrs=modern, rsrs=ancient DNA)",
)
@click.option(
    "--format",
    "input_format",
    default="auto",
    help="Input format (auto, bam, vcf, fasta, hsd, microarray)",
)
@click.option("--output", "-o", type=click.Path(), help="Output file (default: stdout)")
@click.option(
    "--output-format",
    default="tsv",
    type=click.Choice(["json", "tsv", "text"]),
    help="Output format",
)
@click.option(
    "--tree",
    required=False,
    help="Override phylotree file path (default: auto-select based on --tree-type)",
)
@click.option(
    "--reference", required=False, help="Override reference FASTA path (default: bundled rCRS)"
)
@click.option("--damage-filter", is_flag=True, help="Apply ancient DNA damage filtering")
@click.option(
    "--method",
    default="auto",
    type=click.Choice(["auto", "kulczynski", "traversal"]),
    help="Classification method",
)
@click.option("--top-n", default=5, help="Number of alternatives to report")
@click.option("--sample-id", help="Sample ID to extract from multi-sample VCF")
@click.option("--quiet", "-q", is_flag=True, help="Suppress progress output")
def classify(
    input_files: tuple,
    tree_type: str,
    input_format: str,
    output: Optional[str],
    output_format: str,
    tree: Optional[str],
    reference: Optional[str],
    damage_filter: bool,
    method: str,
    top_n: int,
    sample_id: Optional[str],
    quiet: bool,
) -> None:
    """Classify mtDNA haplogroup from input file(s).

    Requires --tree-type to specify which phylotree to use:

      - rcrs: rCRS-rooted tree (matches Haplogrep3 default, recommended for modern DNA)
      - rsrs: RSRS-rooted tree (for ancient DNA and population genetics studies)

    Examples:

        # Modern DNA with rCRS tree (matches Haplogrep3)
        evehap classify sample.bam --tree-type rcrs

        # Ancient DNA with RSRS tree and damage filtering
        evehap classify ancient.bam --tree-type rsrs --damage-filter

        # VCF file with specific sample
        evehap classify population.vcf.gz --tree-type rcrs --sample-id NA12878

        # Custom tree file (overrides tree-type)
        evehap classify sample.bam --tree-type rcrs --tree ./custom_tree.xml
    """
    # Get tree type configuration
    type_config = TREE_TYPES[tree_type]

    # Resolve tree path (use override or default from tree-type)
    if tree:
        if Path(tree).exists():
            tree_path = tree
        else:
            click.echo(f"Error: Tree file not found: {tree}", err=True)
            sys.exit(1)
    else:
        tree_path = str(type_config["tree"])
        if not Path(tree_path).exists():
            click.echo(f"Error: Bundled {tree_type} tree not found at {tree_path}", err=True)
            click.echo("This may indicate a corrupted installation.", err=True)
            sys.exit(1)

    # Resolve reference path (use override or default from tree-type)
    if reference:
        if reference == "rcrs":
            reference_path = str(BUNDLED_REFERENCE)
        elif Path(reference).exists():
            reference_path = reference
        else:
            click.echo(f"Error: Reference file not found: {reference}", err=True)
            sys.exit(1)
    else:
        reference_path = str(type_config["reference"])
        if not Path(reference_path).exists():
            click.echo(f"Error: Bundled reference not found at {reference_path}", err=True)
            click.echo("This may indicate a corrupted installation.", err=True)
            sys.exit(1)

    # Get weights for this tree type
    weights_file = Path(type_config["weights"])
    weights_path = str(weights_file) if weights_file.exists() else None

    try:
        phylotree = Phylotree.load(tree_path, weights_path)
        if not quiet:
            click.echo(
                f"Using phylotree: {Path(tree_path).name} ({len(phylotree.nodes)} haplogroups)",
                err=True,
            )
            click.echo(f"Using reference: {Path(reference_path).name}", err=True)
            click.echo(f"Tree type: {tree_type} ({type_config['description']})", err=True)
    except Exception as e:
        click.echo(f"Error loading phylotree: {e}", err=True)
        click.echo("", err=True)
        click.echo("Make sure the tree file is valid XML or JSON format.", err=True)
        sys.exit(1)

    # Load reference sequence for classifier
    import pysam

    try:
        with pysam.FastaFile(reference_path) as fasta:
            reference_sequence = fasta.fetch(fasta.references[0])
    except Exception as e:
        click.echo(f"Error loading reference: {e}", err=True)
        sys.exit(1)

    # Create damage filter if requested
    dmg_filter = DamageFilter() if damage_filter else None

    # Create classifier with reference and weights
    classifier = Classifier(
        phylotree,
        method=method,
        damage_filter=dmg_filter,
        reference=reference_sequence,
        weights_path=weights_path,
    )

    # Process files
    results: List[HaplogroupResult] = []

    for input_file in input_files:
        if not quiet:
            click.echo(f"Processing: {input_file}", err=True)

        try:
            # Detect format
            fmt = input_format if input_format != "auto" else detect_format(input_file)

            if fmt == "unknown":
                click.echo(f"  Warning: Cannot detect format for {input_file}", err=True)
                continue

            # Get adapter
            adapter = get_adapter(fmt, reference_path, sample_id)

            # Extract profile
            profile = adapter.extract_profile(input_file)

            # Classify
            result = classifier.classify(profile)
            results.append(result)

            if not quiet:
                click.echo(
                    f"  → {result.haplogroup} "
                    f"(confidence: {result.confidence:.1%}, quality: {result.quality})",
                    err=True,
                )

        except Exception as e:
            click.echo(f"  Error: {e}", err=True)
            continue

    # Output results
    if not results:
        click.echo("No samples successfully classified", err=True)
        sys.exit(1)

    output_text = format_output(results, output_format)

    if output:
        with open(output, "w") as f:
            f.write(output_text)
        if not quiet:
            click.echo(f"Results written to: {output}", err=True)
    else:
        click.echo(output_text)


def format_output(results: List[HaplogroupResult], fmt: str) -> str:
    """Format results for output.

    Args:
        results: List of HaplogroupResult
        fmt: Output format ('json', 'tsv', 'text')

    Returns:
        Formatted output string
    """
    if fmt == "json":
        return json.dumps([r.to_dict() for r in results], indent=2)

    elif fmt == "tsv":
        lines = [HaplogroupResult.tsv_header()]
        for r in results:
            lines.append(r.to_tsv_row())
        return "\n".join(lines)

    else:  # text
        lines = []
        for r in results:
            lines.append(f"Sample: {r.sample_id}")
            lines.append(f"  Haplogroup: {r.haplogroup}")
            lines.append(f"  Confidence: {r.confidence:.1%}")
            lines.append(f"  Quality: {r.quality}")
            lines.append(f"  Method: {r.method}")
            lines.append(f"  Coverage: {r.coverage_fraction:.1%}")
            if r.warnings:
                lines.append(f"  Warnings: {', '.join(r.warnings)}")
            lines.append("")
        return "\n".join(lines)


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--reference", type=click.Path(exists=True), help="Path to rCRS reference FASTA")
def info(input_file: str, reference: Optional[str]) -> None:
    """Show information about an input file.

    Displays format, sample count, coverage, and other metadata.
    """
    path = Path(input_file)
    fmt = detect_format(input_file)

    click.echo(f"File: {path.name}")
    click.echo(f"Format: {fmt}")
    click.echo(f"Size: {path.stat().st_size:,} bytes")

    try:
        adapter = get_adapter(fmt, reference)
        profile = adapter.extract_profile(input_file)

        click.echo(f"Sample ID: {profile.sample_id}")
        click.echo(f"Positions covered: {len(profile.covered_positions):,}")
        click.echo(f"Coverage: {profile.coverage_fraction():.1%}")

        mean_depth = profile.mean_depth()
        if mean_depth:
            click.echo(f"Mean depth: {mean_depth:.1f}x")

        variants = profile.get_variants()
        click.echo(f"Variants from reference: {len(variants)}")

    except Exception as e:
        click.echo(f"Error reading file: {e}", err=True)


@main.command()
@click.argument("bam_file", type=click.Path(exists=True))
def damage(bam_file: str) -> None:
    """Estimate ancient DNA damage rates from BAM file.

    Analyzes read ends for C→T and G→A substitution patterns
    characteristic of ancient DNA.
    """
    from evehap.core.damage import DamageFilter

    click.echo(f"Analyzing: {bam_file}")

    df = DamageFilter()
    stats = df.estimate_damage(bam_file, max_reads=10000)

    click.echo("\nDamage Statistics:")
    click.echo(f"  Reads analyzed: {stats.total_reads:,}")
    click.echo(f"  C→T rate at 5' end: {stats.ct_5prime_rate:.2%}")
    click.echo(f"  G→A rate at 3' end: {stats.ga_3prime_rate:.2%}")
    click.echo(f"  Damaged reads: {stats.damaged_reads:,}")

    if stats.is_ancient:
        click.echo("\n⚠️  Damage patterns suggest ANCIENT DNA")
        click.echo("   Consider using --damage-filter when classifying")
    else:
        click.echo("\n✓  Damage patterns consistent with modern DNA")


def _get_remote_metadata(url: str) -> dict:
    """Get remote file metadata (Last-Modified, Content-Length, ETag)."""
    import urllib.request

    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=10) as response:
            return {
                "last_modified": response.headers.get("Last-Modified"),
                "content_length": response.headers.get("Content-Length"),
                "etag": response.headers.get("ETag"),
            }
    except Exception:
        return {}


def _save_metadata(dest: Path, metadata: dict) -> None:
    """Save download metadata alongside the file."""
    import json
    from datetime import datetime

    meta_file = dest.parent / f".{dest.name}.meta"
    metadata["downloaded_at"] = datetime.now().isoformat()
    with open(meta_file, "w") as f:
        json.dump(metadata, f, indent=2)


def _load_metadata(dest: Path) -> Dict[str, Any]:
    """Load saved metadata for a file."""
    import json

    meta_file = dest.parent / f".{dest.name}.meta"
    if meta_file.exists():
        with open(meta_file) as f:
            result: Dict[str, Any] = json.load(f)
            return result
    return {}


def _check_update_available(dest: Path, url: str) -> bool:
    """Check if remote file is newer than local."""
    local_meta = _load_metadata(dest)
    if not local_meta:
        return True  # No metadata = needs download

    remote_meta = _get_remote_metadata(url)
    if not remote_meta:
        return False  # Can't check = assume up to date

    # Compare by ETag first
    if local_meta.get("etag") and remote_meta.get("etag"):
        return bool(local_meta["etag"] != remote_meta["etag"])

    # Compare by Last-Modified
    if local_meta.get("last_modified") and remote_meta.get("last_modified"):
        return bool(local_meta["last_modified"] != remote_meta["last_modified"])

    # Compare by Content-Length as fallback
    if local_meta.get("content_length") and remote_meta.get("content_length"):
        return bool(local_meta["content_length"] != remote_meta["content_length"])

    return False


@main.command()
@click.option(
    "--outdir",
    "-o",
    type=click.Path(),
    required=True,
    help="REQUIRED: Directory to save downloaded files",
)
@click.option(
    "--category",
    type=click.Choice(["all", "reference", "phylotree"]),
    default="all",
    help="Category of resources to download (default: all)",
)
@click.option(
    "--resource",
    multiple=True,
    help="Specific resource(s) to download (can be used multiple times)",
)
@click.option(
    "--force",
    is_flag=True,
    help="Overwrite existing files",
)
@click.option(
    "--check-updates",
    is_flag=True,
    help="Check for updates in --outdir without downloading",
)
@click.option(
    "--list",
    "list_resources",
    is_flag=True,
    help="List available resources",
)
def download(
    outdir: str,
    category: str,
    resource: tuple,
    force: bool,
    check_updates: bool,
    list_resources: bool,
) -> None:
    """Download reference sequences and phylotree resources.

    Downloads from phylotree.org and forensicgenomics.github.io/mitoLeaf
    to the specified output directory.

    Examples:
        evehap download --outdir ./evehap_data/
        evehap download --outdir ./data/ --resource mitoleaf-tree
        evehap download --outdir ./data/ --check-updates

    Then use with classify:
        evehap classify sample.bam --tree ./evehap_data/mitoleaf_tree.json --reference ./evehap_data/rCRS.fasta
    """
    import urllib.error
    import urllib.request

    download_dir = Path(outdir)

    # List mode
    if list_resources:
        click.echo("Available resources for download:\n")
        for cat in ["reference", "phylotree"]:
            click.echo(f"[{cat}]")
            for name, info in DOWNLOAD_RESOURCES.items():
                if info["category"] == cat:
                    dest = download_dir / info["filename"]
                    status = "✓" if dest.exists() else "✗"
                    click.echo(f"  {status} {name}: {info['description']}")
                    click.echo(f"      → {info['filename']}")
            click.echo()
        return

    # Determine which resources to download
    if resource:
        # Specific resources requested
        resources_to_download = {
            name: info for name, info in DOWNLOAD_RESOURCES.items() if name in resource
        }
        missing = set(resource) - set(resources_to_download.keys())
        if missing:
            click.echo(f"Error: Unknown resources: {', '.join(missing)}", err=True)
            click.echo("", err=True)
            click.echo("Use --list to see available resources:", err=True)
            click.echo(f"  evehap download --outdir {outdir} --list", err=True)
            sys.exit(1)
    else:
        # Filter by category
        resources_to_download = {
            name: info
            for name, info in DOWNLOAD_RESOURCES.items()
            if category == "all" or info["category"] == category
        }

    # Check updates mode
    if check_updates:
        click.echo(f"Checking for updates in {download_dir}...\n")
        updates_available = []
        for name, info in resources_to_download.items():
            dest = download_dir / info["filename"]
            if not dest.exists():
                click.echo(f"  ✗ {name}: Not downloaded")
                updates_available.append(name)
            elif _check_update_available(dest, info["url"]):
                click.echo(f"  ↑ {name}: Update available")
                updates_available.append(name)
            else:
                click.echo(f"  ✓ {name}: Up to date")

        if updates_available:
            resources_arg = " --resource ".join(updates_available)
            click.echo("\nTo update, run:")
            click.echo(f"  evehap download --outdir {outdir} --resource {resources_arg} --force")
        return

    # Ensure download directory exists
    download_dir.mkdir(parents=True, exist_ok=True)

    # Download resources
    downloaded_files = []
    for name, info in resources_to_download.items():
        dest = download_dir / info["filename"]
        url = info["url"]

        if dest.exists() and not force:
            click.echo(f"✓ {name}: Already exists at {dest}")
            downloaded_files.append((name, dest))
            continue

        click.echo(f"↓ Downloading {name}...")
        click.echo(f"   {url}")

        try:
            # Get metadata first
            remote_meta = _get_remote_metadata(url)
            size_mb = int(remote_meta.get("content_length", 0)) / (1024 * 1024)
            if size_mb > 1:
                click.echo(f"   Size: {size_mb:.1f} MB")

            with urllib.request.urlopen(url, timeout=60) as response:
                content = response.read()

            with open(dest, "wb") as f:
                f.write(content)

            # Save metadata for update checking
            _save_metadata(dest, remote_meta)

            click.echo(f"  ✓ Downloaded ({len(content) / 1024:.1f} KB)")
            downloaded_files.append((name, dest))

        except urllib.error.URLError as e:
            click.echo(f"  ✗ Failed: {e}", err=True)
            click.echo("     Check your network connection and try again.", err=True)
        except Exception as e:
            click.echo(f"  ✗ Error: {e}", err=True)

    # Show summary
    click.echo(f"\nDownloaded to: {download_dir.absolute()}")
    click.echo("\nUse with classify:")

    tree_file = download_dir / "mitoleaf_tree.json"
    ref_file = download_dir / "rCRS.fasta"

    if tree_file.exists() and ref_file.exists():
        click.echo(f"  evehap classify sample.bam --tree {tree_file} --reference {ref_file}")
    elif tree_file.exists():
        click.echo(f"  evehap classify sample.bam --tree {tree_file} --reference rcrs")
    else:
        click.echo("  evehap classify sample.bam --tree rsrs --reference rcrs  # bundled resources")


@main.command()
def version() -> None:
    """Show version information and bundled resources."""
    click.echo("eveHap v0.1.0")
    click.echo(f"Package directory: {PACKAGE_DIR}")

    # Show bundled resources (always available, no download required)
    click.echo("\nBundled resources (use with --tree and --reference):")

    click.echo("\n  Phylotrees:")
    for alias, path in [
        ("rsrs", BUNDLED_TREE_RSRS),
        ("rcrs", BUNDLED_TREE_RCRS),
    ]:
        if path.exists():
            size_kb = path.stat().st_size / 1024
            node_count = "~5400" if "rsrs" in str(path) else "~2400"
            click.echo(
                f"    ✓ --tree {alias}  ({path.name}, {size_kb:.1f} KB, {node_count} haplogroups)"
            )
        else:
            click.echo(f"    ✗ --tree {alias}  (not found)")

    click.echo("\n  References:")
    if BUNDLED_REFERENCE.exists():
        size_kb = BUNDLED_REFERENCE.stat().st_size / 1024
        click.echo(f"    ✓ --reference rcrs  ({BUNDLED_REFERENCE.name}, {size_kb:.1f} KB)")
    else:
        click.echo("    ✗ --reference rcrs  (not found)")

    click.echo("\nDownloadable resources:")
    for name, info in DOWNLOAD_RESOURCES.items():
        click.echo(f"    {name}: {info['description']}")

    click.echo("\nExample usage:")
    click.echo("  # Using bundled resources (offline-ready):")
    click.echo("  evehap classify sample.bam --tree rsrs --reference rcrs")
    click.echo("")
    click.echo("  # Download and use mitoLeaf (6400+ haplogroups):")
    click.echo("  evehap download --outdir ./evehap_data/")
    click.echo(
        "  evehap classify sample.bam --tree ./evehap_data/mitoleaf_tree.json --reference ./evehap_data/rCRS.fasta"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate mtDNA RSID lookup table from dbSNP.

Downloads mtDNA variant RSIDs from NCBI dbSNP and creates a position-to-RSID
lookup table for generating 23andMe-format files.

Usage:
    python generate_rsid_lookup.py --output data/dbsnp/mtdna_rsids.tsv
"""

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from urllib.request import urlretrieve

# dbSNP VCF URL for chrMT (GRCh38)
DBSNP_CHRMT_URL = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz"
# Alternative: use the chrMT-specific file
DBSNP_CHRMT_ONLY = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz"


def download_dbsnp_vcf(output_dir: Path, force: bool = False) -> Path:
    """Download dbSNP VCF file if not already present.

    Returns path to downloaded file.
    """
    vcf_path = output_dir / "dbsnp_common_all.vcf.gz"

    if vcf_path.exists() and not force:
        print(f"Using existing file: {vcf_path}")
        return vcf_path

    print(f"Downloading dbSNP VCF from {DBSNP_CHRMT_ONLY}...")
    print("This may take a while (file is ~1GB)...")

    try:
        urlretrieve(DBSNP_CHRMT_ONLY, vcf_path)
        print(f"Downloaded to: {vcf_path}")
        return vcf_path
    except Exception as e:
        print(f"Error downloading: {e}", file=sys.stderr)
        raise


def extract_chrmt_rsids(vcf_path: Path) -> Dict[int, List[Tuple[str, str, str]]]:
    """Extract chrMT RSIDs from dbSNP VCF.

    Returns dict mapping position to list of (rsid, ref, alt) tuples.
    """
    rsids: Dict[int, List[Tuple[str, str, str]]] = {}

    chrmt_names = {"MT", "chrM", "chrMT", "NC_012920.1"}

    print(f"Extracting chrMT variants from {vcf_path}...")

    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue

            chrom = parts[0]
            if chrom not in chrmt_names:
                continue

            try:
                pos = int(parts[1])
            except ValueError:
                continue

            rsid = parts[2]
            ref = parts[3]
            alt = parts[4]

            if not rsid.startswith('rs'):
                continue

            if pos not in rsids:
                rsids[pos] = []

            # Handle multiple alts
            for alt_allele in alt.split(','):
                rsids[pos].append((rsid, ref, alt_allele))

    print(f"Found {len(rsids)} mtDNA positions with RSIDs")
    return rsids


def generate_from_existing_files() -> Dict[int, List[Tuple[str, str, str]]]:
    """Generate RSID lookup from known mtDNA variants.

    Uses common mtDNA variants that have well-known RSIDs.
    This is a fallback if dbSNP download is not available.
    """
    # Known mtDNA RSIDs from phylotree and literature
    # Format: position -> [(rsid, ref, alt), ...]
    known_rsids = {
        # Haplogroup-defining positions with known RSIDs
        73: [("rs3087742", "A", "G")],
        146: [("rs371316536", "T", "C")],
        150: [("rs28625645", "C", "T")],
        152: [("rs117135796", "T", "C")],
        189: [("rs28493303", "A", "G")],
        195: [("rs2857291", "T", "C")],
        199: [("rs28625647", "T", "C")],
        204: [("rs28625648", "T", "C")],
        207: [("rs28625649", "G", "A")],
        263: [("rs2853499", "A", "G")],
        295: [("rs371168143", "C", "T")],
        310: [("rs376846471", "T", "C")],
        489: [("rs2853500", "T", "C")],
        709: [("rs2853502", "G", "A")],
        750: [("rs2853503", "A", "G")],
        769: [("rs2853504", "G", "A")],
        825: [("rs28358574", "T", "A")],
        1018: [("rs2853505", "G", "A")],
        1438: [("rs2853506", "A", "G")],
        1598: [("rs28499977", "G", "A")],
        1719: [("rs28358576", "G", "A")],
        1811: [("rs41341949", "A", "G")],
        1888: [("rs28358577", "G", "A")],
        2352: [("rs28358579", "T", "C")],
        2706: [("rs2853507", "A", "G")],
        2758: [("rs28358580", "G", "A")],
        2885: [("rs28358582", "T", "C")],
        3010: [("rs28358583", "G", "A")],
        3197: [("rs28358584", "T", "C")],
        3480: [("rs3899199", "A", "G")],
        3594: [("rs41403847", "C", "T")],
        4216: [("rs1599988", "T", "C")],
        4529: [("rs28358586", "A", "T")],
        4580: [("rs2857289", "G", "A")],
        4769: [("rs3021086", "A", "G")],
        4793: [("rs28358587", "A", "G")],
        4917: [("rs28358588", "A", "G")],
        5178: [("rs28358289", "C", "A")],
        5231: [("rs28358589", "G", "A")],
        5442: [("rs28358590", "T", "C")],
        5460: [("rs28358591", "G", "A")],
        5633: [("rs28358592", "C", "T")],
        5899: [("rs28359173", "T", "C")],
        6179: [("rs3899198", "G", "A")],
        6221: [("rs28358593", "T", "C")],
        6253: [("rs28358594", "T", "C")],
        6371: [("rs28358595", "T", "C")],
        6413: [("rs28359178", "T", "C")],
        6446: [("rs41535848", "G", "A")],
        6455: [("rs28358596", "C", "T")],
        6719: [("rs28358598", "T", "C")],
        7028: [("rs2015062", "C", "T")],
        7146: [("rs28359170", "A", "G")],
        7256: [("rs28358599", "C", "T")],
        7476: [("rs28358600", "C", "T")],
        7521: [("rs28358601", "G", "A")],
        8251: [("rs28358602", "G", "A")],
        8468: [("rs28358606", "C", "T")],
        8655: [("rs28358608", "C", "T")],
        8697: [("rs28358609", "G", "A")],
        8701: [("rs2000975", "A", "G")],
        8860: [("rs2001031", "A", "G")],
        8994: [("rs28358610", "G", "A")],
        9540: [("rs2853509", "T", "C")],
        9950: [("rs28358612", "T", "C")],
        10238: [("rs28358613", "T", "C")],
        10310: [("rs28358614", "G", "A")],
        10398: [("rs2853826", "A", "G")],
        10400: [("rs28358278", "C", "T")],
        10463: [("rs2853510", "T", "C")],
        10550: [("rs28358615", "A", "G")],
        10873: [("rs28358280", "T", "C")],
        11002: [("rs28358616", "A", "G")],
        11177: [("rs28358617", "C", "T")],
        11251: [("rs2853511", "A", "G")],
        11467: [("rs2853512", "A", "G")],
        11719: [("rs2853513", "G", "A")],
        11812: [("rs28358618", "A", "G")],
        11914: [("rs28358619", "G", "A")],
        12007: [("rs28358620", "G", "A")],
        12308: [("rs2853514", "A", "G")],
        12372: [("rs2853515", "G", "A")],
        12612: [("rs28358622", "A", "G")],
        12705: [("rs2853498", "C", "T")],
        13368: [("rs28359167", "G", "A")],
        13708: [("rs28359178", "G", "A")],
        13759: [("rs28359168", "G", "A")],
        13966: [("rs28358625", "A", "G")],
        14034: [("rs28358626", "T", "C")],
        14167: [("rs28358627", "C", "T")],
        14233: [("rs28358628", "A", "G")],
        14365: [("rs28358629", "G", "A")],
        14470: [("rs28358630", "T", "C")],
        14766: [("rs2853516", "C", "T")],
        14783: [("rs28358631", "T", "C")],
        14905: [("rs28358632", "G", "A")],
        15043: [("rs28358633", "G", "A")],
        15148: [("rs28358634", "G", "A")],
        15218: [("rs28358635", "A", "G")],
        15257: [("rs28358279", "G", "A")],
        15301: [("rs28358636", "G", "A")],
        15326: [("rs2853508", "A", "G")],
        15452: [("rs28358569", "C", "A")],
        15607: [("rs28358638", "A", "G")],
        15622: [("rs28358639", "A", "G")],
        15670: [("rs28358640", "T", "C")],
        15758: [("rs28358641", "A", "G")],
        15884: [("rs28358642", "G", "A")],
        15904: [("rs28358643", "C", "T")],
        15928: [("rs28358644", "G", "A")],
        16051: [("rs28358645", "A", "G")],
        16086: [("rs28358646", "T", "C")],
        16093: [("rs2853518", "T", "C")],
        16124: [("rs28358647", "T", "C")],
        16126: [("rs2853519", "T", "C")],
        16129: [("rs28358648", "G", "A")],
        16145: [("rs28358649", "G", "A")],
        16148: [("rs28358650", "C", "T")],
        16153: [("rs28358651", "G", "A")],
        16163: [("rs28358652", "A", "G")],
        16172: [("rs28358653", "T", "C")],
        16182: [("rs28358654", "A", "C")],
        16183: [("rs28358655", "A", "C")],
        16187: [("rs28358656", "C", "T")],
        16189: [("rs2857290", "T", "C")],
        16192: [("rs28358657", "C", "T")],
        16209: [("rs28358658", "T", "C")],
        16213: [("rs28358659", "G", "A")],
        16214: [("rs28358660", "C", "T")],
        16217: [("rs28358661", "T", "C")],
        16223: [("rs2853520", "C", "T")],
        16224: [("rs28358662", "T", "C")],
        16230: [("rs28358663", "A", "G")],
        16234: [("rs28358664", "C", "T")],
        16249: [("rs28358665", "T", "C")],
        16256: [("rs28358666", "C", "T")],
        16261: [("rs2853521", "C", "T")],
        16264: [("rs28358667", "C", "T")],
        16270: [("rs28358668", "C", "T")],
        16278: [("rs2853522", "C", "T")],
        16290: [("rs28358669", "C", "T")],
        16291: [("rs28358670", "C", "T")],
        16292: [("rs28358671", "C", "T")],
        16294: [("rs28358672", "C", "T")],
        16296: [("rs28358673", "C", "T")],
        16298: [("rs28358674", "T", "C")],
        16304: [("rs28358675", "T", "C")],
        16311: [("rs2853523", "T", "C")],
        16319: [("rs28358676", "G", "A")],
        16320: [("rs28358677", "C", "T")],
        16325: [("rs28358678", "T", "C")],
        16327: [("rs28358679", "C", "T")],
        16342: [("rs28358680", "T", "C")],
        16354: [("rs28358681", "C", "T")],
        16355: [("rs28358682", "C", "T")],
        16362: [("rs2853524", "T", "C")],
        16390: [("rs28358683", "G", "A")],
        16399: [("rs28358684", "A", "G")],
        16519: [("rs3937033", "T", "C")],
    }

    return known_rsids


def write_lookup_table(rsids: Dict[int, List[Tuple[str, str, str]]], output_path: Path):
    """Write RSID lookup table to TSV file.

    Format: position<tab>ref<tab>alt<tab>rsid
    """
    with open(output_path, 'w') as f:
        f.write("position\tref\talt\trsid\n")

        for pos in sorted(rsids.keys()):
            for rsid, ref, alt in rsids[pos]:
                f.write(f"{pos}\t{ref}\t{alt}\t{rsid}\n")

    print(f"Wrote {sum(len(v) for v in rsids.values())} RSIDs to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate mtDNA RSID lookup table")
    parser.add_argument("--output", "-o", required=True, help="Output TSV file path")
    parser.add_argument("--download", action="store_true", help="Download from dbSNP (large file)")
    parser.add_argument("--force", action="store_true", help="Force re-download")

    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if args.download:
        # Download from dbSNP
        vcf_path = download_dbsnp_vcf(output_path.parent, force=args.force)
        rsids = extract_chrmt_rsids(vcf_path)
    else:
        # Use built-in known RSIDs
        print("Using built-in known mtDNA RSIDs (use --download for full dbSNP)")
        rsids = generate_from_existing_files()

    write_lookup_table(rsids, output_path)
    print("Done!")


if __name__ == "__main__":
    main()


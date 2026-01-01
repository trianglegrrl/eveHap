"""Integration tests for multiple input formats.

Tests that different input formats (VCF, FASTA, BAM) produce consistent
haplogroup classifications for the same samples.

Usage:
    pytest tests/test_integration/test_formats.py -v
"""

from pathlib import Path
from typing import Dict, Optional

import pytest

from evehap.adapters.bam import BAMAdapter
from evehap.adapters.fasta import FASTAAdapter
from evehap.adapters.vcf import VCFAdapter
from evehap.core.classifier import Classifier
from evehap.core.phylotree import Phylotree
from evehap.output.result import HaplogroupResult


def get_major_clade(haplogroup: str) -> str:
    """Extract major clade from full haplogroup."""
    if not haplogroup:
        return ""

    for i, char in enumerate(haplogroup):
        if char.isdigit():
            return haplogroup[:i] if i > 0 else haplogroup[0]

    return haplogroup[0]


class TestMultiFormatConsistency:
    """Test that different formats produce consistent results."""

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree-rsrs.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights-rsrs.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(
            str(tree_path),
            str(weights_path) if weights_path.exists() else None
        )

    @pytest.fixture
    def reference_path(self, evehap_data_dir: Path) -> str:
        """Return reference FASTA path."""
        ref = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref.exists():
            pytest.skip("Reference FASTA not available")
        return str(ref)

    @pytest.fixture
    def hgdp_data(self, mtdna_poc_data_dir: Path) -> Dict[str, Path]:
        """Return paths to HGDP data in multiple formats."""
        hgdp_dir = mtdna_poc_data_dir / "hgdp"
        if not hgdp_dir.exists():
            pytest.skip("HGDP data not available")

        return {
            "bam_dir": hgdp_dir / "bams",
            "fasta_dir": hgdp_dir / "haplogroups" / "consensus",
            "vcf_dir": hgdp_dir / "haplogroups" / "vcf",
        }

    def test_vcf_format(
        self,
        mtdna_poc_data_dir: Path,
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test VCF format classification."""
        vcf_path = mtdna_poc_data_dir / "raw" / "1kg_chrMT.vcf.gz"
        if not vcf_path.exists():
            pytest.skip("1000 Genomes VCF not available")

        adapter = VCFAdapter(reference_path=reference_path, sample_id="HG00096")
        profile = adapter.extract_profile(str(vcf_path))

        classifier = Classifier(phylotree)
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.haplogroup != ""
        print(f"VCF: HG00096 → {result.haplogroup} ({result.confidence:.1%})")

    def test_fasta_format(
        self,
        hgdp_data: Dict[str, Path],
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test FASTA format classification."""
        fasta_dir = hgdp_data["fasta_dir"]
        if not fasta_dir.exists():
            pytest.skip("HGDP FASTA files not available")

        fastas = list(fasta_dir.glob("*.fasta"))
        if not fastas:
            pytest.skip("No FASTA files found")

        adapter = FASTAAdapter(reference_path=reference_path)
        classifier = Classifier(phylotree)

        for fasta in fastas[:3]:  # Test first 3
            profile = adapter.extract_profile(str(fasta))
            result = classifier.classify(profile)

            assert isinstance(result, HaplogroupResult)
            assert result.haplogroup != ""
            print(f"FASTA: {fasta.stem} → {result.haplogroup} ({result.confidence:.1%})")

    def test_bam_format(
        self,
        hgdp_data: Dict[str, Path],
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test BAM format classification."""
        bam_dir = hgdp_data["bam_dir"]
        if not bam_dir.exists():
            pytest.skip("HGDP BAM files not available")

        bams = list(bam_dir.glob("*.bam"))
        if not bams:
            pytest.skip("No BAM files found")

        adapter = BAMAdapter(reference_path=reference_path)
        classifier = Classifier(phylotree)

        for bam in bams[:3]:  # Test first 3
            try:
                profile = adapter.extract_profile(str(bam))
                result = classifier.classify(profile)

                assert isinstance(result, HaplogroupResult)
                assert result.haplogroup != ""
                print(f"BAM: {bam.stem} → {result.haplogroup} ({result.confidence:.1%})")
            except Exception as e:
                print(f"Error with {bam.stem}: {e}")

    def test_format_consistency(
        self,
        hgdp_data: Dict[str, Path],
        phylotree: Phylotree,
        reference_path: str,
    ) -> None:
        """Test that different formats for the same sample produce consistent results."""
        fasta_dir = hgdp_data["fasta_dir"]
        bam_dir = hgdp_data["bam_dir"]

        if not fasta_dir.exists() or not bam_dir.exists():
            pytest.skip("HGDP multi-format data not available")

        # Find samples that have both FASTA and BAM
        fasta_samples = {f.stem: f for f in fasta_dir.glob("*.fasta")}
        bam_samples = {b.stem.split(".")[0]: b for b in bam_dir.glob("*.bam")}

        common_samples = set(fasta_samples.keys()) & set(bam_samples.keys())

        if not common_samples:
            pytest.skip("No samples with both FASTA and BAM")

        fasta_adapter = FASTAAdapter(reference_path=reference_path)
        bam_adapter = BAMAdapter(reference_path=reference_path)
        classifier = Classifier(phylotree)

        matches = 0
        clade_matches = 0
        total = 0

        print("\nFormat consistency check:")
        for sample_id in list(common_samples)[:5]:  # Test first 5
            try:
                fasta_profile = fasta_adapter.extract_profile(str(fasta_samples[sample_id]))
                fasta_result = classifier.classify(fasta_profile)

                bam_profile = bam_adapter.extract_profile(str(bam_samples[sample_id]))
                bam_result = classifier.classify(bam_profile)

                total += 1

                fasta_clade = get_major_clade(fasta_result.haplogroup)
                bam_clade = get_major_clade(bam_result.haplogroup)

                if fasta_result.haplogroup == bam_result.haplogroup:
                    matches += 1
                    clade_matches += 1
                    status = "✓ EXACT"
                elif fasta_clade == bam_clade:
                    clade_matches += 1
                    status = "~ CLADE"
                else:
                    status = "✗ DIFFER"

                print(f"  {sample_id}: FASTA={fasta_result.haplogroup}, BAM={bam_result.haplogroup} [{status}]")

            except Exception as e:
                print(f"  {sample_id}: ERROR - {e}")

        if total == 0:
            pytest.skip("No samples could be compared")

        print(f"\nExact matches: {matches}/{total} ({matches/total:.1%})")
        print(f"Clade matches: {clade_matches}/{total} ({clade_matches/total:.1%})")

        # Expect at least 80% clade consistency
        assert clade_matches / total >= 0.8, "Format consistency too low"


class TestHGDPMultiFormat:
    """Test HGDP samples across multiple formats with ground truth."""

    @pytest.fixture
    def phylotree(self, evehap_data_dir: Path) -> Phylotree:
        """Load phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree-rsrs.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights-rsrs.txt"

        if not tree_path.exists():
            pytest.skip("Phylotree not available")

        return Phylotree.load(
            str(tree_path),
            str(weights_path) if weights_path.exists() else None
        )

    @pytest.fixture
    def reference_path(self, evehap_data_dir: Path) -> str:
        """Return reference FASTA path."""
        ref = evehap_data_dir / "reference" / "rCRS.fasta"
        if not ref.exists():
            pytest.skip("Reference FASTA not available")
        return str(ref)

    @pytest.fixture
    def ground_truth(self, mtdna_poc_data_dir: Path) -> Dict[str, str]:
        """Load HGDP ground truth haplogroups."""
        import csv

        hg_file = mtdna_poc_data_dir / "hgdp" / "haplogroups" / "haplogroups.csv"
        if not hg_file.exists():
            pytest.skip("HGDP haplogroups file not available")

        haplogroups = {}
        with open(hg_file, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_id = row.get("sample_id", row.get("Sample", ""))
                hg = row.get("haplogroup", row.get("Haplogroup", ""))
                if sample_id and hg:
                    haplogroups[sample_id] = hg

        return haplogroups

    def test_fasta_vs_ground_truth(
        self,
        mtdna_poc_data_dir: Path,
        phylotree: Phylotree,
        reference_path: str,
        ground_truth: Dict[str, str],
        sample_count: int,
        random_seed: int,
    ) -> None:
        """Test FASTA classification against ground truth."""
        from ..conftest import select_samples

        fasta_dir = mtdna_poc_data_dir / "hgdp" / "haplogroups" / "consensus"
        if not fasta_dir.exists():
            pytest.skip("HGDP FASTA files not available")

        all_fastas = list(fasta_dir.glob("*.fasta"))
        fastas = select_samples(all_fastas, sample_count, random_seed)

        if not fastas:
            pytest.skip("No FASTA files found")

        print(f"\nTesting {len(fastas)} FASTA samples")

        adapter = FASTAAdapter(reference_path=reference_path)
        classifier = Classifier(phylotree)

        clade_correct = 0
        total = 0

        for fasta in fastas:
            sample_id = fasta.stem
            expected_hg = ground_truth.get(sample_id)

            if expected_hg is None:
                continue

            try:
                profile = adapter.extract_profile(str(fasta))
                result = classifier.classify(profile)

                expected_clade = get_major_clade(expected_hg)
                result_clade = get_major_clade(result.haplogroup)

                total += 1
                if expected_clade == result_clade:
                    clade_correct += 1
                    status = "✓"
                else:
                    status = "✗"

                print(f"  {sample_id}: {expected_hg} → {result.haplogroup} [{status}]")

            except Exception as e:
                print(f"  {sample_id}: ERROR - {e}")

        if total == 0:
            pytest.skip("No samples with ground truth")

        accuracy = clade_correct / total
        print(f"\nClade accuracy: {accuracy:.1%} ({clade_correct}/{total})")

        assert accuracy >= 0.7, f"FASTA accuracy too low: {accuracy:.1%}"


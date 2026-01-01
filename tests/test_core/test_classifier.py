"""Tests for evehap.core.classifier module."""

from pathlib import Path
from typing import Optional

import pytest

from evehap.core.classifier import (
    Classifier,
    KulczynskiScorer,
    TreeTraverser,
    TraversalResult,
)
from evehap.core.phylotree import Phylotree, PhyloNode, Mutation
from evehap.core.profile import AlleleProfile, AlleleObservation
from evehap.output.result import HaplogroupResult


class TestKulczynskiScorer:
    """Tests for Kulczynski similarity scoring."""

    @pytest.fixture
    def simple_tree(self) -> Phylotree:
        """Create a simple phylotree for testing."""
        tree = Phylotree()

        # Create root
        root = PhyloNode(haplogroup="root")
        tree.root = root
        tree.nodes["root"] = root

        # Create L (African)
        L = PhyloNode(
            haplogroup="L",
            parent=root,
            defining_mutations=[Mutation(position=73, ancestral="A", derived="G")]
        )
        root.children.append(L)
        tree.nodes["L"] = L

        # Create N (non-African)
        N = PhyloNode(
            haplogroup="N",
            parent=root,
            defining_mutations=[
                Mutation(position=73, ancestral="A", derived="G"),
                Mutation(position=263, ancestral="A", derived="G"),
            ]
        )
        root.children.append(N)
        tree.nodes["N"] = N

        # Create H (European)
        H = PhyloNode(
            haplogroup="H",
            parent=N,
            defining_mutations=[
                Mutation(position=750, ancestral="A", derived="G"),
            ]
        )
        N.children.append(H)
        tree.nodes["H"] = H

        return tree

    def test_scorer_init(self, simple_tree: Phylotree) -> None:
        """Test scorer initialization."""
        scorer = KulczynskiScorer(simple_tree)
        assert scorer.phylotree is simple_tree

    def test_score_perfect_match(self, simple_tree: Phylotree) -> None:
        """Test scoring when profile perfectly matches haplogroup."""
        # Create profile with all H-defining mutations
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=750, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        scorer = KulczynskiScorer(simple_tree)
        h_node = simple_tree.get_haplogroup("H")
        score = scorer.score(profile, h_node)

        # Should be high score (close to 1.0)
        assert score > 0.8

    def test_score_no_match(self, simple_tree: Phylotree) -> None:
        """Test scoring when profile doesn't match haplogroup."""
        # Create profile with reference alleles (no mutations)
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"A": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"A": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=750, ref_allele="A", alleles={"A": 1.0}, source="test"
        ))

        scorer = KulczynskiScorer(simple_tree)
        h_node = simple_tree.get_haplogroup("H")
        score = scorer.score(profile, h_node)

        # Should be low score
        assert score < 0.3

    def test_score_partial_match(self, simple_tree: Phylotree) -> None:
        """Test scoring when profile partially matches haplogroup."""
        # Create profile with some H mutations but not all
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"A": 1.0}, source="test"  # Missing
        ))
        profile.add_observation(AlleleObservation(
            position=750, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        scorer = KulczynskiScorer(simple_tree)
        h_node = simple_tree.get_haplogroup("H")
        score = scorer.score(profile, h_node)

        # Should be medium score
        assert 0.3 < score < 0.9

    def test_score_all(self, simple_tree: Phylotree) -> None:
        """Test scoring all haplogroups."""
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=750, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        scorer = KulczynskiScorer(simple_tree)
        results = scorer.score_all(profile, top_n=3)

        assert len(results) > 0
        assert isinstance(results, list)
        assert all(isinstance(r, tuple) and len(r) == 2 for r in results)
        # Should be sorted by score descending
        scores = [r[1] for r in results]
        assert scores == sorted(scores, reverse=True)


class TestTreeTraverser:
    """Tests for tree traversal classification."""

    @pytest.fixture
    def simple_tree(self) -> Phylotree:
        """Create a simple phylotree for testing."""
        tree = Phylotree()

        root = PhyloNode(haplogroup="root")
        tree.root = root
        tree.nodes["root"] = root

        # L branch
        L = PhyloNode(
            haplogroup="L",
            parent=root,
            defining_mutations=[Mutation(position=73, ancestral="A", derived="G")]
        )
        root.children.append(L)
        tree.nodes["L"] = L

        # L1 sub-branch
        L1 = PhyloNode(
            haplogroup="L1",
            parent=L,
            defining_mutations=[Mutation(position=263, ancestral="A", derived="G")]
        )
        L.children.append(L1)
        tree.nodes["L1"] = L1

        # M branch
        M = PhyloNode(
            haplogroup="M",
            parent=root,
            defining_mutations=[
                Mutation(position=73, ancestral="A", derived="G"),
                Mutation(position=750, ancestral="A", derived="G"),
            ]
        )
        root.children.append(M)
        tree.nodes["M"] = M

        return tree

    def test_traverser_init(self, simple_tree: Phylotree) -> None:
        """Test traverser initialization."""
        traverser = TreeTraverser(simple_tree)
        assert traverser.phylotree is simple_tree
        assert traverser.max_conflicts == 3
        assert traverser.min_support_ratio == 0.7

    def test_traverse_high_support(self, simple_tree: Phylotree) -> None:
        """Test traversal with clear support for a branch."""
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        traverser = TreeTraverser(simple_tree)
        result = traverser.traverse(profile)

        assert isinstance(result, TraversalResult)
        assert result.haplogroup in ("L", "L1")  # Should traverse L branch
        assert len(result.path) > 0

    def test_traverse_stops_on_conflict(self, simple_tree: Phylotree) -> None:
        """Test traversal stops when conflicts exceed threshold."""
        # Profile with conflicting mutations
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"A": 1.0}, source="test"  # Ancestral
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=750, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        traverser = TreeTraverser(simple_tree, max_conflicts=0)
        result = traverser.traverse(profile)

        # Should stop early or return root
        assert result.stopped_early or result.haplogroup == "root"

    def test_traverse_sparse_data(self, simple_tree: Phylotree) -> None:
        """Test traversal with sparse coverage."""
        profile = AlleleProfile(sample_id="test", source_format="test")
        # Only cover position 73
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        traverser = TreeTraverser(simple_tree)
        result = traverser.traverse(profile)

        # Should still return a result
        assert isinstance(result, TraversalResult)
        assert result.uncovered_count >= 0


class TestClassifier:
    """Tests for unified Classifier."""

    @pytest.fixture
    def phylotree_path(self, evehap_data_dir: Path) -> Path:
        """Return path to phylotree."""
        tree_path = evehap_data_dir / "phylotree" / "tree.xml"
        if not tree_path.exists():
            pytest.skip("Phylotree not available")
        return tree_path

    @pytest.fixture
    def loaded_phylotree(self, phylotree_path: Path) -> Phylotree:
        """Load the phylotree."""
        weights_path = phylotree_path.parent / "weights.txt"
        return Phylotree.load(
            str(phylotree_path),
            str(weights_path) if weights_path.exists() else None
        )

    def test_classifier_init(self, loaded_phylotree: Phylotree) -> None:
        """Test classifier initialization."""
        classifier = Classifier(loaded_phylotree)
        assert classifier.phylotree is loaded_phylotree
        assert classifier.method == "auto"

    def test_classify_high_coverage(self, loaded_phylotree: Phylotree) -> None:
        """Test classification with high coverage profile."""
        # Create a profile with many positions covered
        profile = AlleleProfile(sample_id="test", source_format="test")

        # Add observations for all positions 1-16569 (high coverage)
        # For simplicity, fill with reference and a few key mutations
        for pos in range(1, 1000):  # First 1000 positions
            profile.add_observation(AlleleObservation(
                position=pos, ref_allele="A", alleles={"A": 1.0}, source="test"
            ))

        # Add H-typical mutations
        for pos in [73, 263, 750, 1438, 4769, 8860, 15326]:
            if pos <= 1000:
                profile.add_observation(AlleleObservation(
                    position=pos, ref_allele="A", alleles={"G": 1.0}, source="test"
                ))

        classifier = Classifier(loaded_phylotree, method="kulczynski")
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.sample_id == "test"
        assert result.haplogroup != ""

    def test_classify_low_coverage(self, loaded_phylotree: Phylotree) -> None:
        """Test classification with low coverage profile."""
        profile = AlleleProfile(sample_id="test", source_format="test")

        # Only a few positions covered
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))
        profile.add_observation(AlleleObservation(
            position=263, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        classifier = Classifier(loaded_phylotree, method="traversal")
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.method == "traversal"

    def test_classify_auto_method(self, loaded_phylotree: Phylotree) -> None:
        """Test auto method selection."""
        profile = AlleleProfile(sample_id="test", source_format="test")
        profile.add_observation(AlleleObservation(
            position=73, ref_allele="A", alleles={"G": 1.0}, source="test"
        ))

        classifier = Classifier(loaded_phylotree, method="auto")
        result = classifier.classify(profile)

        assert isinstance(result, HaplogroupResult)
        assert result.method in ("kulczynski", "traversal")



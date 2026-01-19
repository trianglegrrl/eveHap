"""Tests for evehap.core.phylotree module."""

from pathlib import Path

import pytest

from evehap.core.phylotree import Mutation, PhyloNode, Phylotree
from evehap.core.profile import AlleleObservation, AlleleProfile


class TestMutation:
    """Tests for Mutation dataclass."""

    def test_create_mutation(self) -> None:
        """Test creating a Mutation."""
        mut = Mutation(
            position=73,
            ancestral="A",
            derived="G",
            weight=1.5,
            is_backmutation=False,
        )
        assert mut.position == 73
        assert mut.ancestral == "A"
        assert mut.derived == "G"
        assert mut.weight == 1.5
        assert mut.is_backmutation is False

    def test_default_values(self) -> None:
        """Test default values."""
        mut = Mutation(position=73, ancestral="A", derived="G")
        assert mut.weight == 1.0
        assert mut.is_backmutation is False

    def test_repr(self) -> None:
        """Test string representation."""
        mut = Mutation(position=73, ancestral="A", derived="G")
        assert str(mut) == "73G"

    def test_matches_present(self) -> None:
        """Test matches returns True when mutation is present."""
        mut = Mutation(position=73, ancestral="A", derived="G")
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.add_observation(AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0}))
        assert mut.matches(profile) is True

    def test_matches_absent(self) -> None:
        """Test matches returns False when mutation is absent."""
        mut = Mutation(position=73, ancestral="A", derived="G")
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.add_observation(AlleleObservation(position=73, ref_allele="A", alleles={"A": 1.0}))
        assert mut.matches(profile) is False

    def test_matches_not_covered(self) -> None:
        """Test matches returns None when position not covered."""
        mut = Mutation(position=73, ancestral="A", derived="G")
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert mut.matches(profile) is None


class TestPhyloNode:
    """Tests for PhyloNode dataclass."""

    def test_create_node(self) -> None:
        """Test creating a PhyloNode."""
        node = PhyloNode(haplogroup="H2a2a1")
        assert node.haplogroup == "H2a2a1"
        assert node.parent is None
        assert node.children == []
        assert node.defining_mutations == []

    def test_path_from_root_single(self) -> None:
        """Test path from root for single node."""
        node = PhyloNode(haplogroup="root")
        path = node.get_path_from_root()
        assert len(path) == 1
        assert path[0].haplogroup == "root"

    def test_path_from_root_chain(self) -> None:
        """Test path from root for chain of nodes."""
        root = PhyloNode(haplogroup="root")
        l = PhyloNode(haplogroup="L", parent=root)
        l3 = PhyloNode(haplogroup="L3", parent=l)
        n = PhyloNode(haplogroup="N", parent=l3)

        root.children = [l]
        l.children = [l3]
        l3.children = [n]

        path = n.get_path_from_root()
        assert len(path) == 4
        assert [p.haplogroup for p in path] == ["root", "L", "L3", "N"]

    def test_get_all_expected_mutations(self) -> None:
        """Test getting all expected mutations from root to node."""
        mut1 = Mutation(position=73, ancestral="A", derived="G")
        mut2 = Mutation(position=263, ancestral="A", derived="G")
        mut3 = Mutation(position=8860, ancestral="A", derived="G")

        root = PhyloNode(haplogroup="root", defining_mutations=[])
        h = PhyloNode(haplogroup="H", parent=root, defining_mutations=[mut1])
        h2 = PhyloNode(haplogroup="H2", parent=h, defining_mutations=[mut2])
        h2a = PhyloNode(haplogroup="H2a", parent=h2, defining_mutations=[mut3])

        root.children = [h]
        h.children = [h2]
        h2.children = [h2a]

        mutations = h2a.get_all_expected_mutations()
        assert len(mutations) == 3
        assert mutations[0].position == 73
        assert mutations[1].position == 263
        assert mutations[2].position == 8860

    def test_distance_to_self(self) -> None:
        """Test distance to self is 0."""
        node = PhyloNode(haplogroup="H")
        assert node.distance_to(node) == 0

    def test_distance_to_parent(self) -> None:
        """Test distance to parent is 1."""
        root = PhyloNode(haplogroup="root")
        h = PhyloNode(haplogroup="H", parent=root)
        root.children = [h]

        assert h.distance_to(root) == 1
        assert root.distance_to(h) == 1

    def test_distance_to_sibling(self) -> None:
        """Test distance between siblings is 2."""
        root = PhyloNode(haplogroup="root")
        h = PhyloNode(haplogroup="H", parent=root)
        l = PhyloNode(haplogroup="L", parent=root)
        root.children = [h, l]

        assert h.distance_to(l) == 2


class TestPhylotree:
    """Tests for Phylotree class."""

    def test_create_empty_tree(self) -> None:
        """Test creating an empty Phylotree."""
        tree = Phylotree()
        assert tree.root is None
        assert len(tree.nodes) == 0
        assert len(tree.mutations) == 0
        assert len(tree.hotspots) == 0

    def test_get_haplogroup_exists(self) -> None:
        """Test getting an existing haplogroup."""
        tree = Phylotree()
        node = PhyloNode(haplogroup="H2a2a1")
        tree.nodes["H2a2a1"] = node

        result = tree.get_haplogroup("H2a2a1")
        assert result == node

    def test_get_haplogroup_not_exists(self) -> None:
        """Test getting a non-existent haplogroup."""
        tree = Phylotree()
        result = tree.get_haplogroup("NONEXISTENT")
        assert result is None

    def test_get_defining_positions(self) -> None:
        """Test getting all defining positions."""
        tree = Phylotree()
        tree.mutations[73] = [Mutation(73, "A", "G")]
        tree.mutations[263] = [Mutation(263, "A", "G")]
        tree.mutations[8860] = [Mutation(8860, "A", "G")]

        positions = tree.get_defining_positions()
        assert positions == {73, 263, 8860}

    def test_get_mutations_at_position(self) -> None:
        """Test getting mutations at a position."""
        tree = Phylotree()
        mut1 = Mutation(73, "A", "G")
        mut2 = Mutation(73, "A", "T")  # Different mutation at same position
        tree.mutations[73] = [mut1, mut2]

        mutations = tree.get_mutations_at(73)
        assert len(mutations) == 2
        assert mut1 in mutations
        assert mut2 in mutations

    def test_get_mutations_at_empty_position(self) -> None:
        """Test getting mutations at position with no mutations."""
        tree = Phylotree()
        mutations = tree.get_mutations_at(73)
        assert mutations == []

    def test_is_hotspot_true(self) -> None:
        """Test is_hotspot returns True for hotspot."""
        tree = Phylotree()
        tree.hotspots = {310, 315, 16519}

        assert tree.is_hotspot(310) is True
        assert tree.is_hotspot(315) is True
        assert tree.is_hotspot(16519) is True

    def test_is_hotspot_false(self) -> None:
        """Test is_hotspot returns False for non-hotspot."""
        tree = Phylotree()
        tree.hotspots = {310, 315, 16519}

        assert tree.is_hotspot(73) is False
        assert tree.is_hotspot(263) is False


class TestPhylotreeLoad:
    """Tests for Phylotree.load() functionality."""

    @pytest.fixture
    def evehap_data_dir(self) -> Path:
        """Return eveHap data directory."""
        return Path(__file__).parent.parent.parent / "data"

    def test_load_tree_xml(self, evehap_data_dir: Path) -> None:
        """Test loading phylotree from tree.xml."""
        tree_path = evehap_data_dir / "phylotree" / "tree.xml"

        if not tree_path.exists():
            pytest.skip("tree.xml not available")

        tree = Phylotree.load(str(tree_path))

        # Should have loaded haplogroups
        assert len(tree.nodes) > 0

        # Should have a root
        assert tree.root is not None

        # Should have common haplogroups
        assert tree.get_haplogroup("H") is not None or tree.get_haplogroup("H2a2a1") is not None

    def test_load_with_weights(self, evehap_data_dir: Path) -> None:
        """Test loading phylotree with weights file."""
        tree_path = evehap_data_dir / "phylotree" / "tree.xml"
        weights_path = evehap_data_dir / "phylotree" / "weights.txt"

        if not tree_path.exists() or not weights_path.exists():
            pytest.skip("tree.xml or weights.txt not available")

        tree = Phylotree.load(str(tree_path), str(weights_path))

        # Mutations should have weights assigned
        for pos, mutations in tree.mutations.items():
            for mut in mutations:
                # Most mutations should have weight > 0
                assert mut.weight >= 0

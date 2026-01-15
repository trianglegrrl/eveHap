"""Phylotree parser and data structures for eveHap.

Parses Haplogrep3-format phylotree XML files and provides tree navigation.
"""

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple

if TYPE_CHECKING:
    from evehap.core.profile import AlleleProfile

# rCRS reference for determining ancestral alleles
# This is a minimal mapping for the most common positions
# Full reference is loaded from FASTA when available
RCRS_REFERENCE: Dict[int, str] = {}


def _parse_mutation_string(mut_str: str) -> Tuple[int, str, str, bool]:
    """Parse a mutation string like '73G', '315.1C', '522-523d'.

    Args:
        mut_str: Mutation string in phylotree format

    Returns:
        Tuple of (position, ancestral, derived, is_backmutation)

    Raises:
        ValueError: If mutation string cannot be parsed
    """
    mut_str = mut_str.strip()

    # Handle back mutations (prefixed with !)
    is_backmutation = mut_str.startswith("!")
    if is_backmutation:
        mut_str = mut_str[1:]

    # Handle deletions like '522-523d' or '3107d'
    if mut_str.endswith("d"):
        # Deletion
        if "-" in mut_str:
            # Range deletion like '522-523d'
            parts = mut_str[:-1].split("-")
            position = int(parts[0])
        else:
            # Single deletion like '3107d'
            position = int(mut_str[:-1])
        return (position, "N", "-", is_backmutation)

    # Handle insertions like '315.1C' or '309.1CC'
    if "." in mut_str:
        # Insertion
        parts = mut_str.split(".")
        position = int(parts[0])
        # The derived allele is everything after the dot's number
        match = re.match(r"(\d+)([A-Z]+)", parts[1])
        if match:
            derived = match.group(2)
        else:
            derived = parts[1]
        return (position, "-", derived, is_backmutation)

    # Standard substitution like '73G' or '16519C'
    match = re.match(r"(\d+)([A-CGTUN]+)$", mut_str, re.IGNORECASE)
    if match:
        position = int(match.group(1))
        derived = match.group(2).upper()
        # Get ancestral from rCRS if available, otherwise use 'N'
        ancestral = RCRS_REFERENCE.get(position, "N")
        return (position, ancestral, derived, is_backmutation)

    raise ValueError(f"Cannot parse mutation string: {mut_str}")


@dataclass
class Mutation:
    """A phylogenetically informative mutation.

    Attributes:
        position: 1-based mtDNA position
        ancestral: Ancestral allele
        derived: Derived allele
        weight: Phylogenetic weight (rarer mutations = higher weight)
        is_backmutation: True if this reverts to ancestral state
    """

    position: int
    ancestral: str
    derived: str
    weight: float = 1.0
    is_backmutation: bool = False

    def matches(self, profile: "AlleleProfile") -> Optional[bool]:
        """Check if profile has this mutation.

        Args:
            profile: AlleleProfile to check

        Returns:
            True if mutation present, False if ancestral, None if not covered
        """
        return profile.is_derived(self.position, self.derived)

    def __repr__(self) -> str:
        """Return string representation like '73G'."""
        return f"{self.position}{self.derived}"


@dataclass
class PhyloNode:
    """A node in the phylotree (represents a haplogroup).

    Attributes:
        haplogroup: Haplogroup name (e.g., "H2a2a1")
        parent: Parent node (None for root)
        children: List of child nodes
        defining_mutations: Mutations that define this branch
    """

    haplogroup: str
    parent: Optional["PhyloNode"] = None
    children: List["PhyloNode"] = field(default_factory=list)
    defining_mutations: List[Mutation] = field(default_factory=list)

    def get_path_from_root(self) -> List["PhyloNode"]:
        """Return path from root to this node."""
        path = []
        node: Optional[PhyloNode] = self
        while node is not None:
            path.append(node)
            node = node.parent
        return list(reversed(path))

    def get_all_expected_mutations(self) -> List[Mutation]:
        """Return all mutations expected for this haplogroup (cumulative from root)."""
        mutations = []
        for node in self.get_path_from_root():
            mutations.extend(node.defining_mutations)
        return mutations

    def distance_to(self, other: "PhyloNode") -> int:
        """Phylogenetic distance (number of branches) to another node."""
        # Find common ancestor
        self_path = {n.haplogroup for n in self.get_path_from_root()}
        other_path = other.get_path_from_root()

        # Find where paths diverge
        common = None
        for node in other_path:
            if node.haplogroup in self_path:
                common = node
                break

        if common is None:
            raise ValueError(f"No common ancestor between {self.haplogroup} and {other.haplogroup}")

        # Count edges from self to common + common to other
        self_dist = 0
        current: Optional[PhyloNode] = self
        while current and current.haplogroup != common.haplogroup:
            self_dist += 1
            current = current.parent

        other_dist = 0
        current = other
        while current and current.haplogroup != common.haplogroup:
            other_dist += 1
            current = current.parent

        return self_dist + other_dist


class Phylotree:
    """Complete phylotree structure.

    Attributes:
        root: Root node of the tree
        nodes: Dictionary mapping haplogroup name to node
        mutations: Dictionary mapping position to mutations at that position
        hotspots: Set of known hotspot positions
    """

    def __init__(self) -> None:
        """Initialize an empty Phylotree."""
        self.root: Optional[PhyloNode] = None
        self.nodes: Dict[str, PhyloNode] = {}
        self.mutations: Dict[int, List[Mutation]] = {}
        self.hotspots: Set[int] = set()

    @classmethod
    def load(cls, tree_path: str, weights_path: Optional[str] = None) -> "Phylotree":
        """Load phylotree from file (auto-detects format).

        Args:
            tree_path: Path to tree file (.xml or .json)
            weights_path: Optional path to weights.txt file (XML only)

        Returns:
            Loaded Phylotree

        Raises:
            FileNotFoundError: If tree file doesn't exist
            ValueError: If tree file cannot be parsed
        """
        tree_file = Path(tree_path)
        if not tree_file.exists():
            raise FileNotFoundError(f"Tree file not found: {tree_path}")

        # Auto-detect format
        if tree_file.suffix.lower() == ".json":
            return cls.load_mitoleaf(tree_path)
        else:
            return cls.load_xml(tree_path, weights_path)

    @classmethod
    def load_mitoleaf(cls, tree_path: str) -> "Phylotree":
        """Load phylotree from mitoLeaf JSON format.

        mitoLeaf format has nested nodes with:
        - name: haplogroup name
        - HG: defining mutations as space-separated string
        - children: list of child nodes

        Args:
            tree_path: Path to tree.json file

        Returns:
            Loaded Phylotree
        """
        import json

        tree_file = Path(tree_path)
        if not tree_file.exists():
            raise FileNotFoundError(f"Tree file not found: {tree_path}")

        with open(tree_file) as f:
            data = json.load(f)

        phylotree = cls()

        # Create synthetic root
        phylotree.root = PhyloNode(haplogroup="root")
        phylotree.nodes["root"] = phylotree.root

        # Parse the mitoLeaf tree
        phylotree._parse_mitoleaf_node(data, phylotree.root)

        return phylotree

    def _parse_mitoleaf_node(
        self,
        node_data: dict,
        parent: PhyloNode,
    ) -> PhyloNode:
        """Parse a mitoLeaf JSON node recursively.

        Args:
            node_data: Dictionary with name, HG, children
            parent: Parent PhyloNode

        Returns:
            Parsed PhyloNode
        """
        haplogroup_name = node_data.get("name", "unknown")

        # Create node
        node = PhyloNode(haplogroup=haplogroup_name, parent=parent)
        parent.children.append(node)
        self.nodes[haplogroup_name] = node

        # Parse defining mutations from HG string
        hg_string = node_data.get("HG", "")
        if hg_string:
            for mut_str in hg_string.split():
                try:
                    pos, anc, der, is_back = _parse_mutation_string(mut_str)
                    mutation = Mutation(
                        position=pos,
                        ancestral=anc,
                        derived=der,
                        weight=1.0,
                        is_backmutation=is_back,
                    )
                    node.defining_mutations.append(mutation)

                    # Add to position-based index
                    if pos not in self.mutations:
                        self.mutations[pos] = []
                    self.mutations[pos].append(mutation)
                except ValueError:
                    # Skip unparseable mutations
                    pass

        # Recursively parse children
        for child_data in node_data.get("children", []):
            self._parse_mitoleaf_node(child_data, node)

        return node

    @classmethod
    def load_xml(cls, tree_path: str, weights_path: Optional[str] = None) -> "Phylotree":
        """Load phylotree from Haplogrep3 XML format.

        Args:
            tree_path: Path to tree.xml file
            weights_path: Optional path to weights.txt file

        Returns:
            Loaded Phylotree

        Raises:
            FileNotFoundError: If tree file doesn't exist
            ValueError: If tree file cannot be parsed
        """
        tree_file = Path(tree_path)
        if not tree_file.exists():
            raise FileNotFoundError(f"Tree file not found: {tree_path}")

        phylotree = cls()

        # Load weights if provided
        weights: Dict[str, float] = {}
        if weights_path:
            weights = cls._load_weights(weights_path)

        # Parse XML tree
        try:
            xml_tree = ET.parse(tree_path)
            root_element = xml_tree.getroot()
        except ET.ParseError as e:
            raise ValueError(f"Failed to parse tree XML: {e}") from e

        # The root element should be <phylotree>
        if root_element.tag != "phylotree":
            raise ValueError(f"Expected <phylotree> root element, got <{root_element.tag}>")

        # Create a synthetic root node
        phylotree.root = PhyloNode(haplogroup="root")
        phylotree.nodes["root"] = phylotree.root

        # Parse all haplogroups recursively
        # The XML structure has nested <haplogroup> elements
        # We need to find the top-level haplogroups first
        for child_element in root_element:
            if child_element.tag == "haplogroup":
                phylotree._parse_haplogroup_element(child_element, phylotree.root, weights)

        return phylotree

    def _parse_haplogroup_element(
        self,
        element: ET.Element,
        parent: PhyloNode,
        weights: Dict[str, float],
    ) -> PhyloNode:
        """Parse a <haplogroup> element and its children.

        Args:
            element: XML element representing a haplogroup
            parent: Parent PhyloNode
            weights: Dictionary of mutation weights

        Returns:
            Parsed PhyloNode
        """
        haplogroup_name = element.get("name", "unknown")

        # Create node
        node = PhyloNode(haplogroup=haplogroup_name, parent=parent)
        parent.children.append(node)
        self.nodes[haplogroup_name] = node

        # Parse defining mutations from <details><poly>...</poly></details>
        details = element.find("details")
        if details is not None:
            for poly in details.findall("poly"):
                if poly.text:
                    try:
                        pos, anc, der, is_back = _parse_mutation_string(poly.text)
                        weight = weights.get(f"{pos}{der}", 1.0)
                        mutation = Mutation(
                            position=pos,
                            ancestral=anc,
                            derived=der,
                            weight=weight,
                            is_backmutation=is_back,
                        )
                        node.defining_mutations.append(mutation)

                        # Add to position-based index
                        if pos not in self.mutations:
                            self.mutations[pos] = []
                        self.mutations[pos].append(mutation)
                    except ValueError:
                        # Skip unparseable mutations
                        pass

        # Recursively parse child haplogroups
        for child_element in element:
            if child_element.tag == "haplogroup":
                self._parse_haplogroup_element(child_element, node, weights)

        return node

    @staticmethod
    def _load_weights(weights_path: str) -> Dict[str, float]:
        """Load mutation weights from weights.txt file.

        The weights file format is:
        mutation<tab>weight1<tab>weight2<tab>weight3<tab>weight4<tab>count
        e.g.: 152C	1.0	1.0	1.0	1.0	275.0

        Args:
            weights_path: Path to weights.txt file

        Returns:
            Dictionary mapping mutation string to weight
        """
        weights: Dict[str, float] = {}
        weights_file = Path(weights_path)

        if not weights_file.exists():
            return weights

        with open(weights_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    mutation = parts[0]
                    try:
                        # Use the first weight column
                        weight = float(parts[1])
                        weights[mutation] = weight
                    except ValueError:
                        pass

        return weights

    def get_haplogroup(self, name: str) -> Optional[PhyloNode]:
        """Get node by haplogroup name.

        Args:
            name: Haplogroup name (e.g., "H2a2a1")

        Returns:
            PhyloNode or None if not found
        """
        return self.nodes.get(name)

    def get_root(self) -> Optional[PhyloNode]:
        """Get the root node of the phylotree.

        Returns:
            Root PhyloNode or None if tree is empty
        """
        return self.root

    def get_defining_positions(self) -> Set[int]:
        """Get all positions that define haplogroups.

        Returns:
            Set of all mtDNA positions that have defining mutations
        """
        return set(self.mutations.keys())

    def get_mutations_at(self, position: int) -> List[Mutation]:
        """Get all mutations at a position.

        Args:
            position: 1-based mtDNA position

        Returns:
            List of Mutations at this position
        """
        return self.mutations.get(position, [])

    def is_hotspot(self, position: int) -> bool:
        """Check if position is a known hotspot.

        Args:
            position: 1-based mtDNA position

        Returns:
            True if this is a hypermutable hotspot position
        """
        return position in self.hotspots

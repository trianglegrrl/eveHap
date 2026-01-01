"""Classification algorithms for eveHap.

Implements two classification strategies:
1. KulczynskiScorer: Best for high-coverage data (>80% mtDNA covered)
2. TreeTraverser: Best for low-coverage/ancient DNA data

The unified Classifier automatically selects the appropriate strategy.

The Kulczynski scorer uses mutation weights (phylogenetic weights) as in Haplogrep3:
- Each mutation has a weight based on rarity (common mutations = low weight, rare = high)
- The Kulczynski formula uses WEIGHTED sums, not raw counts
- This naturally prevents basal haplogroups from dominating
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

from evehap.core.phylotree import Phylotree, PhyloNode, Mutation
from evehap.core.profile import AlleleProfile
from evehap.output.result import HaplogroupResult

if TYPE_CHECKING:
    from evehap.core.damage import DamageFilter


def load_mutation_weights(weights_path: Optional[str] = None) -> Dict[str, float]:
    """Load mutation weights from Haplogrep3-format weights file.

    Format: MUTATION<tab>WEIGHT1<tab>WEIGHT2<tab>...
    We use the first weight column.

    Args:
        weights_path: Path to weights file. If None, uses bundled default.

    Returns:
        Dict mapping mutation strings (e.g., "73G", "263G") to weights.
    """
    weights: Dict[str, float] = {}

    if weights_path is None:
        # Use bundled weights file
        # Try package data directory first, then project data directory
        module_dir = Path(__file__).parent.parent
        project_dir = Path(__file__).parent.parent.parent

        # Try both RSRS and rCRS weights files
        for weights_name in ["weights-rsrs.txt", "weights-rcrs.txt", "weights.txt"]:
            for base_dir in [module_dir, project_dir]:
                candidate = base_dir / "data" / "phylotree" / weights_name
                if candidate.exists():
                    weights_path = str(candidate)
                    break
            if weights_path:
                break

    path = Path(weights_path)
    if not path.exists():
        return weights

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) >= 2:
                mutation = parts[0].strip()
                # Remove backmutation marker (!) for lookup
                mutation_clean = mutation.rstrip("!")
                try:
                    weight = float(parts[1])
                    weights[mutation_clean] = weight
                    # Also store with backmutation marker
                    if mutation.endswith("!"):
                        weights[mutation] = weight
                except ValueError:
                    continue

    return weights


@dataclass
class TraversalResult:
    """Result of tree traversal classification.

    Attributes:
        path: Path from root to final node
        haplogroup: Final haplogroup assignment
        confidence: Confidence score (0.0-1.0)
        derived_count: Total derived alleles found
        ancestral_count: Total unexpected ancestral alleles
        uncovered_count: Positions not covered in profile
        stopped_early: True if traversal stopped before reaching a leaf
    """

    path: List[PhyloNode] = field(default_factory=list)
    haplogroup: str = ""
    confidence: float = 0.0
    derived_count: int = 0
    ancestral_count: int = 0
    uncovered_count: int = 0
    stopped_early: bool = False


class KulczynskiScorer:
    """Calculate Kulczynski similarity scores using weighted mutations.

    The Kulczynski measure (as implemented in Haplogrep3) is:
        score = 0.5 * (HaplogroupWeight + SampleWeight)

    Where:
        HaplogroupWeight = WeightFoundPolys / WeightExpectedPolys
        SampleWeight = WeightFoundPolys / WeightSamplePolys

    Each mutation has a phylogenetic weight:
        - Common/hotspot mutations have LOW weights (~1.0-2.0)
        - Rare/informative mutations have HIGH weights (~8.0-10.0)

    This naturally prevents basal haplogroups from dominating, as they
    typically have common mutations with low weights.
    """

    # Default weight for mutations not in the weights file
    DEFAULT_WEIGHT = 5.0  # Middle of the 1-10 range

    def __init__(
        self,
        phylotree: Phylotree,
        reference: Optional[str] = None,
        weights_path: Optional[str] = None,
    ) -> None:
        """Initialize scorer with phylotree.

        Args:
            phylotree: Phylotree reference
            reference: Optional reference sequence (rCRS) for filling missing positions.
                      If provided, missing positions will be assumed to match reference.
            weights_path: Optional path to weights file. If None, uses bundled default.
        """
        self.phylotree = phylotree
        self.reference = reference
        self._cache: Dict[str, Set[int]] = {}
        self.weights = load_mutation_weights(weights_path)

    def _get_mutation_weight(self, position: int, derived: str) -> float:
        """Get weight for a mutation.

        Looks up the mutation in format like "73G", "263G", etc.

        Args:
            position: Mutation position (1-based)
            derived: Derived allele

        Returns:
            Weight for the mutation, or DEFAULT_WEIGHT if not found.
        """
        # Standard format: position + derived allele
        mutation_str = f"{position}{derived}"
        if mutation_str in self.weights:
            return self.weights[mutation_str]

        # Some weights file entries have format like "8281-8289d" for deletions
        # For deletions (derived="-"), check for deletion range format
        if derived == "-":
            # Try single position deletion
            del_str = f"{position}d"
            if del_str in self.weights:
                return self.weights[del_str]

        return self.DEFAULT_WEIGHT

    def _get_expected_positions(self, node: PhyloNode) -> Set[Tuple[int, str]]:
        """Get all expected (position, derived) pairs for a haplogroup.

        Args:
            node: Haplogroup node

        Returns:
            Set of (position, derived_allele) tuples
        """
        cache_key = node.haplogroup
        if cache_key not in self._cache:
            mutations = node.get_all_expected_mutations()
            expected = set()
            for m in mutations:
                if not m.is_backmutation:
                    expected.add((m.position, m.derived))
                else:
                    # Backmutation removes the derived allele
                    expected.discard((m.position, m.ancestral))
            self._cache[cache_key] = expected
        return self._cache[cache_key]

    def precompute_sample_info(
        self, profile: AlleleProfile
    ) -> Tuple[Set[Tuple[int, str]], float]:
        """Pre-compute sample variants and weight once for all haplogroups.

        This is a MAJOR optimization - avoids iterating 15k+ observations
        for each of 6k+ haplogroups.

        Args:
            profile: Sample profile

        Returns:
            Tuple of (sample_variants set, weight_sample_all)
        """
        weight_sample_all = 0.0
        sample_variants: Set[Tuple[int, str]] = set()

        for pos, obs in profile.observations.items():
            if obs.is_variant:
                sample_variants.add((pos, obs.major_allele))
                # Only add to sample_weight if phylogenetically informative
                if pos in self.phylotree.mutations:
                    tree_muts = self.phylotree.mutations[pos]
                    if any(m.derived == obs.major_allele for m in tree_muts):
                        mutation_weight = self._get_mutation_weight(pos, obs.major_allele)
                        weight_sample_all += mutation_weight

        return sample_variants, weight_sample_all

    def score(
        self,
        profile: AlleleProfile,
        haplogroup: PhyloNode,
        sample_info: Optional[Tuple[Set[Tuple[int, str]], float]] = None,
    ) -> float:
        """Calculate Kulczynski score for a haplogroup using weighted mutations.

        Implements Haplogrep3's Kulczynski formula:
            score = 0.5 * (HaplogroupWeight + SampleWeight)

        Where:
            HaplogroupWeight = WeightFoundPolys / WeightExpectedPolys
            SampleWeight = WeightFoundPolys / WeightSamplePolys

        Args:
            profile: Sample profile
            haplogroup: Haplogroup node to score
            sample_info: Pre-computed (sample_variants, weight_sample_all) for efficiency

        Returns:
            Score between 0.0 (no match) and 1.0 (perfect match)
        """
        expected = self._get_expected_positions(haplogroup)

        if not expected:
            return 0.0

        # Use pre-computed sample info if provided (major optimization)
        if sample_info is not None:
            sample_variants, weight_sample_all = sample_info
        else:
            # Fallback: compute here (slower for batch scoring)
            sample_variants, weight_sample_all = self.precompute_sample_info(profile)

        # For VCF data (assumes_full_coverage), positions not in observations are reference
        # We need to add weight for expected mutations that match reference
        uses_reference = (
            hasattr(profile, 'assumes_full_coverage')
            and profile.assumes_full_coverage
            and self.reference is not None
        )

        # Calculate WEIGHTED sums for expected mutations
        weight_expected = 0.0  # Sum of weights for expected mutations in covered range
        weight_found = 0.0     # Sum of weights for found mutations

        for pos, derived in expected:
            mutation_weight = self._get_mutation_weight(pos, derived)

            obs = profile.get_observation(pos)
            if obs is not None:
                # Position is covered in profile
                weight_expected += mutation_weight
                if obs.major_allele == derived:
                    weight_found += mutation_weight
            elif uses_reference and pos <= len(self.reference):
                # VCF mode: position not in profile means reference allele
                weight_expected += mutation_weight
                ref_allele = self.reference[pos - 1].upper()  # 0-indexed
                if ref_allele == derived:
                    weight_found += mutation_weight
                    # For VCF mode, also add to sample_all if expected=ref
                    # (reference positions don't add to sample weight)
            # If position is not covered and not VCF mode, skip it

        if weight_expected == 0:
            return 0.0

        # Haplogrep3 Kulczynski formula
        haplogroup_weight = weight_found / weight_expected

        # Sample weight: what fraction of sample's variants are explained?
        # For this, we count weight of variants that match expected mutations
        weight_found_sample = 0.0
        for pos, derived in expected:
            if (pos, derived) in sample_variants:
                mutation_weight = self._get_mutation_weight(pos, derived)
                weight_found_sample += mutation_weight

        if weight_sample_all > 0:
            sample_weight = weight_found_sample / weight_sample_all
        else:
            # No variants = reference sequence = sample_weight based on haplogroup coverage
            sample_weight = haplogroup_weight

        # Kulczynski measure with adaptive weighting
        # For microarray data (sparse observations, many uninformative variants),
        # weight sample_weight more heavily since we can't achieve 100% haplogroup match
        # For VCF/BAM data (full coverage), use standard 50/50 weighting
        if len(profile.observations) < 1000:
            # Sparse data (microarray) - weight sample_weight 2:1
            score = 0.33 * haplogroup_weight + 0.67 * sample_weight
        else:
            # Full coverage data (VCF/BAM)
            score = 0.5 * (haplogroup_weight + sample_weight)

        # Add a tiny tie-breaking bonus for phylogenetic depth
        depth = len(haplogroup.get_all_expected_mutations())
        depth_bonus = depth * 0.00001
        score += depth_bonus

        # Clamp to [0, 1]
        return max(0.0, min(1.0, score))

    def score_all(
        self, profile: AlleleProfile, top_n: int = 10, use_pruning: bool = True
    ) -> List[Tuple[str, float]]:
        """Score haplogroups and return top N.

        Uses weighted Kulczynski scoring like Haplogrep3. Supports two modes:

        1. Exhaustive (use_pruning=False): Scores all haplogroups
        2. Pruning (use_pruning=True): Uses tree-guided pruning for faster scoring

        Args:
            profile: Sample profile
            top_n: Number of top results to return
            use_pruning: If True, uses tree-pruning optimization

        Returns:
            List of (haplogroup, score) tuples sorted by score descending
        """
        if use_pruning:
            return self.score_with_pruning(profile, top_n=top_n)

        # Pre-compute sample info ONCE (major optimization)
        sample_info = self.precompute_sample_info(profile)

        scores: List[Tuple[str, float, int]] = []  # (hg, score, depth)

        # Score all non-root nodes
        for hg_name, node in self.phylotree.nodes.items():
            if hg_name == "root":
                continue

            score = self.score(profile, node, sample_info=sample_info)
            if score > 0:
                # Track phylogenetic depth (number of expected mutations)
                depth = len(node.get_all_expected_mutations())
                scores.append((hg_name, score, depth))

        if not scores:
            return []

        # Sort by score descending, then by depth (more mutations = deeper)
        # This tie-breaks near-equal scores in favor of more specific haplogroups
        scores.sort(key=lambda x: (x[1], x[2]), reverse=True)

        return [(hg, score) for hg, score, _ in scores[:top_n]]

    def score_with_pruning(
        self,
        profile: AlleleProfile,
        min_score: float = 0.15,
        top_n: int = 10,
    ) -> List[Tuple[str, float]]:
        """Score using tree-guided pruning instead of exhaustive search.

        Uses a low min_score (0.15) because for sparse data like microarrays,
        intermediate nodes on the path to the correct haplogroup can score low
        even when the terminal haplogroup scores high.

        Strategy:
        1. Start at root, score all immediate children
        2. Keep exploring children of nodes that score above min_score
        3. Use a priority queue to explore most promising branches first
        4. Early termination when we've explored enough high-scoring nodes

        This typically scores 50-200 haplogroups instead of 5000+, providing
        10-100x speedup for classification.

        Note: For rCRS trees where the root has few children (single reference
        haplogroup), pruning may not be effective. In such cases, we fall back
        to exhaustive search.

        Args:
            profile: Sample profile
            min_score: Minimum score to continue exploring children
            top_n: Number of top results to return

        Returns:
            List of (haplogroup, score) tuples sorted by score descending
        """
        import heapq

        # Find root node
        root_node = self.phylotree.get_root()
        if root_node is None:
            # Fallback to exhaustive if no root found
            return self.score_all(profile, top_n=top_n, use_pruning=False)

        # Walk down single-child chains to find the "real" root for pruning
        # (e.g., root -> mtMRCA -> L0/L1'2'3'4'5'6 in RSRS trees)
        start_node = root_node
        while len(start_node.children) == 1:
            start_node = start_node.children[0]

        # If we've descended to a leaf or very sparse tree, fall back to exhaustive
        if len(start_node.children) < 2:
            return self.score_all(profile, top_n=top_n, use_pruning=False)

        # Pre-compute sample info ONCE (major optimization)
        sample_info = self.precompute_sample_info(profile)

        scores: List[Tuple[str, float]] = []
        scored_count = 0

        # Priority queue: (negative_score, haplogroup, node)
        # Use negative score because heapq is a min-heap
        candidates: List[Tuple[float, str, PhyloNode]] = []
        visited: Set[str] = set()

        # Mark all ancestors as visited (we don't need to score them)
        node = start_node
        while node is not None:
            visited.add(node.haplogroup)
            node = node.parent

        # Start with start_node's children
        for child in start_node.children:
            heapq.heappush(candidates, (0.0, child.haplogroup, child))

        # Best scores found so far (for early termination)
        best_scores: List[Tuple[float, str]] = []

        while candidates:
            # Get most promising node (highest score = most negative priority)
            neg_parent_score, hg_name, node = heapq.heappop(candidates)

            if hg_name in visited:
                continue
            visited.add(hg_name)

            # Score this node using pre-computed sample info
            score = self.score(profile, node, sample_info=sample_info)
            scored_count += 1

            if score > 0:
                scores.append((hg_name, score))
                heapq.heappush(best_scores, (score, hg_name))
                if len(best_scores) > top_n * 2:
                    heapq.heappop(best_scores)

            # Explore children if score is promising
            if score >= min_score:
                for child in node.children:
                    if child.haplogroup not in visited:
                        # Priority = negative score (higher score = higher priority)
                        heapq.heappush(candidates, (-score, child.haplogroup, child))

            # Early termination: if we have enough very high-quality results
            # Only terminate early if we have truly excellent matches
            if len(best_scores) >= top_n:
                min_best = best_scores[0][0]  # Smallest score in top N
                # Only terminate if all top_n scores are > 0.8 and we've scored many nodes
                if min_best > 0.8 and scored_count > 500:
                    # We have excellent results, safe to stop
                    break

        if not scores:
            return []

        # Sort by score descending
        scores.sort(key=lambda x: x[1], reverse=True)

        return scores[:top_n]


class TreeTraverser:
    """Tree traversal classification for sparse data.

    For low-coverage ancient DNA, traverses the tree evaluating
    support at each branch. Uses a greedy algorithm that:
    1. Starts at root
    2. For each child, counts derived vs ancestral alleles
    3. Takes the child with best support ratio
    4. Stops when no children have sufficient support
    """

    def __init__(
        self,
        phylotree: Phylotree,
        max_conflicts: int = 3,
        min_support_ratio: float = 0.7,
    ) -> None:
        """Initialize traverser.

        Args:
            phylotree: Phylotree to traverse
            max_conflicts: Max ancestral alleles before stopping
            min_support_ratio: Min derived/(derived+ancestral) to continue
        """
        self.phylotree = phylotree
        self.max_conflicts = max_conflicts
        self.min_support_ratio = min_support_ratio

    def _evaluate_branch(
        self, profile: AlleleProfile, node: PhyloNode
    ) -> Tuple[int, int, int]:
        """Evaluate how well profile supports a branch.

        Args:
            profile: Sample profile
            node: Node to evaluate

        Returns:
            Tuple of (derived_count, ancestral_count, uncovered_count)
        """
        derived = 0
        ancestral = 0
        uncovered = 0

        for mutation in node.defining_mutations:
            result = mutation.matches(profile)

            if result is None:
                uncovered += 1
            elif result:
                derived += 1
            else:
                ancestral += 1

        return derived, ancestral, uncovered

    def traverse(self, profile: AlleleProfile) -> TraversalResult:
        """Traverse tree to find best haplogroup placement.

        Starting from root, evaluates each branch's support.
        Stops when conflicts exceed threshold or no children have support.

        Args:
            profile: Sample profile

        Returns:
            TraversalResult with path, final haplogroup, and confidence
        """
        result = TraversalResult()

        if self.phylotree.root is None:
            result.haplogroup = "unknown"
            return result

        current = self.phylotree.root
        result.path.append(current)

        total_derived = 0
        total_ancestral = 0
        total_uncovered = 0

        while current.children:
            best_child = None
            best_score = -1.0
            best_derived = 0
            best_ancestral = 0
            best_uncovered = 0

            for child in current.children:
                derived, ancestral, uncovered = self._evaluate_branch(profile, child)

                # Skip if too many conflicts
                if ancestral > self.max_conflicts:
                    continue

                # Calculate support ratio
                total_checked = derived + ancestral
                if total_checked > 0:
                    support_ratio = derived / total_checked
                elif len(child.defining_mutations) == 0:
                    # Node has no defining mutations - allow traversal
                    # (e.g., mtMRCA which is just a structural node)
                    support_ratio = 1.0
                else:
                    # No covered positions - allow traversal if uncovered
                    support_ratio = 0.5 if uncovered > 0 else 0.0

                # Only consider branches with sufficient support
                if support_ratio >= self.min_support_ratio or (
                    derived > 0 and ancestral == 0
                ) or len(child.defining_mutations) == 0:
                    if support_ratio > best_score or (
                        support_ratio == best_score and derived > best_derived
                    ):
                        best_child = child
                        best_score = support_ratio
                        best_derived = derived
                        best_ancestral = ancestral
                        best_uncovered = uncovered

            if best_child is None:
                # No valid children - stop here
                result.stopped_early = len(current.children) > 0
                break

            # Move to best child
            current = best_child
            result.path.append(current)
            total_derived += best_derived
            total_ancestral += best_ancestral
            total_uncovered += best_uncovered

        result.haplogroup = current.haplogroup
        result.derived_count = total_derived
        result.ancestral_count = total_ancestral
        result.uncovered_count = total_uncovered

        # Calculate confidence
        total_evaluated = total_derived + total_ancestral
        if total_evaluated > 0:
            result.confidence = total_derived / total_evaluated
        elif total_uncovered > 0:
            result.confidence = 0.5  # Low confidence due to lack of data
        else:
            result.confidence = 0.0

        return result


class Classifier:
    """Unified haplogroup classifier combining multiple strategies.

    Automatically selects between Kulczynski scoring (for high-coverage data)
    and tree traversal (for low-coverage/ancient DNA).
    """

    # Coverage threshold for method selection
    # Set to 0 to always use Kulczynski because:
    # 1. Kulczynski scoring works well with sparse data (even 2% coverage)
    # 2. Traversal method has issues with rCRS trees (returns wrong haplogroup)
    # 3. Haplogrep3 uses Kulczynski for all inputs
    HIGH_COVERAGE_THRESHOLD = 0.0

    def __init__(
        self,
        phylotree: Phylotree,
        method: str = "auto",
        damage_filter: Optional["DamageFilter"] = None,
        reference: Optional[str] = None,
        weights_path: Optional[str] = None,
    ) -> None:
        """Initialize classifier.

        Args:
            phylotree: Phylotree reference
            method: 'kulczynski', 'traversal', or 'auto'
            damage_filter: Optional damage filter for ancient DNA
            reference: Optional reference sequence (rCRS). If provided, missing
                      positions in VCF profiles will be assumed to match reference.
            weights_path: Optional path to mutation weights file. If None, uses bundled default.
        """
        self.phylotree = phylotree
        self.method = method
        self.damage_filter = damage_filter
        self.reference = reference

        self._scorer = KulczynskiScorer(phylotree, reference=reference, weights_path=weights_path)
        self._traverser = TreeTraverser(phylotree)

    def _select_method(self, profile: AlleleProfile) -> str:
        """Select classification method based on profile characteristics.

        Args:
            profile: Sample profile

        Returns:
            Method name ('kulczynski' or 'traversal')
        """
        if self.method != "auto":
            return self.method

        # Use coverage to select method
        coverage = profile.coverage_fraction()
        if coverage >= self.HIGH_COVERAGE_THRESHOLD:
            return "kulczynski"
        else:
            return "traversal"

    def classify(self, profile: AlleleProfile) -> HaplogroupResult:
        """Classify sample and return haplogroup result.

        If method='auto', chooses strategy based on coverage:
        - High coverage (>=80%): Kulczynski scoring
        - Low coverage (<80%): Tree traversal

        Args:
            profile: Sample profile

        Returns:
            HaplogroupResult with classification details
        """
        # Apply damage filter if available
        if self.damage_filter:
            profile = self.damage_filter.filter_profile(profile)

        # Select method
        method = self._select_method(profile)

        # Run classification
        if method == "kulczynski":
            return self._classify_kulczynski(profile)
        else:
            return self._classify_traversal(profile)

    def _classify_kulczynski(self, profile: AlleleProfile) -> HaplogroupResult:
        """Classify using Kulczynski scoring.

        Args:
            profile: Sample profile

        Returns:
            HaplogroupResult
        """
        result = HaplogroupResult(
            sample_id=profile.sample_id,
            method="kulczynski",
            coverage_fraction=profile.coverage_fraction(),
            mean_depth=profile.mean_depth(),
        )

        # Get top scores - use exhaustive search like Haplogrep3
        # (pruning doesn't work well for rCRS trees where the topology is inverted)
        scores = self._scorer.score_all(profile, top_n=10, use_pruning=False)

        if not scores:
            result.haplogroup = "unknown"
            result.confidence = 0.0
            result.quality = "low"
            return result

        # Best match
        best_hg, best_score = scores[0]
        result.haplogroup = best_hg
        result.confidence = best_score

        # Alternative haplogroups
        result.alternative_haplogroups = scores[1:5]

        # Get expected mutations for the best haplogroup
        node = self.phylotree.get_haplogroup(best_hg)
        if node:
            expected_mutations = node.get_all_expected_mutations()
            result.expected_mutations = expected_mutations

            # Categorize mutations
            for mutation in expected_mutations:
                match_result = mutation.matches(profile)
                if match_result is True:
                    result.found_mutations.append(mutation)
                elif match_result is False:
                    result.missing_mutations.append(mutation)

        # Determine quality
        if best_score >= 0.9:
            result.quality = "high"
        elif best_score >= 0.7:
            result.quality = "medium"
        else:
            result.quality = "low"

        # Add warnings
        self._add_warnings(result, profile)

        return result

    def _classify_traversal(self, profile: AlleleProfile) -> HaplogroupResult:
        """Classify using tree traversal.

        Args:
            profile: Sample profile

        Returns:
            HaplogroupResult
        """
        result = HaplogroupResult(
            sample_id=profile.sample_id,
            method="traversal",
            coverage_fraction=profile.coverage_fraction(),
            mean_depth=profile.mean_depth(),
        )

        # Traverse tree
        traversal = self._traverser.traverse(profile)

        result.haplogroup = traversal.haplogroup
        result.confidence = traversal.confidence

        # Get expected mutations
        node = self.phylotree.get_haplogroup(traversal.haplogroup)
        if node:
            expected_mutations = node.get_all_expected_mutations()
            result.expected_mutations = expected_mutations

            # Categorize mutations based on traversal
            for mutation in expected_mutations:
                match_result = mutation.matches(profile)
                if match_result is True:
                    result.found_mutations.append(mutation)
                elif match_result is False:
                    result.missing_mutations.append(mutation)

        # Determine quality
        if traversal.confidence >= 0.9 and not traversal.stopped_early:
            result.quality = "high"
        elif traversal.confidence >= 0.6:
            result.quality = "medium"
        else:
            result.quality = "low"

        # Add warnings
        self._add_warnings(result, profile)

        if traversal.stopped_early:
            result.warnings.append("MISSING_KEY_MUTATIONS")

        return result

    def _add_warnings(self, result: HaplogroupResult, profile: AlleleProfile) -> None:
        """Add QC warnings to result.

        Args:
            result: HaplogroupResult to add warnings to
            profile: Source profile
        """
        # Check coverage
        if profile.coverage_fraction() < 0.5:
            result.warnings.append("LOW_COVERAGE")

        # Check depth
        mean_depth = profile.mean_depth()
        if mean_depth is not None and mean_depth < 10:
            result.warnings.append("LOW_DEPTH")

        # Check for heteroplasmy
        heteroplasmic_count = sum(
            1 for obs in profile.observations.values()
            if obs.is_heteroplasmic
        )
        if heteroplasmic_count > 0:
            result.warnings.append("HETEROPLASMY")


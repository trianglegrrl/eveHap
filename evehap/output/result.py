"""Result data structures for eveHap.

Stub module - to be implemented.
"""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from evehap.core.damage import DamageStats
    from evehap.core.phylotree import Mutation


# QC warning codes and messages
QC_WARNINGS = {
    "LOW_COVERAGE": "Coverage below 50% - reduced confidence",
    "LOW_DEPTH": "Mean depth below 10x - may have errors",
    "HETEROPLASMY": "Heteroplasmic positions detected",
    "MISSING_KEY_MUTATIONS": "Key defining mutations not covered",
    "CONFLICTING_MUTATIONS": "Multiple conflicting mutations found",
    "DAMAGE_DETECTED": "Ancient DNA damage patterns detected",
    "HOTSPOT_VARIANT": "Variant at known hotspot position",
}


@dataclass
class HaplogroupResult:
    """Complete result of haplogroup classification.

    Attributes:
        sample_id: Sample identifier
        haplogroup: Best haplogroup assignment
        confidence: Confidence score (0.0-1.0)
        quality: Quality level ('high', 'medium', 'low')
        expected_mutations: Mutations expected for haplogroup
        found_mutations: Mutations found in sample
        missing_mutations: Expected but not found
        extra_mutations: Found but not expected (private mutations)
        method: Classification method ('kulczynski' or 'traversal')
        alternative_haplogroups: Other candidates with scores
        coverage_fraction: Fraction of mtDNA covered
        mean_depth: Mean read depth
        warnings: QC warnings
        damage_filtered: True if damage filter applied
        damage_stats: Damage statistics if available
    """

    sample_id: str = ""
    haplogroup: str = ""
    confidence: float = 0.0
    quality: str = "low"

    expected_mutations: List["Mutation"] = field(default_factory=list)
    found_mutations: List["Mutation"] = field(default_factory=list)
    missing_mutations: List["Mutation"] = field(default_factory=list)
    extra_mutations: List["Mutation"] = field(default_factory=list)

    method: str = ""
    alternative_haplogroups: List[Tuple[str, float]] = field(default_factory=list)

    coverage_fraction: float = 0.0
    mean_depth: Optional[float] = None
    warnings: List[str] = field(default_factory=list)

    damage_filtered: bool = False
    damage_stats: Optional["DamageStats"] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization.

        Returns:
            Dictionary representation of the result
        """
        return {
            "sample_id": self.sample_id,
            "haplogroup": self.haplogroup,
            "confidence": self.confidence,
            "quality": self.quality,
            "method": self.method,
            "coverage_fraction": self.coverage_fraction,
            "mean_depth": self.mean_depth,
            "expected_mutations": [str(m) for m in self.expected_mutations],
            "found_mutations": [str(m) for m in self.found_mutations],
            "missing_mutations": [str(m) for m in self.missing_mutations],
            "extra_mutations": [str(m) for m in self.extra_mutations],
            "alternative_haplogroups": self.alternative_haplogroups,
            "warnings": self.warnings,
            "damage_filtered": self.damage_filtered,
        }

    def to_tsv_row(self) -> str:
        """Convert to TSV row string.

        Returns:
            Tab-separated values string
        """
        return "\t".join(
            [
                self.sample_id,
                self.haplogroup,
                f"{self.confidence:.4f}",
                self.quality,
                f"{self.coverage_fraction:.4f}",
                str(self.mean_depth) if self.mean_depth else "",
                self.method,
                ";".join(self.warnings) if self.warnings else "",
            ]
        )

    @staticmethod
    def tsv_header() -> str:
        """Get TSV header row.

        Returns:
            Tab-separated header string
        """
        return "\t".join(
            [
                "sample_id",
                "haplogroup",
                "confidence",
                "quality",
                "coverage",
                "depth",
                "method",
                "warnings",
            ]
        )

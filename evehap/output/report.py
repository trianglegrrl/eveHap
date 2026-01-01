"""Report generation for eveHap.

Stub module - to be implemented.
"""

from dataclasses import dataclass
from typing import List, TYPE_CHECKING

if TYPE_CHECKING:
    from evehap.output.result import HaplogroupResult


@dataclass
class QCReport:
    """Quality control report for a batch of samples.

    Attributes:
        results: List of HaplogroupResults
        total_samples: Total samples processed
        successful: Number of successful classifications
        warnings_count: Number of samples with warnings
    """

    results: List["HaplogroupResult"]

    @property
    def total_samples(self) -> int:
        """Total number of samples."""
        return len(self.results)

    @property
    def successful(self) -> int:
        """Number of samples with confidence > 0."""
        return sum(1 for r in self.results if r.confidence > 0)

    @property
    def warnings_count(self) -> int:
        """Number of samples with warnings."""
        return sum(1 for r in self.results if r.warnings)

    def to_text(self) -> str:
        """Generate human-readable text report.

        Returns:
            Formatted text report
        """
        lines = [
            "eveHap Classification Report",
            "=" * 40,
            f"Total samples: {self.total_samples}",
            f"Successful: {self.successful}",
            f"With warnings: {self.warnings_count}",
            "",
            "Results:",
            "-" * 40,
        ]

        for result in self.results:
            lines.append(
                f"{result.sample_id}: {result.haplogroup} "
                f"(confidence: {result.confidence:.2f}, quality: {result.quality})"
            )
            if result.warnings:
                for warning in result.warnings:
                    lines.append(f"  ! {warning}")

        return "\n".join(lines)


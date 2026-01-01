"""FASTQ adapter for eveHap.

Stub module - to be implemented.
"""

from typing import TYPE_CHECKING

from evehap.adapters.base import InputAdapter

if TYPE_CHECKING:
    from evehap.core.profile import AlleleProfile


class FASTQAdapter(InputAdapter):
    """Adapter for raw FASTQ read files."""

    def __init__(self, reference_path: str, aligner: str = "minimap2") -> None:
        """Initialize FASTQ adapter.

        Args:
            reference_path: Path to reference FASTA
            aligner: Alignment tool to use
        """
        self.reference_path = reference_path
        self.aligner = aligner

    def can_handle(self, file_path: str) -> bool:
        """Check if file is a FASTQ."""
        return file_path.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "FASTQ"

    def extract_profile(self, file_path: str, **kwargs: object) -> "AlleleProfile":
        """Extract profile from FASTQ.

        Aligns reads to reference and extracts allele observations.

        Args:
            file_path: Path to FASTQ file

        Returns:
            AlleleProfile extracted from the FASTQ
        """
        raise NotImplementedError("FASTQAdapter.extract_profile() not yet implemented")


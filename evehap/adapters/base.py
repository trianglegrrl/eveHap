"""Base adapter interface for eveHap."""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from evehap.core.profile import AlleleProfile


class InputAdapter(ABC):
    """Abstract base class for input format adapters.

    All adapters must implement this interface to convert their
    input format to an AlleleProfile.
    """

    @abstractmethod
    def can_handle(self, file_path: str) -> bool:
        """Check if this adapter can handle the given file.

        Args:
            file_path: Path to the input file

        Returns:
            True if this adapter can read the file
        """
        pass

    @abstractmethod
    def extract_profile(self, file_path: str, **kwargs: object) -> "AlleleProfile":
        """Extract AlleleProfile from input file.

        Args:
            file_path: Path to the input file
            **kwargs: Additional adapter-specific options

        Returns:
            AlleleProfile extracted from the file
        """
        pass

    @property
    @abstractmethod
    def format_name(self) -> str:
        """Return the format name (e.g., 'BAM', 'VCF').

        Returns:
            Human-readable format name
        """
        pass


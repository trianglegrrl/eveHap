"""HSD format adapter for eveHap.

Parses Haplogrep HSD (Haplogroup Sequence Data) format files.
HSD files contain tab-delimited polymorphism data with notation like:
  - 73G (substitution at position 73 to G)
  - 315.1C (insertion after position 315)
  - 523d or 523DEL (deletion at position 523)
"""

import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from evehap.adapters.base import InputAdapter
from evehap.core.profile import AlleleObservation, AlleleProfile

if TYPE_CHECKING:
    from evehap.core.damage import DamageFilter

# Regular expressions for polymorphism parsing
# Substitution: optional ref base, position, alt base (e.g., "73G" or "A73G")
SUBST_PATTERN = re.compile(r"^([ACGT])?(\d+)([ACGT])$", re.IGNORECASE)

# Insertion: position.index, bases (e.g., "315.1C" or "309.1CC")
INSERT_PATTERN = re.compile(r"^(\d+)\.(\d+)([ACGT]+)$", re.IGNORECASE)

# Deletion: position(s) + d/del/DEL (e.g., "523d", "523DEL", "8281-8289d")
DELETE_PATTERN = re.compile(r"^(\d+)(?:-(\d+))?[dD](?:EL)?$", re.IGNORECASE)


def parse_polymorphism(poly_str: str) -> Tuple[int, Optional[str], Optional[str], str]:
    """Parse a polymorphism notation string.

    Args:
        poly_str: Polymorphism string (e.g., "73G", "315.1C", "523d")

    Returns:
        Tuple of (position, ref_allele, alt_allele, type)
        where type is "substitution", "insertion", or "deletion"

    Raises:
        ValueError: If polymorphism cannot be parsed
    """
    poly_str = poly_str.strip()

    # Try substitution pattern
    match = SUBST_PATTERN.match(poly_str)
    if match:
        ref = match.group(1).upper() if match.group(1) else None
        pos = int(match.group(2))
        alt = match.group(3).upper()
        return pos, ref, alt, "substitution"

    # Try insertion pattern
    match = INSERT_PATTERN.match(poly_str)
    if match:
        pos = int(match.group(1))
        # index = int(match.group(2))  # e.g., 1 in "315.1"
        alt = match.group(3).upper()
        return pos, None, alt, "insertion"

    # Try deletion pattern
    match = DELETE_PATTERN.match(poly_str)
    if match:
        start_pos = int(match.group(1))
        # end_pos = int(match.group(2)) if match.group(2) else start_pos
        return start_pos, None, None, "deletion"

    raise ValueError(f"Cannot parse polymorphism: {poly_str}")


class HSDAdapter(InputAdapter):
    """Adapter for Haplogrep HSD format.

    HSD files are tab-delimited with columns:
    - ID: Sample identifier
    - Range: Covered region (e.g., "1-16569")
    - Haplogroup: (optional) Known haplogroup
    - Polymorphisms: Space-separated polymorphism notations

    Polymorphisms describe differences from the rCRS reference.
    """

    def __init__(
        self,
        reference_path: Optional[str] = None,
        fill_reference: bool = False,
    ) -> None:
        """Initialize HSD adapter.

        Args:
            reference_path: Path to rCRS reference FASTA (required if fill_reference=True)
            fill_reference: Whether to fill in reference positions (default: False)
        """
        self.reference_path = reference_path
        self.fill_reference = fill_reference

    def can_handle(self, file_path: str) -> bool:
        """Check if file is an HSD file."""
        return file_path.endswith(".hsd")

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "HSD"

    def list_samples(self, file_path: str) -> List[str]:
        """List all samples in the HSD file.

        Args:
            file_path: Path to HSD file

        Returns:
            List of sample IDs
        """
        hsd_path = Path(file_path)
        if not hsd_path.exists():
            raise FileNotFoundError(f"HSD file not found: {file_path}")

        samples = []
        first_line = True
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if not parts:
                    continue

                sample_id = parts[0].strip()

                # Skip header line
                if first_line:
                    first_line = False
                    if sample_id.lower() in ("id", "sampleid", "sample_id", "sample"):
                        continue

                samples.append(sample_id)

        return samples

    def extract_profile(
        self,
        file_path: str,
        damage_filter: Optional["DamageFilter"] = None,
        sample_id: Optional[str] = None,
        **kwargs: Any,
    ) -> AlleleProfile:
        """Extract profile from HSD file.

        Parses polymorphism notation and converts to AlleleObservations.
        Optionally fills in reference positions for unspecified locations.

        Args:
            file_path: Path to HSD file
            sample_id: Sample ID to extract (None = first sample)
            damage_filter: Optional damage filter (not typically used for HSD)

        Returns:
            AlleleProfile extracted from the HSD

        Raises:
            FileNotFoundError: If HSD file doesn't exist
            ValueError: If sample not found in HSD
        """
        hsd_path = Path(file_path)
        if not hsd_path.exists():
            raise FileNotFoundError(f"HSD file not found: {file_path}")

        # Parse HSD file
        samples_data = self._parse_hsd(file_path)

        if not samples_data:
            raise ValueError(f"No samples found in HSD file: {file_path}")

        # Select target sample
        if sample_id is None:
            target_id = list(samples_data.keys())[0]
        else:
            if sample_id not in samples_data:
                available = list(samples_data.keys())[:5]
                raise ValueError(
                    f"Sample '{sample_id}' not found in HSD. " f"Available: {available}..."
                )
            target_id = sample_id

        sample_data = samples_data[target_id]

        # Create profile
        profile = AlleleProfile(
            sample_id=target_id,
            source_format=self.format_name,
        )

        # Load reference if needed
        ref_sequence: Optional[str] = None
        if self.fill_reference and self.reference_path:
            ref_sequence = self._load_reference()

        # Parse polymorphisms
        variant_positions: Dict[int, AlleleObservation] = {}

        for poly in sample_data.get("polymorphisms", []):
            try:
                pos, ref, alt, ptype = parse_polymorphism(poly)

                if ptype == "substitution" and alt:
                    # Determine reference allele
                    if ref is None and ref_sequence and 0 < pos <= len(ref_sequence):
                        ref = ref_sequence[pos - 1].upper()
                    elif ref is None:
                        ref = "N"  # Unknown reference

                    obs = AlleleObservation(
                        position=pos,
                        ref_allele=ref,
                        alleles={alt: 1.0},
                        depth=None,
                        quality=None,
                        source=self.format_name,
                    )
                    variant_positions[pos] = obs

                elif ptype == "insertion":
                    # For insertions, record at the anchor position
                    # The insertion is after this position
                    if ref_sequence and 0 < pos <= len(ref_sequence):
                        ref = ref_sequence[pos - 1].upper()
                    else:
                        ref = "N"

                    # Mark as insertion by storing the inserted bases
                    obs = AlleleObservation(
                        position=pos,
                        ref_allele=ref,
                        alleles={f"+{alt}": 1.0},  # Mark as insertion
                        depth=None,
                        quality=None,
                        source=self.format_name,
                    )
                    variant_positions[pos] = obs

                elif ptype == "deletion":
                    # Mark as deletion
                    if ref_sequence and 0 < pos <= len(ref_sequence):
                        ref = ref_sequence[pos - 1].upper()
                    else:
                        ref = "N"

                    obs = AlleleObservation(
                        position=pos,
                        ref_allele=ref,
                        alleles={"-": 1.0},  # Mark as deletion
                        depth=None,
                        quality=None,
                        source=self.format_name,
                    )
                    variant_positions[pos] = obs

            except ValueError:
                # Skip unparseable polymorphisms
                continue

        # Fill in reference positions if requested
        if self.fill_reference and ref_sequence:
            for pos_1based in range(1, len(ref_sequence) + 1):
                if pos_1based in variant_positions:
                    profile.add_observation(variant_positions[pos_1based])
                else:
                    ref_base = ref_sequence[pos_1based - 1].upper()
                    if ref_base in "ACGT":
                        obs = AlleleObservation(
                            position=pos_1based,
                            ref_allele=ref_base,
                            alleles={ref_base: 1.0},
                            depth=None,
                            quality=None,
                            source="reference",
                        )
                        profile.add_observation(obs)
        else:
            # Just add variant positions
            for obs in variant_positions.values():
                profile.add_observation(obs)

        return profile

    def _parse_hsd(self, file_path: str) -> Dict[str, Dict]:
        """Parse HSD file into dictionary.

        Args:
            file_path: Path to HSD file

        Returns:
            Dictionary mapping sample IDs to their data
        """
        samples: Dict[str, Dict] = {}

        with open(file_path) as f:
            first_line = True

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")

                # Check if this is a header line
                if first_line:
                    first_line = False
                    # Header contains "ID" or starts with "SampleId" etc.
                    first_col = parts[0].strip().lower()
                    if first_col in ("id", "sampleid", "sample_id", "sample"):
                        continue  # Skip header
                    # Otherwise treat as data

                if len(parts) < 1:
                    continue

                sample_id = parts[0].strip()

                # Parse range (column 2)
                range_str = parts[1].strip() if len(parts) > 1 else ""

                # Parse haplogroup (column 3) - may be empty
                haplogroup = parts[2].strip() if len(parts) > 2 else ""

                # Parse polymorphisms (column 4+)
                # Polymorphisms are space-separated within column 4
                polymorphisms = []
                for i in range(3, len(parts)):
                    poly_str = parts[i].strip()
                    if poly_str:
                        # Split by spaces to get individual polymorphisms
                        for poly in poly_str.split():
                            if poly:
                                polymorphisms.append(poly)

                samples[sample_id] = {
                    "range": range_str,
                    "haplogroup": haplogroup,
                    "polymorphisms": polymorphisms,
                }

        return samples

    def _load_reference(self) -> str:
        """Load reference sequence from FASTA.

        Returns:
            Reference sequence string
        """
        if not self.reference_path:
            raise ValueError("Reference path not set")

        ref_path = Path(self.reference_path)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference_path}")

        # Simple FASTA parser
        sequence_lines = []
        with open(self.reference_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                sequence_lines.append(line)

        return "".join(sequence_lines)

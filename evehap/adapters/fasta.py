"""FASTA adapter for eveHap.

Extracts mtDNA allele profiles from FASTA consensus sequences.
Compares each position to rCRS reference to identify variants.
"""

from pathlib import Path
from typing import Dict, Optional, TYPE_CHECKING

import pysam

from evehap.adapters.base import InputAdapter
from evehap.core.profile import AlleleObservation, AlleleProfile, MTDNA_LENGTH

if TYPE_CHECKING:
    from evehap.core.damage import DamageFilter

# IUPAC ambiguity codes for nucleotides
# Maps each ambiguity code to the bases it represents
IUPAC_CODES: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",  # Treat U as T
    "R": "AG",   # purine
    "Y": "CT",   # pyrimidine
    "S": "GC",   # strong
    "W": "AT",   # weak
    "K": "GT",   # keto
    "M": "AC",   # amino
    "B": "CGT",  # not A
    "D": "AGT",  # not C
    "H": "ACT",  # not G
    "V": "ACG",  # not T
    "N": "",     # any/unknown - treat as missing
    "-": "",     # gap - treat as missing
    ".": "",     # gap - treat as missing
}


class FASTAAdapter(InputAdapter):
    """Adapter for FASTA consensus sequences.

    Extracts allele profiles by comparing each position in the
    sample sequence to the rCRS reference. Handles IUPAC ambiguity
    codes as potential heteroplasmy.
    """

    def __init__(
        self,
        reference_path: str,
        sequence_index: int = 0,
    ) -> None:
        """Initialize FASTA adapter.

        Args:
            reference_path: Path to rCRS reference FASTA
            sequence_index: Which sequence to extract from multi-FASTA (default: 0)
        """
        self.reference_path = reference_path
        self.sequence_index = sequence_index

    def can_handle(self, file_path: str) -> bool:
        """Check if file is a FASTA."""
        return file_path.endswith((".fasta", ".fa", ".fna"))

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "FASTA"

    def extract_profile(
        self,
        file_path: str,
        damage_filter: Optional["DamageFilter"] = None,
        **kwargs: object,
    ) -> AlleleProfile:
        """Extract profile from FASTA.

        Compares each position to rCRS reference to identify variants.
        Handles IUPAC ambiguity codes as heteroplasmy.

        Args:
            file_path: Path to FASTA file
            damage_filter: Optional damage filter (not typically used for FASTA)

        Returns:
            AlleleProfile extracted from the FASTA

        Raises:
            FileNotFoundError: If FASTA or reference file doesn't exist
            ValueError: If FASTA is empty or invalid
        """
        fasta_path = Path(file_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

        ref_path = Path(self.reference_path)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference_path}")

        # Load reference sequence
        ref_sequence = self._load_reference()

        # Load sample sequence
        sample_id, sample_sequence = self._load_sample_sequence(file_path)

        # Create profile
        profile = AlleleProfile(
            sample_id=sample_id,
            source_format=self.format_name,
        )

        # Compare each position
        max_len = min(len(ref_sequence), len(sample_sequence))

        for pos_0based in range(max_len):
            pos_1based = pos_0based + 1

            ref_base = ref_sequence[pos_0based].upper()
            sample_base = sample_sequence[pos_0based].upper()

            # Skip if reference is not a standard base
            if ref_base not in "ACGT":
                continue

            # Handle deletion marker from bcftools consensus --mark-del
            if sample_base == "-":
                # This is a deletion at this position
                alleles = {"-": 1.0}
            else:
                # Resolve IUPAC code
                resolved_bases = IUPAC_CODES.get(sample_base, "")

                if not resolved_bases:
                    # Unknown base - skip this position
                    continue

                # Create observation
                if len(resolved_bases) == 1:
                    # Unambiguous base
                    alleles = {resolved_bases: 1.0}
                else:
                    # Ambiguous - split equally
                    freq = 1.0 / len(resolved_bases)
                    alleles = {base: freq for base in resolved_bases}

            obs = AlleleObservation(
                position=pos_1based,
                ref_allele=ref_base,
                alleles=alleles,
                depth=None,  # FASTA has no depth info
                quality=None,
                source=self.format_name,
            )

            profile.add_observation(obs)

        return profile

    def _load_reference(self) -> str:
        """Load reference sequence from FASTA.

        Returns:
            Reference sequence string
        """
        sequences = self._parse_fasta(self.reference_path)
        if not sequences:
            raise ValueError(f"No sequences found in reference: {self.reference_path}")

        # Return first sequence
        return list(sequences.values())[0]

    def _load_sample_sequence(self, file_path: str) -> tuple:
        """Load sample sequence from FASTA.

        Args:
            file_path: Path to FASTA file

        Returns:
            Tuple of (sample_id, sequence)
        """
        sequences = self._parse_fasta(file_path)
        if not sequences:
            raise ValueError(f"No sequences found in: {file_path}")

        # Get sequence by index
        sample_ids = list(sequences.keys())
        if self.sequence_index >= len(sample_ids):
            sample_id = sample_ids[0]
        else:
            sample_id = sample_ids[self.sequence_index]

        return sample_id, sequences[sample_id]

    def _parse_fasta(self, file_path: str) -> Dict[str, str]:
        """Parse FASTA file into dictionary.

        Args:
            file_path: Path to FASTA file

        Returns:
            Dictionary mapping sequence IDs to sequences
        """
        sequences: Dict[str, str] = {}
        current_id: Optional[str] = None
        current_seq: list = []

        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    # New sequence
                    if current_id is not None:
                        sequences[current_id] = "".join(current_seq)

                    # Parse header - take first word after >
                    header = line[1:].strip()
                    current_id = header.split()[0] if header else "unknown"
                    current_seq = []
                else:
                    # Sequence line
                    current_seq.append(line)

        # Don't forget last sequence
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

        return sequences

    def list_sequences(self, file_path: str) -> list:
        """List all sequence IDs in a FASTA file.

        Args:
            file_path: Path to FASTA file

        Returns:
            List of sequence IDs
        """
        fasta_path = Path(file_path)
        if not fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

        sequences = self._parse_fasta(file_path)
        return list(sequences.keys())


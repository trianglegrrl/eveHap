"""FASTA adapter for eveHap.

Extracts mtDNA allele profiles from FASTA consensus sequences.
Uses alignment-aware comparison to handle insertions/deletions
relative to the rCRS reference.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

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

    Extracts allele profiles by aligning sample sequences to the rCRS
    reference and extracting variants from the alignment. This properly
    handles insertions and deletions in the sample sequence.

    Uses BioPython's PairwiseAligner for global alignment when needed,
    with a fast path for sequences that are already aligned (exact length
    match with few differences).
    """

    # Threshold for naive comparison: if fewer differences than this,
    # assume the sequence is pre-aligned and skip alignment
    NAIVE_DIFF_THRESHOLD = 100

    def __init__(
        self,
        reference_path: str,
        sequence_index: int = 0,
        force_alignment: bool = False,
    ) -> None:
        """Initialize FASTA adapter.

        Args:
            reference_path: Path to rCRS reference FASTA
            sequence_index: Which sequence to extract from multi-FASTA (default: 0)
            force_alignment: If True, always use alignment even for exact-length sequences
        """
        self.reference_path = reference_path
        self.sequence_index = sequence_index
        self.force_alignment = force_alignment
        self._ref_sequence: Optional[str] = None
        self._aligner = None

    def can_handle(self, file_path: str) -> bool:
        """Check if file is a FASTA."""
        return file_path.endswith((".fasta", ".fa", ".fna"))

    @property
    def format_name(self) -> str:
        """Return format name."""
        return "FASTA"

    def _get_aligner(self):
        """Get or create the pairwise aligner.

        Uses settings optimized for mtDNA comparison:
        - Global alignment (end-to-end)
        - Strong gap penalties (indels are rare in mtDNA)
        """
        if self._aligner is None:
            from Bio.Align import PairwiseAligner
            self._aligner = PairwiseAligner()
            # Global alignment
            self._aligner.mode = 'global'
            # Scoring: favor matches, penalize gaps heavily
            self._aligner.match_score = 2
            self._aligner.mismatch_score = -1
            self._aligner.open_gap_score = -10
            self._aligner.extend_gap_score = -0.5
        return self._aligner

    def _get_reference(self) -> str:
        """Get reference sequence, loading if needed."""
        if self._ref_sequence is None:
            self._ref_sequence = self._load_reference()
        return self._ref_sequence

    def _count_naive_differences(self, ref: str, sample: str) -> int:
        """Count differences using naive position-by-position comparison.

        Used to determine if alignment is needed.
        """
        max_len = min(len(ref), len(sample))
        return sum(1 for i in range(max_len) if ref[i].upper() != sample[i].upper())

    def _needs_alignment(self, ref_sequence: str, sample_sequence: str) -> bool:
        """Determine if alignment is needed.

        Returns False (skip alignment) if:
        - force_alignment is False
        - Sequence lengths match exactly
        - Naive comparison shows few differences

        This provides a fast path for pre-aligned sequences.
        """
        if self.force_alignment:
            return True

        # If lengths differ, alignment is definitely needed
        if len(sample_sequence) != len(ref_sequence):
            return True

        # For exact-length sequences, check naive difference count
        naive_diffs = self._count_naive_differences(ref_sequence, sample_sequence)
        return naive_diffs > self.NAIVE_DIFF_THRESHOLD

    def _align_sequences(self, ref: str, sample: str) -> Tuple[str, str]:
        """Align sample to reference using global pairwise alignment.

        Args:
            ref: Reference sequence
            sample: Sample sequence

        Returns:
            Tuple of (aligned_ref, aligned_sample) with gaps inserted
        """
        aligner = self._get_aligner()

        # Get best alignment
        alignments = aligner.align(ref.upper(), sample.upper())
        if not alignments:
            raise ValueError("Failed to align sequences")

        # Get the best alignment
        best = alignments[0]

        # Extract aligned sequences
        # The alignment object can be indexed to get the aligned sequences
        aligned_ref = ""
        aligned_sample = ""

        # Walk through alignment blocks
        ref_pos = 0
        sample_pos = 0

        for (ref_start, ref_end), (sample_start, sample_end) in zip(
            best.aligned[0], best.aligned[1]
        ):
            # Add gaps for any unaligned prefix in reference
            while ref_pos < ref_start:
                aligned_ref += ref[ref_pos]
                aligned_sample += "-"
                ref_pos += 1

            # Add gaps for any unaligned prefix in sample
            while sample_pos < sample_start:
                aligned_ref += "-"
                aligned_sample += sample[sample_pos]
                sample_pos += 1

            # Add aligned block
            block_len = ref_end - ref_start
            aligned_ref += ref[ref_start:ref_end]
            aligned_sample += sample[sample_start:sample_end]
            ref_pos = ref_end
            sample_pos = sample_end

        # Add any trailing unaligned portions
        while ref_pos < len(ref):
            aligned_ref += ref[ref_pos]
            aligned_sample += "-"
            ref_pos += 1

        while sample_pos < len(sample):
            aligned_ref += "-"
            aligned_sample += sample[sample_pos]
            sample_pos += 1

        return aligned_ref, aligned_sample

    def _extract_variants_from_alignment(
        self,
        aligned_ref: str,
        aligned_sample: str,
        sample_id: str,
    ) -> AlleleProfile:
        """Extract variants from aligned sequences.

        Walks through the alignment tracking reference positions,
        and creates observations for each covered reference position.

        Args:
            aligned_ref: Aligned reference sequence (with gaps)
            aligned_sample: Aligned sample sequence (with gaps)
            sample_id: Sample identifier

        Returns:
            AlleleProfile with observations at reference positions
        """
        profile = AlleleProfile(
            sample_id=sample_id,
            source_format=self.format_name,
        )

        ref_pos = 0  # 0-based reference position

        for i in range(len(aligned_ref)):
            ref_base = aligned_ref[i].upper()
            sample_base = aligned_sample[i].upper()

            # If reference has a gap, this is an insertion in sample
            # We don't track insertions as observations (they don't map to ref positions)
            if ref_base == "-":
                continue

            # Increment reference position (1-based for output)
            ref_pos += 1
            pos_1based = ref_pos

            # Skip if reference is not a standard base
            if ref_base not in "ACGT":
                continue

            # If sample has a gap, this is a deletion
            if sample_base == "-":
                alleles = {"-": 1.0}
            else:
                # Resolve IUPAC code
                resolved_bases = IUPAC_CODES.get(sample_base, "")

                if not resolved_bases:
                    # Unknown base - skip this position
                    continue

                # Create observation
                if len(resolved_bases) == 1:
                    alleles = {resolved_bases: 1.0}
                else:
                    # Ambiguous - split equally
                    freq = 1.0 / len(resolved_bases)
                    alleles = {base: freq for base in resolved_bases}

            obs = AlleleObservation(
                position=pos_1based,
                ref_allele=ref_base,
                alleles=alleles,
                depth=None,
                quality=None,
                source=self.format_name,
            )
            profile.add_observation(obs)

        return profile

    def _extract_variants_naive(
        self,
        ref_sequence: str,
        sample_sequence: str,
        sample_id: str,
    ) -> AlleleProfile:
        """Extract variants using naive position-by-position comparison.

        Used as fast path when sequences are already aligned (same length,
        few differences).

        Args:
            ref_sequence: Reference sequence
            sample_sequence: Sample sequence
            sample_id: Sample identifier

        Returns:
            AlleleProfile with observations
        """
        profile = AlleleProfile(
            sample_id=sample_id,
            source_format=self.format_name,
        )

        max_len = min(len(ref_sequence), len(sample_sequence))

        for pos_0based in range(max_len):
            pos_1based = pos_0based + 1

            ref_base = ref_sequence[pos_0based].upper()
            sample_base = sample_sequence[pos_0based].upper()

            # Skip if reference is not a standard base
            if ref_base not in "ACGT":
                continue

            # Handle deletion marker
            if sample_base == "-":
                alleles = {"-": 1.0}
            else:
                # Resolve IUPAC code
                resolved_bases = IUPAC_CODES.get(sample_base, "")

                if not resolved_bases:
                    continue

                if len(resolved_bases) == 1:
                    alleles = {resolved_bases: 1.0}
                else:
                    freq = 1.0 / len(resolved_bases)
                    alleles = {base: freq for base in resolved_bases}

            obs = AlleleObservation(
                position=pos_1based,
                ref_allele=ref_base,
                alleles=alleles,
                depth=None,
                quality=None,
                source=self.format_name,
            )
            profile.add_observation(obs)

        return profile

    def extract_profile(
        self,
        file_path: str,
        damage_filter: Optional["DamageFilter"] = None,
        **kwargs: object,
    ) -> AlleleProfile:
        """Extract profile from FASTA.

        Uses alignment-aware comparison to properly handle insertions
        and deletions in the sample sequence. Falls back to naive
        comparison for pre-aligned sequences (exact length match with
        few differences) for efficiency.

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

        # Load sequences
        ref_sequence = self._get_reference()
        sample_id, sample_sequence = self._load_sample_sequence(file_path)

        # Decide whether to use alignment
        if self._needs_alignment(ref_sequence, sample_sequence):
            # Use alignment-aware extraction
            aligned_ref, aligned_sample = self._align_sequences(
                ref_sequence, sample_sequence
            )
            return self._extract_variants_from_alignment(
                aligned_ref, aligned_sample, sample_id
            )
        else:
            # Use fast naive extraction
            return self._extract_variants_naive(
                ref_sequence, sample_sequence, sample_id
            )

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

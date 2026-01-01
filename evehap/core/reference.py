"""Reference sequence handling for eveHap.

Handles translation between rCRS and RSRS reference coordinates.

rCRS (revised Cambridge Reference Sequence) = H2a2a1 haplogroup
RSRS (Reconstructed Sapiens Reference Sequence) = ancestral mtDNA

Most sequencing data is aligned to rCRS, but the RSRS-based phylotree
has the complete structure. This module provides translation between them.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple

# Package data directory
PACKAGE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = PACKAGE_DIR / "data"

# Default reference paths
RCRS_PATH = DATA_DIR / "reference" / "rCRS.fasta"
RSRS_PATH = DATA_DIR / "reference" / "RSRS.fasta"

# Cached sequences
_RCRS_SEQ: Optional[str] = None
_RSRS_SEQ: Optional[str] = None
_DIFF_POSITIONS: Optional[Dict[int, Tuple[str, str]]] = None


def _load_fasta(path: Path) -> str:
    """Load sequence from FASTA file."""
    if not path.exists():
        return ""
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def get_rcrs_sequence() -> str:
    """Get rCRS reference sequence (cached)."""
    global _RCRS_SEQ
    if _RCRS_SEQ is None:
        _RCRS_SEQ = _load_fasta(RCRS_PATH)
    return _RCRS_SEQ


def get_rsrs_sequence() -> str:
    """Get RSRS reference sequence (cached)."""
    global _RSRS_SEQ
    if _RSRS_SEQ is None:
        _RSRS_SEQ = _load_fasta(RSRS_PATH)
    return _RSRS_SEQ


def get_reference_differences() -> Dict[int, Tuple[str, str]]:
    """Get positions where rCRS and RSRS differ.

    Returns:
        Dictionary mapping position (1-based) to (RSRS_allele, rCRS_allele)
    """
    global _DIFF_POSITIONS
    if _DIFF_POSITIONS is None:
        rcrs = get_rcrs_sequence()
        rsrs = get_rsrs_sequence()

        _DIFF_POSITIONS = {}
        for i in range(min(len(rcrs), len(rsrs))):
            r = rcrs[i].upper()
            s = rsrs[i].upper()
            if r != s and r in "ACGT" and s in "ACGT":
                _DIFF_POSITIONS[i + 1] = (s, r)  # (RSRS, rCRS)

    return _DIFF_POSITIONS


def rcrs_to_rsrs_allele(position: int, rcrs_allele: str) -> str:
    """Convert an allele from rCRS-relative to RSRS-relative.

    At positions where rCRS and RSRS differ:
    - If sample has rCRS allele -> this is derived relative to RSRS
    - If sample differs from rCRS -> check if it matches RSRS (ancestral)

    Args:
        position: 1-based mtDNA position
        rcrs_allele: Allele observed in sample (relative to rCRS alignment)

    Returns:
        The allele (unchanged - the allele is the allele)

    Note: This function is mainly for documentation. The actual allele
    doesn't change - what changes is the interpretation of ref vs alt.
    """
    return rcrs_allele.upper()


def get_rsrs_allele_at(position: int) -> str:
    """Get RSRS reference allele at position.

    Args:
        position: 1-based mtDNA position

    Returns:
        RSRS allele at position, or 'N' if unavailable
    """
    rsrs = get_rsrs_sequence()
    if rsrs and 0 < position <= len(rsrs):
        return rsrs[position - 1].upper()
    return "N"


def get_rcrs_allele_at(position: int) -> str:
    """Get rCRS reference allele at position.

    Args:
        position: 1-based mtDNA position

    Returns:
        rCRS allele at position, or 'N' if unavailable
    """
    rcrs = get_rcrs_sequence()
    if rcrs and 0 < position <= len(rcrs):
        return rcrs[position - 1].upper()
    return "N"


def translate_profile_reference(
    profile: "AlleleProfile",
    from_ref: str = "rcrs",
    to_ref: str = "rsrs",
) -> "AlleleProfile":
    """Translate profile reference alleles from one reference to another.

    This updates the ref_allele field in each observation to reflect
    the new reference sequence, without changing the observed alleles.

    Args:
        profile: AlleleProfile to translate
        from_ref: Source reference ('rcrs' or 'rsrs')
        to_ref: Target reference ('rcrs' or 'rsrs')

    Returns:
        New AlleleProfile with updated reference alleles
    """
    from evehap.core.profile import AlleleProfile, AlleleObservation

    if from_ref == to_ref:
        return profile

    # Create new profile
    new_profile = AlleleProfile(
        sample_id=profile.sample_id,
        source_format=profile.source_format,
    )
    new_profile.metadata = profile.metadata.copy()
    new_profile.metadata["original_reference"] = from_ref
    new_profile.metadata["current_reference"] = to_ref

    # Get reference getter for target
    if to_ref == "rsrs":
        get_ref = get_rsrs_allele_at
    else:
        get_ref = get_rcrs_allele_at

    # Translate observations
    for pos, obs in profile.observations.items():
        new_ref = get_ref(pos)

        new_obs = AlleleObservation(
            position=obs.position,
            ref_allele=new_ref,
            alleles=obs.alleles.copy(),
            depth=obs.depth,
            quality=obs.quality,
            source=obs.source,
        )
        new_profile.add_observation(new_obs)

    return new_profile



import pandas as pd
from typing import Tuple
from src.config.logging import get_logger

logger = get_logger(__name__)

from src.utils.classification_utils import (check_exact_match, check_contained_match, check_similar_match)

"""
Submission workflow
- first a sequence is run through the classification workflow
- this workflow inherits the dataclasses from the classification workflow
- use this information to do the following:
1. compare with existing ships
2. check 
assign accession number
2.

"""

def assign_accession(
    sequence: str, existing_ships: pd.DataFrame = None, threshold: float = 0.95
) -> Tuple[str, bool]:
    """Assign an accession to a new sequence.

    Args:
        sequence: New sequence to assign accession to
        existing_ships: DataFrame of existing ships (optional, will fetch if None)
        threshold: Similarity threshold for "almost identical" matches

    Returns:
        Tuple[str, bool]: (accession, needs_review)
            - accession: assigned accession number
            - needs_review: True if sequence needs manual review
    """

    logger.debug("Starting accession assignment process")

    logger.debug("Step 1: Checking for exact matches using MD5 hash...")
    exact_match = check_exact_match(sequence, existing_ships)
    if exact_match:
        logger.debug(f"Found exact match: {exact_match}")
        return exact_match, False

    logger.debug("Step 2: Checking for contained matches...")
    container_match = check_contained_match(
        fasta=sequence,
        existing_ships=existing_ships,
        min_coverage=0.95,
        min_identity=0.95,
    )
    if container_match:
        logger.debug(f"Found containing match: {container_match}")
        return container_match, True  # Flag for review since it's truncated

    logger.debug(f"Step 3: Checking for similar matches (threshold={threshold})...")
    similar_match = check_similar_match(sequence, existing_ships, threshold)
    if similar_match:
        logger.debug(f"Found similar match: {similar_match}")
        return similar_match, True  # Flag for review due to high similarity

    logger.debug("No matches found, generating new accession...")
    new_accession = generate_new_accession(existing_ships)
    logger.debug(f"Generated new accession: {new_accession}")
    return new_accession, False

def generate_new_accession(existing_ships: pd.DataFrame) -> str:
    """Generate a new unique accession number."""
    # Extract existing accession numbers
    existing_nums = [
        int(acc.replace("SBS", "").split(".")[0])
        for acc in existing_ships["accession_tag"]
        if acc.startswith("SBS")
    ]

    # Check if we have existing accessions
    if not existing_nums:
        error_msg = "Problem with loading existing ships. No existing SBS accessions found in database."
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Find next available number
    next_num = max(existing_nums) + 1
    logger.debug(f"Last used accession number: SBS{max(existing_nums):06d}")
    logger.debug(f"Assigning new accession number: SBS{next_num:06d}")

    return f"SBS{next_num:06d}"

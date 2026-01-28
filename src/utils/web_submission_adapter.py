#!/usr/bin/env python3
"""
Web submission adapter for Starbase submission portal.

This module provides web-friendly adapter functions that bridge the gap between
the Dash web form and the comprehensive submission_utils.py processor.

It handles:
- Converting web form data to the submission_utils schema
- Providing user-friendly validation messages
- Duplicate detection with helpful feedback
- Database insertion with proper error handling
"""

import base64
import datetime
from typing import Dict, Any, Optional
from Bio import SeqIO
from io import StringIO

from src.config.logging import get_logger
from src.utils.submission_utils import (
    SubmissionProcessor,
    validate_submission,
    check_sequence_duplicate,
)

logger = get_logger(__name__)


class WebValidationError(Exception):
    """Web-friendly validation error."""

    def __init__(self, message: str, field: str = None):
        self.message = message
        self.field = field
        super().__init__(self.message)


def parse_fasta_sequences(contents: str, filename: str) -> Dict[str, str]:
    """
    Parse FASTA sequences from uploaded file contents.

    Args:
        contents: Base64-encoded file contents from Dash upload
        filename: Original filename

    Returns:
        Dict mapping sequence IDs to sequences

    Raises:
        WebValidationError: If FASTA parsing fails
    """
    try:
        # Decode base64 content
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string).decode("utf-8")

        # Parse FASTA
        sequences = {}
        fasta_io = StringIO(decoded)

        for record in SeqIO.parse(fasta_io, "fasta"):
            sequences[record.id] = str(record.seq)

        if not sequences:
            raise WebValidationError(
                "No valid sequences found in FASTA file", "fasta_file"
            )

        logger.info(f"Parsed {len(sequences)} sequences from {filename}")
        return sequences

    except Exception as e:
        logger.error(f"Error parsing FASTA: {str(e)}")
        raise WebValidationError(f"Failed to parse FASTA file: {str(e)}", "fasta_file")


def validate_submission_data(
    seq_contents: str,
    seq_filename: str,
    uploader: str,
    evidence: str,
    genus: str,
    species: str,
    hostchr: str,
    shipstart: int,
    shipend: int,
) -> Dict[str, Any]:
    """
    Validate web form submission data.

    Args:
        seq_contents: Base64-encoded FASTA file contents
        seq_filename: FASTA filename
        uploader: Email of curator
        evidence: Evidence/method used
        genus: Genus name
        species: Species name
        hostchr: Host chromosome/contig ID
        shipstart: Start coordinate
        shipend: End coordinate

    Returns:
        Dict with validated data

    Raises:
        WebValidationError: If validation fails
    """
    errors = []

    # Check required fields
    if not seq_contents:
        errors.append("FASTA file is required")
    if not seq_filename:
        errors.append("FASTA filename is required")
    if not uploader:
        errors.append("Curator email is required")
    if not evidence:
        errors.append("Evidence/annotation method is required")
    if not genus:
        errors.append("Genus is required")
    if not species:
        errors.append("Species is required")
    if not hostchr:
        errors.append("Host chromosome/contig ID is required")
    if shipstart is None:
        errors.append("Start coordinate is required")
    if shipend is None:
        errors.append("End coordinate is required")

    if errors:
        raise WebValidationError("; ".join(errors))

    # Validate coordinates
    if shipstart <= 0:
        errors.append("Start coordinate must be greater than 0")
    if shipend <= 0:
        errors.append("End coordinate must be greater than 0")
    if shipstart == shipend:
        errors.append("Start and end coordinates cannot be the same")

    if errors:
        raise WebValidationError("; ".join(errors))

    # Parse FASTA
    try:
        sequences = parse_fasta_sequences(seq_contents, seq_filename)

        # For now, we expect single sequence per submission
        if len(sequences) > 1:
            logger.warning(
                f"Multiple sequences found in {seq_filename}, using first one"
            )

        # Get first sequence
        seq_id, sequence = next(iter(sequences.items()))

    except WebValidationError:
        raise
    except Exception as e:
        raise WebValidationError(f"Error processing FASTA file: {str(e)}", "fasta_file")

    # Return validated data
    return {
        "sequence": sequence,
        "seq_id": seq_id,
        "filename": seq_filename,
        "uploader": uploader,
        "evidence": evidence,
        "genus": genus,
        "species": species,
        "hostchr": hostchr,
        "shipstart": int(shipstart),
        "shipend": int(shipend),
    }


def process_submission_data(
    validated_data: Dict[str, Any], strand_radio: int
) -> Dict[str, Any]:
    """
    Process validated submission data and prepare for database insertion.

    Args:
        validated_data: Output from validate_submission_data()
        strand_radio: Strand selection (1=forward, 2=reverse)

    Returns:
        Dict ready for database insertion
    """
    # Convert strand radio to strand symbol
    strand = "+" if strand_radio == 1 else "-"

    # Normalize coordinates (ensure start < end)
    start = validated_data["shipstart"]
    end = validated_data["shipend"]

    if start > end:
        # Swap for storage
        start, end = end, start
        # If user provided reversed coordinates, they likely meant reverse strand
        if strand == "+":
            strand = "-"
            logger.info("Auto-detected reverse strand from coordinate order")

    # Generate starshipID from filename
    starship_id = validated_data["seq_id"]

    # Check for duplicates
    duplicate_info = check_sequence_duplicate(
        validated_data["sequence"], validated_data["genus"], validated_data["species"]
    )

    # Prepare data structure matching submission_utils schema
    processed_data = {
        "sequence": validated_data["sequence"],
        "starshipID": starship_id,
        "evidence": validated_data["evidence"],
        "source": f"web_submission_{datetime.datetime.now().strftime('%Y%m%d')}",
        "genus": validated_data["genus"],
        "species": validated_data["species"],
        "strain": None,
        "contig_id": validated_data["hostchr"],
        "element_start": start,
        "element_end": end,
        "element_strand": strand,
        "curator": validated_data["uploader"],
        "curated_status": "needs_review",  # Default for web submissions
        "notes": validated_data.get("comment", ""),
        "duplicate_info": duplicate_info,
    }

    # Validate using submission_utils validator
    is_valid, validation_errors = validate_submission(processed_data)

    if not is_valid:
        raise WebValidationError("; ".join(validation_errors))

    return processed_data


def perform_database_insertion(
    processed_data: Dict[str, Any],
    anno_contents: Optional[str],
    anno_filename: Optional[str],
    anno_date: Optional[float],
    seq_date: float,
) -> Dict[str, Any]:
    """
    Insert processed submission into database.

    Args:
        processed_data: Output from process_submission_data()
        anno_contents: Optional GFF file contents (base64)
        anno_filename: Optional GFF filename
        anno_date: Optional GFF upload timestamp
        seq_date: FASTA upload timestamp

    Returns:
        Dict with insertion results

    Raises:
        WebValidationError: If insertion fails
    """
    try:
        # Check for duplicates
        duplicate_info = processed_data.get("duplicate_info")

        if duplicate_info and duplicate_info.is_duplicate:
            if not duplicate_info.different_taxon:
                # Exact duplicate in same taxon
                raise WebValidationError(
                    f"This sequence already exists in the database "
                    f"(Accession: {duplicate_info.existing_accession}). "
                    f"Duplicate submissions from the same organism are not allowed.",
                    "sequence",
                )
            else:
                # Same sequence, different taxon - allow with warning
                logger.info(
                    f"Duplicate sequence from different taxon - creating new entry "
                    f"(existing: {duplicate_info.existing_accession})"
                )

        # Parse GFF if provided
        gff_entries = None
        if anno_contents and anno_filename:
            try:
                gff_entries = _parse_gff_contents(anno_contents, anno_filename)
                processed_data["gff_entries"] = gff_entries
            except Exception as e:
                logger.warning(f"Error parsing GFF file: {str(e)}")
                # Don't fail submission if GFF parsing fails

        # Use SubmissionProcessor for database insertion
        processor = SubmissionProcessor(dry_run=False)
        result = processor.process_submission(processed_data)

        if not result["success"]:
            error_msg = "; ".join(result.get("errors", ["Unknown error"]))
            raise WebValidationError(f"Database insertion failed: {error_msg}")

        # Determine if needs review
        needs_review = processed_data.get("curated_status") == "needs_review"

        return {
            "success": True,
            "ship_id": result["ship_id"],
            "accession": result["accession"],
            "needs_review": needs_review,
            "filename": processed_data.get("starshipID"),
            "uploader": processed_data.get("curator"),
            "warnings": result.get("warnings", []),
        }

    except WebValidationError:
        raise
    except Exception as e:
        logger.error(f"Database insertion error: {str(e)}", exc_info=True)
        raise WebValidationError(f"Failed to insert into database: {str(e)}")


def _parse_gff_contents(anno_contents: str, anno_filename: str) -> list:
    """
    Parse GFF file contents into list of entry dicts.

    Args:
        anno_contents: Base64-encoded GFF file contents
        anno_filename: GFF filename

    Returns:
        List of GFF entry dicts
    """
    try:
        # Decode base64
        content_type, content_string = anno_contents.split(",")
        decoded = base64.b64decode(content_string).decode("utf-8")

        gff_entries = []

        for line in decoded.split("\n"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue

            entry = {
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5] if parts[5] != "." else None,
                "strand": parts[6],
                "phase": parts[7] if parts[7] != "." else None,
                "attributes": parts[8] if len(parts) > 8 else None,
            }
            gff_entries.append(entry)

        logger.info(f"Parsed {len(gff_entries)} GFF entries from {anno_filename}")
        return gff_entries

    except Exception as e:
        logger.error(f"Error parsing GFF: {str(e)}")
        raise


def check_duplicate_sequence(sequence: str, genus: str, species: str) -> Dict[str, Any]:
    """
    Check if sequence is a duplicate (web-friendly wrapper).

    Args:
        sequence: DNA sequence
        genus: Genus name
        species: Species name

    Returns:
        Dict with duplicate information
    """
    duplicate_info = check_sequence_duplicate(sequence, genus, species)

    return {
        "is_duplicate": duplicate_info.is_duplicate,
        "existing_accession": duplicate_info.existing_accession,
        "existing_ship_id": duplicate_info.existing_ship_id,
        "different_taxon": duplicate_info.different_taxon,
        "match_type": duplicate_info.match_type,
    }

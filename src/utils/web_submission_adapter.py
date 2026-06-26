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
from typing import Dict, Any, List, Optional, Set
from Bio import SeqIO
from io import StringIO

import pandas as pd

from src.config.logging import get_logger
from src.utils.submission_utils import (
    SubmissionProcessor,
    validate_submission,
    check_sequence_duplicate,
)
from src.utils.seq_utils import clean_sequence, revcomp
from src.utils.classification_utils import assign_accession, generate_md5_hash
from src.database.sql_manager import fetch_ships
from src.database.sql_engine import get_submissions_session
from sqlalchemy import text

logger = get_logger(__name__)

MAX_SHIPS_PER_SUBMISSION = 100
MAX_SEQUENCES_PER_FASTA = 100


class WebValidationError(Exception):
    """Web-friendly validation error."""

    def __init__(self, message: str, field: str = None):
        self.message = message
        self.field = field
        super().__init__(self.message)


def decode_upload_contents(contents: str) -> str:
    """Decode base64 Dash upload contents to plain text."""
    if "," in str(contents) and "base64" in str(contents):
        _ct, content_string = str(contents).split(",", 1)
        return base64.b64decode(content_string).decode("utf-8")
    return str(contents)


def encode_upload_contents(text: str, mime: str = "text/plain") -> str:
    """Encode plain text for Dash upload storage."""
    encoded = base64.b64encode(text.encode("utf-8")).decode("utf-8")
    return f"data:{mime};base64,{encoded}"


def parse_fasta_records(contents: str, filename: str) -> List[Dict[str, str]]:
    """
    Parse FASTA sequences from uploaded file contents.

    Returns:
        Ordered list of {"id": header, "sequence": sequence_string}
    """
    try:
        decoded = decode_upload_contents(contents)
        records = []
        for record in SeqIO.parse(StringIO(decoded), "fasta"):
            records.append({"id": record.id, "sequence": str(record.seq)})

        if not records:
            raise WebValidationError(
                "No valid sequences found in FASTA file", "fasta_file"
            )
        if len(records) > MAX_SEQUENCES_PER_FASTA:
            raise WebValidationError(
                f"FASTA file contains {len(records)} sequences. "
                f"Maximum {MAX_SEQUENCES_PER_FASTA} sequences per file.",
                "fasta_file",
            )

        logger.info(f"Parsed {len(records)} sequences from {filename}")
        return records

    except WebValidationError:
        raise
    except Exception as e:
        logger.error(f"Error parsing FASTA: {str(e)}")
        raise WebValidationError(f"Failed to parse FASTA file: {str(e)}", "fasta_file")


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
    records = parse_fasta_records(contents, filename)
    return {record["id"]: record["sequence"] for record in records}


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
        records = parse_fasta_records(seq_contents, seq_filename)
        if len(records) > 1:
            raise WebValidationError(
                f"FASTA file contains {len(records)} sequences. "
                "Provide location details for each sequence below, or upload one sequence per file.",
                "fasta_file",
            )

        record = records[0]
        seq_id, sequence = record["id"], record["sequence"]

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


def get_gff_seqids(anno_contents: str) -> Set[str]:
    """Return unique seqids referenced in a GFF upload."""
    seqids: Set[str] = set()
    if not anno_contents:
        return seqids

    for line in decode_upload_contents(anno_contents).split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts:
            seqids.add(parts[0])
    return seqids


def filter_gff_for_seqid(anno_contents: Optional[str], seq_id: str) -> Optional[str]:
    """Return GFF upload contents containing only rows for the given seqid."""
    if not anno_contents:
        return None

    header_lines = []
    matching_rows = []
    for line in decode_upload_contents(anno_contents).split("\n"):
        if not line.strip():
            continue
        if line.startswith("#"):
            header_lines.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 8 and parts[0] == seq_id:
            matching_rows.append(line)

    if not matching_rows:
        return None

    filtered_text = "\n".join(header_lines + matching_rows) + "\n"
    return encode_upload_contents(filtered_text)


def encode_single_fasta(seq_id: str, sequence: str) -> str:
    """Encode a single-sequence FASTA for per-ship storage."""
    return encode_upload_contents(f">{seq_id}\n{sequence}\n")


def validate_location_fields(
    hostchr: str,
    shipstart: int,
    shipend: int,
    label: str = "Starship",
) -> None:
    """Validate host location fields for one sequence."""
    errors = []
    if not hostchr:
        errors.append(f"{label}: host contig or scaffold ID is required")
    if shipstart is None:
        errors.append(f"{label}: start coordinate is required")
    if shipend is None:
        errors.append(f"{label}: end coordinate is required")
    if shipstart is not None and shipstart <= 0:
        errors.append(f"{label}: start coordinate must be greater than 0")
    if shipend is not None and shipend <= 0:
        errors.append(f"{label}: end coordinate must be greater than 0")
    if shipstart is not None and shipend is not None and shipstart == shipend:
        errors.append(f"{label}: start and end coordinates cannot be the same")

    if errors:
        raise WebValidationError("; ".join(errors))


def validate_organism_fields(
    genus: str,
    species: str,
    label: str = "Starship",
) -> None:
    """Validate organism fields for one sequence."""
    errors = []
    if not genus:
        errors.append(f"{label}: genus is required")
    if not species:
        errors.append(f"{label}: species is required")
    if errors:
        raise WebValidationError("; ".join(errors))


def _normalize_metadata_strand(value: Any) -> str:
    normalized = str(value).strip().lower()
    if normalized in ("+", "1", "plus"):
        return "+"
    if normalized in ("-", "2", "minus"):
        return "-"
    raise WebValidationError(
        f"Invalid strand value '{value}'. Use +, -, 1, 2, plus, or minus."
    )


def _parse_metadata_coordinate(value: Any) -> Optional[int]:
    if value is None or str(value).strip() == "":
        return None
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def infer_metadata_strand(
    strand_value: Any,
    shipstart: Any = None,
    shipend: Any = None,
) -> str:
    """
    Resolve strand as '+' or '-'.

    Uses explicit strand when provided; otherwise infers reverse strand when
    start > end (same rule as process_submission_data).
    """
    if strand_value is not None and str(strand_value).strip() != "":
        return _normalize_metadata_strand(strand_value)

    start = _parse_metadata_coordinate(shipstart)
    end = _parse_metadata_coordinate(shipend)
    if start is not None and end is not None and start > end:
        return "-"
    return "+"


METADATA_COLUMNS = {
    "ship_header",
    "genus",
    "species",
    "strain",
    "assembly_accession",
    "hostchr",
    "shipstart",
    "shipend",
    "strand",
}


def _validate_metadata_columns(raw_columns: List[str]) -> List[str]:
    """Require exact column names; only ship_header is mandatory."""
    columns = [str(col).strip() for col in raw_columns if str(col).strip()]
    unknown = set(columns) - METADATA_COLUMNS
    if unknown:
        raise WebValidationError(
            "Unknown metadata column(s): "
            + ", ".join(sorted(unknown))
            + ". Allowed columns: "
            + ", ".join(sorted(METADATA_COLUMNS))
        )
    if "ship_header" not in columns:
        raise WebValidationError("Metadata file must include a ship_header column")
    return columns


def parse_location_metadata_tsv(
    contents: str,
    expected_seq_ids: List[str],
    existing_rows: Optional[List[Dict[str, Any]]] = None,
) -> tuple:
    """
    Parse TSV/CSV metadata and merge into location grid rows by ship header.

    Column headers must match exactly: ship_header (required), plus any of
    genus, species, strain, assembly_accession, hostchr, shipstart, shipend,
    strand. Only non-empty cells overwrite the table. Strand is inferred from
    start/end when omitted.

    Returns:
        (merged_rows, imported_row_count)
    """
    from io import StringIO
    from src.components.submission_location_grid import init_location_rows

    if not contents or not contents.strip():
        raise WebValidationError("Metadata file is empty")

    first_line = contents.strip().splitlines()[0]
    delimiter = "\t" if "\t" in first_line else ","

    try:
        df = pd.read_csv(
            StringIO(contents), sep=delimiter, dtype=str, keep_default_na=False
        )
    except Exception as exc:
        raise WebValidationError(f"Could not parse metadata file: {exc}") from exc

    if df.empty:
        raise WebValidationError("Metadata file has no data rows")

    df.columns = _validate_metadata_columns(list(df.columns))

    importable_fields = set(df.columns) - {"ship_header"}
    rows = init_location_rows(expected_seq_ids, existing_rows)
    rows_by_header = {row["ship_header"]: row for row in rows}

    imported_count = 0
    seen_headers = set()

    for _, file_row in df.iterrows():
        header = str(file_row["ship_header"]).strip()
        if not header:
            raise WebValidationError("Metadata file contains a blank ship header")

        if header in seen_headers:
            raise WebValidationError(
                f"Duplicate ship header in metadata file: '{header}'"
            )
        seen_headers.add(header)

        if header not in rows_by_header:
            raise WebValidationError(
                f"Unknown ship header in metadata file: '{header}'"
            )

        target = rows_by_header[header]
        imported_count += 1

        for field in importable_fields - {"strand"}:
            value = str(file_row.get(field, "")).strip()
            if value:
                target[field] = value

        explicit_strand = None
        if "strand" in importable_fields:
            explicit_strand = str(file_row.get("strand", "")).strip() or None

        target["strand"] = infer_metadata_strand(
            explicit_strand,
            target.get("shipstart"),
            target.get("shipend"),
        )

    if imported_count == 0:
        raise WebValidationError("Metadata file contains no usable data rows")

    return [rows_by_header[seq_id] for seq_id in expected_seq_ids], imported_count


def build_submission_entries(
    seq_contents: str,
    seq_filename: str,
    seq_date: Optional[float],
    anno_contents: Optional[str],
    anno_filename: Optional[str],
    anno_date: Optional[float],
    locations: List[Dict[str, Any]],
    uploader: str,
    evidence: str,
    comment: str = "",
    classification: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """
    Build validated submission payloads from one FASTA/GFF pair and location rows.

    Each location row corresponds to one FASTA header (same order as records).
    """
    records = parse_fasta_records(seq_contents, seq_filename)
    if len(locations) != len(records):
        raise WebValidationError(
            f"Expected location details for {len(records)} sequence(s), "
            f"but received {len(locations)}.",
            "location",
        )

    gff_seqids = get_gff_seqids(anno_contents) if anno_contents else set()
    entries = []

    for record, location in zip(records, locations):
        seq_id = record["id"]
        label = f"Sequence {seq_id}"
        validate_organism_fields(
            location.get("genus"),
            location.get("species"),
            label=label,
        )
        validate_location_fields(
            location.get("hostchr"),
            location.get("shipstart"),
            location.get("shipend"),
            label=label,
        )

        per_seq_gff = filter_gff_for_seqid(anno_contents, seq_id)
        if anno_contents and gff_seqids and not per_seq_gff:
            logger.info(f"No GFF rows matched FASTA header '{seq_id}'")

        entries.append(
            {
                "sequence": record["sequence"],
                "seq_id": seq_id,
                "filename": f"{seq_filename} ({seq_id})",
                "seq_contents": encode_single_fasta(seq_id, record["sequence"]),
                "seq_filename": seq_filename,
                "seq_date": seq_date,
                "anno_contents": per_seq_gff,
                "anno_filename": anno_filename if per_seq_gff else None,
                "anno_date": anno_date if per_seq_gff else None,
                "uploader": uploader,
                "evidence": evidence,
                "genus": location["genus"],
                "species": location["species"],
                "strain": location.get("strain"),
                "hostchr": location["hostchr"],
                "shipstart": int(location["shipstart"]),
                "shipend": int(location["shipend"]),
                "strand_radio": location.get("strand_radio", 1),
                "assembly_accession": location.get("assembly_accession"),
                "comment": comment or "",
                "classification": classification if len(entries) == 0 else None,
            }
        )

    return entries


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
        "strain": validated_data.get("strain"),
        "contig_id": validated_data["hostchr"],
        "element_start": start,
        "element_end": end,
        "element_strand": strand,
        "assembly_accession": validated_data.get("assembly_accession"),
        "curator": validated_data["uploader"],
        "curated_status": "needs_review",  # Default for web submissions
        "notes": validated_data.get("comment", ""),
        "duplicate_info": duplicate_info,
    }

    # Map BLAST classification from prefill to processor schema
    classification = validated_data.get("classification")
    if classification:
        processed_data["classification"] = classification
        if classification.get("family"):
            processed_data["ship_family"] = classification["family"]
        if classification.get("navis"):
            processed_data["ship_navis"] = classification["navis"]
        h = classification.get("haplotype")
        if h:
            processed_data["ship_haplotype"] = (
                h.get("haplotype_name") if isinstance(h, dict) else h
            )

    # Validate using submission_utils validator
    is_valid, validation_errors = validate_submission(processed_data)

    if not is_valid:
        raise WebValidationError("; ".join(validation_errors))

    return processed_data


def process_submission_to_staging(
    processed_data: Dict[str, Any],
    existing_ships: Optional[pd.DataFrame] = None,
) -> Dict[str, Any]:
    """
    Process submission for staging only (submissions DB). Does NOT insert into main DB.

    Validates, checks duplicates, assigns display accession. Used for new web submissions.
    Main DB integration (SubmissionProcessor) is for the future promotion flow.

    Args:
        processed_data: Output from process_submission_data()

    Returns:
        Dict with accession, needs_review, filename, uploader (same shape as perform_database_insertion)

    Raises:
        WebValidationError: If exact duplicate same taxon
    """
    duplicate_info = processed_data.get("duplicate_info")

    if duplicate_info and duplicate_info.is_duplicate:
        if not duplicate_info.different_taxon:
            raise WebValidationError(
                f"This sequence already exists in the database "
                f"(Accession: {duplicate_info.existing_accession}). "
                f"Duplicate submissions from the same organism are not allowed.",
                "sequence",
            )
        # Duplicate different taxon: use existing accession for display
        accession_tag = duplicate_info.existing_accession or "Pending"
        needs_review = True
        logger.info(
            f"Duplicate sequence from different taxon - staging only "
            f"(existing: {duplicate_info.existing_accession})"
        )
    else:
        # New sequence or BLAST exact match
        clean_seq = clean_sequence(processed_data["sequence"])
        classification = processed_data.get("classification")

        if (
            classification
            and classification.get("source") == "exact"
            and classification.get("closest_match")
        ):
            accession_tag = classification["closest_match"]
            needs_review = True
            logger.info(f"Using BLAST exact match for display: {accession_tag}")
        else:
            if existing_ships is None:
                existing_ships = build_accession_reference_pool()
            accession_tag, needs_review = assign_accession(clean_seq, existing_ships)

    return {
        "success": True,
        "ship_id": None,
        "accession": accession_tag,
        "needs_review": needs_review,
        "filename": processed_data.get("starshipID"),
        "uploader": processed_data.get("curator"),
        "warnings": [],
    }


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


def parse_submission_fasta(seq_contents: str) -> tuple:
    """Parse plain-text or base64 FASTA stored in submissions DB."""
    text = decode_upload_contents(seq_contents) if seq_contents else ""
    records = list(SeqIO.parse(StringIO(text), "fasta"))
    if not records:
        raise WebValidationError("No valid FASTA sequences in submission", "sequence")
    record = records[0]
    return record.id, str(record.seq)


def fetch_processed_submissions(
    exclude_submission_id: Optional[int] = None,
) -> pd.DataFrame:
    """Load processed staging submissions as an accession reference pool."""
    query = """
        SELECT id, seq_contents, accession_tag
        FROM submissions
        WHERE processing_status = 'processed'
          AND accession_tag IS NOT NULL
    """
    params = {}
    if exclude_submission_id is not None:
        query += " AND id != :exclude_id"
        params["exclude_id"] = exclude_submission_id

    rows = []
    with get_submissions_session() as session:
        result = session.execute(text(query), params)
        for row in result:
            try:
                _seq_id, sequence = parse_submission_fasta(row.seq_contents)
            except WebValidationError:
                continue
            clean_seq = clean_sequence(sequence)
            if not clean_seq or not row.accession_tag:
                continue
            rows.append(
                {
                    "accession_tag": row.accession_tag,
                    "accession_display": row.accession_tag,
                    "sequence": sequence,
                    "md5": generate_md5_hash(clean_seq),
                    "rev_comp_md5": generate_md5_hash(revcomp(clean_seq)),
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                "accession_tag",
                "accession_display",
                "sequence",
                "md5",
                "rev_comp_md5",
            ]
        )
    return pd.DataFrame(rows)


def merge_accession_reference_pools(
    main_ships: pd.DataFrame, staging_ships: pd.DataFrame
) -> pd.DataFrame:
    """Combine main-DB ships and processed staging submissions for assign_accession."""
    frames = [
        df for df in (main_ships, staging_ships) if df is not None and not df.empty
    ]
    if not frames:
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    if "accession_tag" in combined.columns:
        combined = combined.drop_duplicates(subset=["accession_tag"], keep="first")
    return combined


def build_accession_reference_pool(
    exclude_submission_id: Optional[int] = None,
) -> pd.DataFrame:
    """Build merged accession pool from main DB and processed staging rows."""
    main_ships = fetch_ships(with_sequence=True)
    staging_ships = fetch_processed_submissions(
        exclude_submission_id=exclude_submission_id
    )
    return merge_accession_reference_pools(main_ships, staging_ships)


def submission_row_to_validated_data(row: Dict[str, Any]) -> Dict[str, Any]:
    """Convert a submissions DB row into validated_data for process_submission_data."""
    seq_id, sequence = parse_submission_fasta(row["seq_contents"])
    strand = row.get("shipstrand") or "+"
    strand_radio = 1 if strand == "+" else 2

    classification = None
    if row.get("classification_family") or row.get("closest_match"):
        classification = {
            "family": row.get("classification_family"),
            "navis": row.get("classification_navis"),
            "haplotype": row.get("classification_haplotype"),
            "source": row.get("classification_source"),
            "closest_match": row.get("closest_match"),
            "confidence": row.get("classification_confidence"),
        }

    return {
        "sequence": sequence,
        "seq_id": seq_id,
        "filename": row.get("seq_filename") or seq_id,
        "uploader": row.get("uploader") or "",
        "evidence": row.get("evidence") or "",
        "genus": row.get("genus") or "",
        "species": row.get("species") or "",
        "strain": row.get("strain"),
        "hostchr": row.get("hostchr") or "",
        "shipstart": row.get("shipstart"),
        "shipend": row.get("shipend"),
        "strand_radio": strand_radio,
        "assembly_accession": row.get("assembly_accession"),
        "comment": row.get("comment") or "",
        "classification": classification,
    }


def update_staging_submission_after_process(
    sub_id: int, result: Dict[str, Any], processed_data: Dict[str, Any]
) -> None:
    """Persist accession and review flags after admin Process step."""
    classification = processed_data.get("classification") or {}
    with get_submissions_session() as session:
        session.execute(
            text(
                """
                UPDATE submissions
                SET accession_tag = :accession_tag,
                    needs_review = :needs_review,
                    processing_status = 'processed',
                    classification_source = COALESCE(:classification_source, classification_source),
                    classification_family = COALESCE(:classification_family, classification_family),
                    classification_navis = COALESCE(:classification_navis, classification_navis),
                    classification_haplotype = COALESCE(:classification_haplotype, classification_haplotype),
                    closest_match = COALESCE(:closest_match, closest_match),
                    classification_confidence = COALESCE(:classification_confidence, classification_confidence)
                WHERE id = :id
                """
            ),
            {
                "id": sub_id,
                "accession_tag": result.get("accession"),
                "needs_review": 1 if result.get("needs_review") else 0,
                "classification_source": classification.get("source"),
                "classification_family": processed_data.get("ship_family")
                or classification.get("family"),
                "classification_navis": processed_data.get("ship_navis")
                or classification.get("navis"),
                "classification_haplotype": processed_data.get("ship_haplotype")
                or classification.get("haplotype"),
                "closest_match": classification.get("closest_match"),
                "classification_confidence": classification.get("confidence"),
            },
        )
        session.commit()


def process_staging_submission(sub_id: int) -> Dict[str, Any]:
    """
    Run staging checks and accession assignment for one submission row.

    Used by the admin Process action (not at web submit time).
    """
    with get_submissions_session() as session:
        row = session.execute(
            text("SELECT * FROM submissions WHERE id = :id"), {"id": sub_id}
        ).fetchone()

    if not row:
        raise WebValidationError(f"Submission {sub_id} not found")

    row = dict(row._mapping)
    status = row.get("processing_status") or "pending"

    if status == "promoted":
        raise WebValidationError(f"Submission {sub_id} has already been promoted")
    if status == "processed" and row.get("accession_tag"):
        return {
            "success": True,
            "already_processed": True,
            "accession": row.get("accession_tag"),
            "needs_review": bool(row.get("needs_review")),
            "filename": row.get("seq_filename"),
            "uploader": row.get("uploader"),
        }

    validated = submission_row_to_validated_data(row)
    processed = process_submission_data(validated, validated["strand_radio"])
    reference_pool = build_accession_reference_pool(exclude_submission_id=sub_id)
    result = process_submission_to_staging(processed, existing_ships=reference_pool)
    update_staging_submission_after_process(sub_id, result, processed)
    result["success"] = True
    result["already_processed"] = False
    return result


def process_staging_submissions_ordered(sub_ids: List[int]) -> List[Dict[str, Any]]:
    """Process multiple submission rows sequentially (for grouped submits)."""
    results = []
    for sub_id in sorted(sub_ids):
        try:
            results.append({"sub_id": sub_id, **process_staging_submission(sub_id)})
        except WebValidationError as e:
            results.append(
                {"sub_id": sub_id, "success": False, "error": str(e.message)}
            )
        except Exception as e:
            logger.error(f"Failed to process submission {sub_id}: {e}", exc_info=True)
            results.append({"sub_id": sub_id, "success": False, "error": str(e)})
    return results


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

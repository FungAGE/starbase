import base64
import dash
from dash_iconify import DashIconify
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback, ctx, clientside_callback
from dash.exceptions import PreventUpdate

import datetime
from typing import Dict, Any, Optional
import re
from src.database.sql_engine import get_submissions_session
from src.database.models.schema import Submission
from sqlalchemy.exc import SQLAlchemyError
from src.components.ui import create_file_upload, _i, _lt

from src.config.logging import get_logger
from src.config.cache import cache
from src.utils.web_submission_adapter import (
    WebValidationError,
    parse_fasta_records,
    get_gff_seqids,
    build_submission_entries,
    decode_upload_contents,
    parse_location_metadata_tsv,
)
from src.components.submission_location_grid import (
    LOCATION_GRID_ID,
    create_location_grid,
    create_location_status,
    create_location_tsv_upload,
    init_location_rows,
    rows_to_locations,
)

import uuid
import json
from src.utils.email_notifications import (
    send_curator_notification,
    send_submission_confirmation,
)
from src.components.submission_queue import create_submission_queue_banner


logger = get_logger(__name__)


def validate_email(email: str) -> bool:
    """Validate email format."""
    if not email:
        return False
    pattern = r"^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$"
    return re.match(pattern, email) is not None


class SubmissionError(Exception):
    """Base exception for submission-related errors."""

    def __init__(
        self, message: str, error_type: str = "general", user_message: str = None
    ):
        self.message = message
        self.error_type = error_type
        self.user_message = user_message or message
        super().__init__(self.message)


class ValidationError(SubmissionError):
    """Exception for validation errors."""

    def __init__(self, message: str, field: str = None):
        super().__init__(message, "validation", message)
        self.field = field


class ProcessingError(SubmissionError):
    """Exception for processing errors."""

    def __init__(self, message: str, stage: str = None):
        super().__init__(message, "processing", message)
        self.stage = stage


class DatabaseError(SubmissionError):
    """Exception for database errors."""

    def __init__(self, message: str, operation: str = None):
        super().__init__(
            message, "database", "A database error occurred. Please try again."
        )
        self.operation = operation


def handle_submission_error(error: Exception) -> Dict[str, Any]:
    """
    Handle submission errors and return appropriate response.

    Args:
        error: The exception that occurred

    Returns:
        Dict containing modal state, message, and loading status
    """
    logger.error(f"Submission error: {str(error)}", exc_info=True)

    if isinstance(error, ValidationError):
        return {
            "modal_open": False,
            "message": dmc.Alert(
                title="Check your input",
                children=str(error.user_message),
                color="var(--mantine-color-orange-6)",
                variant="light",
            ),
            "loading": False,
        }
    elif isinstance(error, ProcessingError):
        return {
            "modal_open": False,
            "message": dmc.Alert(
                title="Processing failed",
                children=str(error.user_message),
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
            "loading": False,
        }
    elif isinstance(error, DatabaseError):
        return {
            "modal_open": False,
            "message": dmc.Alert(
                title="Something went wrong",
                children=dmc.Stack(
                    [
                        dmc.Text(str(error.user_message)),
                        dmc.Text(
                            "Try again in a moment. If the problem continues, contact support.",
                            size="sm",
                            c="dimmed",
                        ),
                    ],
                    gap="xs",
                ),
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
            "loading": False,
        }
    else:
        return {
            "modal_open": False,
            "message": dmc.Alert(
                title="Something went wrong",
                children=dmc.Stack(
                    [
                        dmc.Text(
                            "We couldn't complete your submission. Please try again."
                        ),
                        dmc.Text(
                            "If the problem continues, contact support.",
                            size="sm",
                            c="dimmed",
                        ),
                    ],
                    gap="xs",
                ),
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
            "loading": False,
        }


def create_submission_status(
    submission_id: str, initial_status: str = "queued"
) -> None:
    """
    Create a new submission status entry.

    Args:
        submission_id: Unique submission identifier
        initial_status: Initial status ("queued", "processing", etc.)
    """
    status_data = {
        "submission_id": submission_id,
        "status": initial_status,
        "created_at": datetime.datetime.now().isoformat(),
        "updated_at": datetime.datetime.now().isoformat(),
        "progress": 0,
        "message": "Submission queued for processing",
        "result": None,
    }
    cache.set(
        f"submission:{submission_id}", json.dumps(status_data), timeout=3600
    )  # 1 hour timeout


def update_submission_status(
    submission_id: str,
    status: str,
    progress: int = None,
    message: str = None,
    result: Dict = None,
) -> None:
    """
    Update submission status.

    Args:
        submission_id: Unique submission identifier
        status: New status
        progress: Progress percentage (0-100)
        message: Status message
        result: Final result data
    """
    cache_key = f"submission:{submission_id}"
    cached_data = cache.get(cache_key)

    if cached_data:
        try:
            status_data = json.loads(cached_data)
        except json.JSONDecodeError:
            status_data = {}
    else:
        status_data = {}

    # Update fields
    status_data.update(
        {"status": status, "updated_at": datetime.datetime.now().isoformat()}
    )

    if progress is not None:
        status_data["progress"] = progress

    if message is not None:
        status_data["message"] = message

    if result is not None:
        status_data["result"] = result

    cache.set(cache_key, json.dumps(status_data), timeout=3600)


def get_submission_status(submission_id: str) -> Optional[Dict]:
    """
    Get submission status.

    Args:
        submission_id: Unique submission identifier

    Returns:
        Status dict or None if not found
    """
    cached_data = cache.get(f"submission:{submission_id}")
    if cached_data:
        try:
            return json.loads(cached_data)
        except json.JSONDecodeError:
            return None
    return None


dash.register_page(__name__)

submission_header = dmc.Title(
    [
        "Submission of Candidate ",
        _i("Starship"),
        " Sequences to ",
        _lt("starbase"),
    ],
    order=1,
    mb="md",
)

submission_info_card = dmc.Paper(
    children=[
        dmc.Stack(
            [
                dmc.Text(
                    "Comparative genomics projects are a collaborative effort. Submit your Starship discoveries to the community to help us build the most comprehensive database of Starship elements.",
                    size="md",
                    mb="md",
                ),
                # TODO: "You'll receive a confirmation email once your submission is processed."
                dmc.Text(
                    "Each submission will be processed and assigned an accession number. It will also go through a manual review before being included in the next database release.",
                    size="md",
                    mb="md",
                ),
                dmc.Alert(
                    "Complete all fields marked with * before submitting.",
                    color="var(--mantine-color-red-6)",
                    variant="light",
                    title="Required fields",
                ),
            ],
            gap="md",
        ),
    ],
    p="xl",
    radius="md",
    withBorder=True,
    mb="xl",
    style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
)

submission_received_modal = dmc.Modal(
    id="submit-modal",
    opened=False,
    centered=True,
    size="lg",
    children=[
        dmc.Stack(
            [
                html.Div(id="output-data-upload"),
                dmc.Group(
                    dmc.Button(
                        "Got it",
                        id="close",
                        variant="light",
                        color="var(--mantine-color-indigo-6)",
                        n_clicks=0,
                    ),
                    justify="flex-end",
                ),
            ],
            gap="md",
        ),
    ],
)

# Section: Upload (side-by-side)
upload_section = dmc.Stack(
    [
        dmc.Title("Upload Files", order=2, mb="md"),
        dmc.Text(
            "Upload a single- or multi-sequence FASTA.",
            size="sm",
            c="dimmed",
            mb="sm",
        ),
        html.Div(id="submit-prefill-info"),
        dmc.Grid(
            [
                dmc.GridCol(
                    dmc.Paper(
                        children=[
                            dmc.Text(
                                [
                                    "Starship Sequence ",
                                    html.Span(
                                        "*",
                                        style={"color": "var(--mantine-color-red-6)"},
                                    ),
                                ],
                                fw=500,
                                mb="sm",
                            ),
                            create_file_upload(
                                upload_id="submit-fasta-upload",
                                output_id="submit-fasta-sequence-upload",
                                accept_types=[".fa", ".fas", ".fasta", ".fna"],
                                placeholder_text="Choose a FASTA file (.fa, .fasta, .fna)",
                            ),
                            dcc.Loading(
                                id="loading-1",
                                type="circle",
                                children=html.Div(id="loading-output-1"),
                            ),
                        ],
                        p="md",
                        radius="md",
                        withBorder=False,
                    ),
                    span={"base": 12, "md": 6},
                ),
                dmc.GridCol(
                    dmc.Paper(
                        children=[
                            dmc.Text("Gene Annotations (GFF3)", fw=500, mb="sm"),
                            create_file_upload(
                                upload_id="submit-upload-gff",
                                output_id="submit-output-gff-upload",
                                accept_types=[".gff", ".gff3", ".tsv"],
                                placeholder_text="Choose a GFF file (.gff, .gff3) — optional",
                            ),
                            dcc.Loading(
                                id="loading-2",
                                type="circle",
                                children=html.Div(id="loading-output-2"),
                            ),
                        ],
                        p="md",
                        radius="md",
                        withBorder=False,
                    ),
                    span={"base": 12, "md": 6},
                ),
            ],
        ),
    ],
    gap="md",
)

# Section: Contact
contact_section = dmc.Stack(
    [
        dmc.Title("Metadata", order=2, mb="md"),
        dmc.Grid(
            [
                dmc.GridCol(
                    dmc.Stack(
                        [
                            dmc.Text(
                                "Contact",
                                size="sm",
                                fw=600,
                                c="dimmed",
                                tt="uppercase",
                            ),
                            dmc.TextInput(
                                id="uploader",
                                label="Your email",
                                placeholder="e.g., you@example.com",
                                required=True,
                                leftSection=DashIconify(icon="fas fa-envelope"),
                            ),
                            dmc.Box(
                                id="email-validation-message",
                                style={"minHeight": "20px"},
                            ),
                            dmc.Select(
                                id="evidence",
                                label="Annotation tool or method",
                                placeholder="Select method",
                                required=True,
                                description="Tool or pipeline used to identify and annotate the Starship",
                                data=[
                                    {"value": "starfish", "label": "starfish"},
                                    {
                                        "value": "manual curation",
                                        "label": "manual curation",
                                    },
                                    {"value": "BLAST", "label": "BLAST"},
                                    {"value": "other", "label": "other"},
                                ],
                            ),
                        ],
                        gap="md",
                    ),
                    span={"base": 12, "md": 8},
                ),
            ],
        ),
    ],
    gap="md",
)

# Section: Location (editable grid — one row per FASTA header)
location_section = dmc.Stack(
    [
        dmc.Title("Location Details", order=2, mb="md"),
        html.Div(id="submit-location-status", children=create_location_status()),
        create_location_tsv_upload(),
        html.Div(id="submit-location-tsv-feedback"),
        create_location_grid(),
    ],
    gap="md",
)

# Section: Notes (progressive disclosure - collapsed by default)
notes_section = dmc.Accordion(
    variant="contained",
    chevronPosition="right",
    children=[
        dmc.AccordionItem(
            value="notes",
            children=[
                dmc.AccordionControl(
                    dmc.Text(
                        "Add optional notes, i.e. linked publications, etc.",
                        size="sm",
                        c="dimmed",
                    ),
                ),
                dmc.AccordionPanel(
                    dmc.Textarea(
                        id="comment",
                        label="Notes",
                        placeholder="e.g., unusual features, annotation notes, or host genome context",
                        description="Any details that would help curators evaluate this submission",
                        minRows=3,
                        autosize=True,
                        maxRows=6,
                    ),
                ),
            ],
        ),
    ],
)

# Sticky submit button
submit_button_section = dmc.Box(
    dmc.Stack(
        [
            dmc.Checkbox(
                id="submit-consent-checkbox",
                label="I understand that accepted submissions will be publicly available in the Starbase database",
                checked=False,
                color="indigo",
            ),
            dmc.Center(
                dmc.Button(
                    "Submit Starship(s)",
                    id="submit-ship",
                    size="lg",
                    variant="filled",
                    color="indigo",
                    loading=False,
                    disabled=True,
                    fullWidth=False,
                ),
            ),
        ],
        gap="md",
        align="center",
    ),
    mt="xl",
    pt="md",
)

submission_body = dmc.Paper(
    children=dmc.Stack(
        [
            upload_section,
            dmc.Divider(),
            contact_section,
            dmc.Divider(),
            location_section,
            dmc.Divider(),
            notes_section,
            submit_button_section,
        ],
        gap="xl",
    ),
    p="xl",
    radius="md",
    withBorder=True,
    style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
)

data_policy_card = dmc.Alert(
    dmc.Stack(
        [
            dmc.Text(
                "This database is an academic resource. Submitted data is used solely to:",
                size="sm",
                fw=500,
            ),
            dmc.List(
                [
                    dmc.ListItem(
                        dmc.Text(
                            "Build and maintain a public database for scientific research",
                            size="sm",
                        )
                    ),
                    dmc.ListItem(
                        dmc.Text(
                            "Attribute your contribution — your name or identifier will be associated with the database entry, not your email address",
                            size="sm",
                        )
                    ),
                    dmc.ListItem(
                        dmc.Text(
                            "Contact you about your submission (confirmation and curation questions)",
                            size="sm",
                        )
                    ),
                ],
                size="sm",
                spacing="xs",
            ),
            dmc.Divider(),
            dmc.Text(
                "Accepted submissions are included in public database releases and will be publicly accessible.",
                size="sm",
                fw=500,
            ),
            dmc.Text(
                "Only submit sequences you are comfortable making publicly available. If your data is from unpublished work, consider whether public release is appropriate at this time.",
                size="sm",
            ),
        ],
        gap="xs",
    ),
    title="How this data is used",
    color="var(--mantine-color-indigo-6)",
    variant="light",
    mb="md",
)

top_section = dmc.Stack(
    [data_policy_card, create_submission_queue_banner(max_items=5)],
    gap="md",
    mb="xl",
)

layout = dmc.Container(
    size="lg",
    children=[
        dcc.Location(id="submit-url", refresh=False),
        dcc.Store(id="submit-fasta-prefill"),
        dcc.Store(id="submit-primary-seq-ids", data=[]),
        dcc.Store(id="submit-location-rows", data=[]),
        submission_header,
        top_section,
        submission_info_card,
        submission_body,
        submission_received_modal,
    ],
    style={
        "margin": "0 auto",
        "padding": "var(--mantine-spacing-xl)",
    },
)


@callback(
    [
        Output("submit-fasta-prefill", "data"),
        Output("comment", "value", allow_duplicate=True),
        Output("evidence", "value", allow_duplicate=True),
        Output("submit-prefill-info", "children"),
    ],
    Input("submit-url", "search"),
    prevent_initial_call=True,
)
def load_blast_prefill(search):
    """Read blast prefill data from cache and populate the form."""
    from urllib.parse import parse_qs
    import json

    if not search:
        raise dash.exceptions.PreventUpdate

    params = parse_qs(search.lstrip("?"))
    blast_id = params.get("blast_id", [None])[0]
    if not blast_id:
        raise dash.exceptions.PreventUpdate

    cached = cache.get(f"submit_prefill:{blast_id}")
    if not cached:
        raise dash.exceptions.PreventUpdate

    prefill_data = json.loads(cached)
    fasta_contents = prefill_data.get("fasta_contents")
    fasta_filename = prefill_data.get("fasta_filename", "sequence_from_blast.fa")
    comment = prefill_data.get("comment", "")
    classification = prefill_data.get("classification")

    if not fasta_contents:
        info_banner = dmc.Alert(
            "We pre-filled metadata from your BLAST search, but the sequence file wasn't found. Upload your FASTA file below.",
            title="Partial pre-fill from BLAST",
            color="var(--mantine-color-orange-6)",
            variant="light",
            mb="sm",
        )
        partial_prefill = {"classification": classification} if classification else None
        return partial_prefill, comment, "BLAST", info_banner

    try:
        records = parse_fasta_records(fasta_contents, fasta_filename)
        n_seqs = len(records)
    except Exception as e:
        logger.error(f"Error decoding prefilled FASTA: {e}")
        n_seqs = None

    info_banner = dmc.Alert(
        dmc.Stack(
            [
                dmc.Text(
                    f"Sequence from BLAST: {fasta_filename}"
                    + (
                        f" ({n_seqs} sequence{'s' if n_seqs != 1 else ''})"
                        if n_seqs
                        else ""
                    ),
                    fw=500,
                ),
                dmc.Text(
                    "Upload a different file below to replace it.",
                    size="sm",
                    c="dimmed",
                ),
            ],
            gap="xs",
        ),
        title="Pre-filled from BLAST",
        color="var(--mantine-color-green-6)",
        variant="light",
        mb="sm",
    )

    prefill_store = {"contents": fasta_contents, "filename": fasta_filename}
    if classification:
        prefill_store["classification"] = classification

    return (
        prefill_store,
        comment,
        "BLAST",
        info_banner,
    )


def _enrich_classification_from_db(classification: dict) -> dict:
    """
    When classification has closest_match but missing family/navis/haplotype,
    look them up from the main DB and merge into the classification dict.
    """
    if not classification or not classification.get("closest_match"):
        return classification

    if (
        classification.get("family")
        and classification.get("navis")
        and classification.get("haplotype")
    ):
        return classification

    try:
        from src.utils.submission_utils import lookup_classification_from_accession

        looked_up = lookup_classification_from_accession(
            classification["closest_match"]
        )
        if looked_up:
            if not classification.get("family") and looked_up.get("family"):
                classification = {**classification, "family": looked_up["family"]}
            if not classification.get("navis") and looked_up.get("navis"):
                classification = {**classification, "navis": looked_up["navis"]}
            if not classification.get("haplotype") and looked_up.get("haplotype"):
                classification = {
                    **classification,
                    "haplotype": looked_up["haplotype"],
                }
    except Exception as e:
        logger.warning(f"Could not lookup classification from DB: {e}")

    return classification


# Function to insert a new submission into the database
def insert_submission(
    seq_contents,
    seq_filename,
    seq_date,
    anno_contents,
    anno_filename,
    anno_date,
    uploader,
    evidence,
    genus,
    species,
    hostchr,
    shipstart,
    shipend,
    shipstrand,
    comment,
    accession=None,
    needs_review=False,
    classification=None,
    assembly_accession=None,
    strain=None,
    submission_group_id=None,
    processing_status="pending",
):
    """
    Insert a new submission record into the submissions table (queue).

    Returns:
        int: The database ID of the inserted submission, or None on failure.
    """
    try:
        # Accept base64 seq_contents (data:xxx;base64,...) or raw sequence string
        if "," in str(seq_contents) and "base64" in str(seq_contents):
            _ct, content_string = seq_contents.split(",", 1)
            seq_decoded = base64.b64decode(content_string).decode("utf-8")
        else:
            seq_decoded = seq_contents if isinstance(seq_contents, str) else ""

        seq_datetime_str = (
            datetime.datetime.fromtimestamp(seq_date).strftime("%Y-%m-%d %H:%M:%S")
            if seq_date
            else datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )

        anno_decoded = None
        anno_filename_val = None
        anno_datetime_str = None
        if anno_contents and anno_date:
            if "," in str(anno_contents) and "base64" in str(anno_contents):
                _ct, content_string = anno_contents.split(",", 1)
                anno_decoded = base64.b64decode(content_string).decode("utf-8")
            anno_filename_val = anno_filename
            anno_datetime_str = datetime.datetime.fromtimestamp(anno_date).strftime(
                "%Y-%m-%d %H:%M:%S"
            )

        with get_submissions_session() as session:
            submission = Submission(
                seq_contents=seq_decoded,
                seq_filename=seq_filename,
                seq_date=seq_datetime_str,
                anno_contents=anno_decoded,
                anno_filename=anno_filename_val,
                anno_date=anno_datetime_str,
                uploader=uploader,
                evidence=evidence,
                genus=genus,
                species=species,
                strain=strain,
                hostchr=hostchr,
                shipstart=shipstart,
                shipend=shipend,
                shipstrand=shipstrand,
                assembly_accession=assembly_accession,
                comment=comment,
                accession_tag=accession,
                needs_review=needs_review,
                submission_group_id=submission_group_id,
                processing_status=processing_status,
                classification_source=classification.get("source")
                if classification
                else None,
                classification_family=classification.get("family")
                if classification
                else None,
                classification_navis=classification.get("navis")
                if classification
                else None,
                classification_haplotype=(
                    classification["haplotype"].get("haplotype_name")
                    if classification
                    and isinstance(classification.get("haplotype"), dict)
                    else (classification.get("haplotype") if classification else None)
                ),
                closest_match=classification.get("closest_match")
                if classification
                else None,
                classification_confidence=classification.get("confidence")
                if classification
                else None,
            )
            session.add(submission)
            session.commit()
            session.refresh(submission)
            logger.debug(
                f"Successfully inserted submission {submission.id} for {seq_filename}"
            )
            return submission.id

    except SQLAlchemyError as e:
        logger.error(f"Database error during submission: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error processing submission: {str(e)}")
        raise


def update_submission_record(
    db_submission_id: int, accession: str = None, needs_review: bool = None
):
    """
    Update a submission record with completion data.

    Args:
        db_submission_id: Database ID of the submission
        accession: Assigned accession tag
        needs_review: Whether the submission needs curator review
    """
    try:
        from sqlalchemy import update
        from src.database.models.schema import Submission

        updates = {}
        if accession is not None:
            updates["accession_tag"] = accession
        if needs_review is not None:
            updates["needs_review"] = needs_review
        if not updates:
            return

        with get_submissions_session() as session:
            stmt = (
                update(Submission)
                .where(Submission.id == db_submission_id)
                .values(**updates)
            )
            session.execute(stmt)
            session.commit()
            logger.debug(f"Updated submission {db_submission_id}: {updates}")
    except Exception as e:
        logger.error(f"Failed to update submission record {db_submission_id}: {str(e)}")


def create_fasta_display(records, filename):
    """Create HTML components to display FASTA file info below the upload box."""
    seq_ids = [record["id"] for record in records]
    children = [
        dmc.Text(f"✓ File: {filename}", size="sm", fw=500),
        dmc.Text(f"Sequences found: {len(records)}", size="sm", c="dimmed"),
    ]
    if len(records) > 1:
        children.append(
            dmc.Text(
                f"Headers: {', '.join(seq_ids[:5])}"
                + ("…" if len(seq_ids) > 5 else ""),
                size="xs",
                c="dimmed",
            )
        )
    children.append(
        dmc.Text(
            "Fill in organism and location details for each sequence in the table below.",
            size="xs",
            c="dimmed",
        )
    )
    return dmc.Alert(
        children=children,
        title="FASTA loaded",
        color="var(--mantine-color-green-6)",
        variant="light",
    )


def upload_box_status(message, *, color="dimmed", fw=500):
    """Compact one-line status shown inside the upload drop zone."""
    return dmc.Text(message, size="sm", c=color, fw=fw, ta="center")


def create_gff_display(anno_filename, row_count, fasta_seq_ids=None, gff_seqids=None):
    """Create HTML components to display GFF file info."""
    children = [
        dmc.Text(f"✓ File: {anno_filename}", size="sm", fw=500),
        dmc.Text(f"Annotation rows: {row_count}", size="sm", c="dimmed"),
    ]
    if fasta_seq_ids and gff_seqids is not None:
        matched = [seq_id for seq_id in fasta_seq_ids if seq_id in gff_seqids]
        unmatched = [seq_id for seq_id in fasta_seq_ids if seq_id not in gff_seqids]
        if matched:
            children.append(
                dmc.Text(
                    f"Matched FASTA headers: {', '.join(matched[:5])}"
                    + ("…" if len(matched) > 5 else ""),
                    size="xs",
                    c="var(--mantine-color-green-7)",
                )
            )
        if unmatched:
            children.append(
                dmc.Text(
                    f"No GFF rows for: {', '.join(unmatched[:5])}"
                    + ("…" if len(unmatched) > 5 else ""),
                    size="xs",
                    c="var(--mantine-color-orange-7)",
                )
            )
    return dmc.Alert(
        children=children,
        title="GFF loaded",
        color="var(--mantine-color-green-6)",
        variant="light",
    )


def _resolve_primary_fasta(fasta_prefill, seq_contents, seq_filename, seq_date):
    """Return primary FASTA upload values, preferring manual upload over BLAST prefill."""
    if not seq_contents and fasta_prefill:
        seq_contents = fasta_prefill.get("contents")
        seq_filename = seq_filename or fasta_prefill.get(
            "filename", "sequence_from_blast.fa"
        )
        seq_date = seq_date or datetime.datetime.now().timestamp()
    return seq_contents, seq_filename, seq_date


def collect_all_submission_entries(
    fasta_prefill,
    seq_contents,
    seq_filename,
    seq_date,
    anno_contents,
    anno_filename,
    anno_date,
    location_rows,
    uploader,
    evidence,
    comment,
):
    """Collect and validate all Starship entries from the location grid."""
    errors = []
    if not uploader:
        errors.append("Curator email is required")
    if not evidence:
        errors.append("Evidence/annotation method is required")
    if errors:
        raise WebValidationError("; ".join(errors))

    classification = (
        fasta_prefill.get("classification")
        if fasta_prefill and fasta_prefill.get("classification")
        else None
    )

    seq_contents, seq_filename, seq_date = _resolve_primary_fasta(
        fasta_prefill, seq_contents, seq_filename, seq_date
    )
    if not seq_contents:
        raise WebValidationError("FASTA file is required")

    records = parse_fasta_records(seq_contents, seq_filename)
    expected_count = len(records)

    if not location_rows or len(location_rows) != expected_count:
        raise WebValidationError(
            "Provide organism and location details for each sequence in the table below.",
            "location",
        )

    locations = rows_to_locations(location_rows)

    return build_submission_entries(
        seq_contents=seq_contents,
        seq_filename=seq_filename,
        seq_date=seq_date,
        anno_contents=anno_contents,
        anno_filename=anno_filename,
        anno_date=anno_date,
        locations=locations,
        uploader=uploader,
        evidence=evidence,
        comment=comment or "",
        classification=classification,
    )


def insert_submission_from_entry(entry: dict, submission_group_id: str) -> int:
    """Insert one validated entry into the submissions queue (no processing)."""
    shipstrand = "+" if (entry.get("strand_radio") or 1) == 1 else "-"
    seq_date = entry.get("seq_date") or datetime.datetime.now().timestamp()
    return insert_submission(
        seq_contents=entry.get("seq_contents"),
        seq_filename=entry.get("seq_filename") or entry.get("filename"),
        seq_date=seq_date,
        anno_contents=entry.get("anno_contents"),
        anno_filename=entry.get("anno_filename"),
        anno_date=entry.get("anno_date"),
        uploader=entry["uploader"],
        evidence=entry["evidence"],
        genus=entry["genus"],
        species=entry["species"],
        strain=entry.get("strain"),
        hostchr=entry["hostchr"],
        shipstart=entry["shipstart"],
        shipend=entry["shipend"],
        shipstrand=shipstrand,
        comment=entry.get("comment") or "",
        accession=None,
        needs_review=True,
        classification=entry.get("classification"),
        assembly_accession=entry.get("assembly_accession"),
        submission_group_id=submission_group_id,
        processing_status="pending",
    )


def build_submission_success_message(group_id: str, ship_count: int, entries: list):
    """Build modal feedback for a grouped submission upload."""
    ship_labels = [
        entry.get("seq_id") or entry.get("filename") or f"Ship {index + 1}"
        for index, entry in enumerate(entries)
    ]
    list_items = [dmc.ListItem(label) for label in ship_labels]

    return dmc.Alert(
        children=[
            dmc.Text(
                f"Uploaded {ship_count} Starship{'s' if ship_count != 1 else ''}.",
                fw=600,
                size="lg",
            ),
            dmc.Text(f"Submission ID: {group_id}", size="sm", c="dimmed"),
            html.Br(),
            dmc.List(list_items, size="sm", spacing="xs")
            if len(list_items) > 1
            else None,
            dmc.Text(
                "Your submission has been received and will be reviewed by our curation team. "
                "You will be contacted by email if we need more information.",
                size="sm",
            ),
        ],
        title="Submission received",
        color="var(--mantine-color-green-6)",
        variant="light",
    )


@callback(
    Output("submit-ship", "disabled"),
    [
        Input("submit-consent-checkbox", "checked"),
        Input("submit-ship", "loading"),
    ],
)
def toggle_submit_button(checked, loading):
    if loading:
        return True
    return not bool(checked)


@callback(
    Output("submit-consent-checkbox", "disabled"),
    Input("submit-ship", "loading"),
)
def disable_consent_during_submit(loading):
    return bool(loading)


@callback(
    [
        Output("submit-modal", "opened"),
        Output("output-data-upload", "children"),
        Output("uploader", "value"),
        Output("evidence", "value", allow_duplicate=True),
        Output("comment", "value", allow_duplicate=True),
        Output("submit-primary-seq-ids", "data", allow_duplicate=True),
        Output("submit-location-rows", "data", allow_duplicate=True),
        Output(LOCATION_GRID_ID, "rowData", allow_duplicate=True),
        Output("submit-location-status", "children", allow_duplicate=True),
    ],
    [
        Input("submit-ship", "n_clicks"),
        Input("close", "n_clicks"),
    ],
    [
        State("submit-modal", "opened"),
        State("submit-fasta-prefill", "data"),
        State("submit-fasta-upload", "contents"),
        State("submit-fasta-upload", "filename"),
        State("submit-fasta-upload", "last_modified"),
        State("submit-upload-gff", "contents"),
        State("submit-upload-gff", "filename"),
        State("submit-upload-gff", "last_modified"),
        State("submit-location-rows", "data"),
        State("uploader", "value"),
        State("evidence", "value"),
        State("comment", "value"),
    ],
    running=[
        (Output("submit-ship", "loading"), True, False),
    ],
    prevent_initial_call=True,
    id="submit-ship-callback",
)
def submit_ship(
    submit_clicks,
    close_clicks,
    modal_opened,
    fasta_prefill,
    seq_contents,
    seq_filename,
    seq_date,
    anno_contents,
    anno_filename,
    anno_date,
    location_rows,
    uploader,
    evidence,
    comment,
):
    """Validate and queue one or more Starship submissions."""
    triggered_id = ctx.triggered_id if ctx.triggered else None
    no_update_tail = (dash.no_update,) * 7

    if triggered_id == "close":
        return (False, "", *no_update_tail)

    if triggered_id != "submit-ship" or not submit_clicks:
        raise PreventUpdate

    if uploader and not validate_email(uploader):
        return (
            True,
            dmc.Alert(
                "Use a valid email address (e.g., name@example.com) so we can send your confirmation.",
                title="Invalid email",
                color="var(--mantine-color-orange-6)",
                variant="light",
            ),
            *no_update_tail,
        )

    try:
        entries = collect_all_submission_entries(
            fasta_prefill=fasta_prefill,
            seq_contents=seq_contents,
            seq_filename=seq_filename,
            seq_date=seq_date,
            anno_contents=anno_contents,
            anno_filename=anno_filename,
            anno_date=anno_date,
            location_rows=location_rows,
            uploader=uploader,
            evidence=evidence,
            comment=comment,
        )

        group_id = str(uuid.uuid4())
        for entry in entries:
            insert_submission_from_entry(entry, group_id)

        try:
            send_curator_notification(
                group_id,
                entries,
                uploader=uploader,
                evidence=evidence,
                comment=comment,
            )
            if uploader:
                send_submission_confirmation(
                    uploader, group_id, ship_count=len(entries)
                )
        except Exception as e:
            logger.error(f"Failed to send email notifications: {str(e)}")

        message = build_submission_success_message(group_id, len(entries), entries)

        return (
            True,
            message,
            "",
            "",
            "",
            [],
            [],
            [],
            create_location_status(),
        )

    except WebValidationError as e:
        validation_error = ValidationError(
            str(e.message) if hasattr(e, "message") else str(e),
            getattr(e, "field", None),
        )
        error_response = handle_submission_error(validation_error)
        return (
            error_response["modal_open"],
            error_response["message"],
            *no_update_tail,
        )

    except SubmissionError as e:
        error_response = handle_submission_error(e)
        return (
            error_response["modal_open"],
            error_response["message"],
            *no_update_tail,
        )

    except Exception as e:
        logger.error(f"Unexpected error in submit_ship: {str(e)}", exc_info=True)
        error_response = handle_submission_error(e)
        return (
            error_response["modal_open"],
            error_response["message"],
            *no_update_tail,
        )


@callback(
    [
        Output("submit-location-rows", "data"),
        Output(LOCATION_GRID_ID, "rowData"),
        Output("submit-primary-seq-ids", "data"),
        Output("submit-location-status", "children"),
        Output(LOCATION_GRID_ID, "style"),
        Output(LOCATION_GRID_ID, "dashGridOptions"),
    ],
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-prefill", "data"),
    ],
    State("submit-location-rows", "data"),
)
def update_primary_location_fields(
    seq_contents, seq_filename, fasta_prefill, existing_rows
):
    """Initialize location grid rows from FASTA headers."""
    from src.components.submission_location_grid import grid_height

    resolved_contents, resolved_filename, _seq_date = _resolve_primary_fasta(
        fasta_prefill, seq_contents, seq_filename, None
    )
    if not resolved_contents:
        return (
            [],
            [],
            [],
            create_location_status(),
            {"width": "100%", "height": "200px"},
            dash.no_update,
        )

    try:
        records = parse_fasta_records(
            resolved_contents, resolved_filename or "upload.fa"
        )
    except WebValidationError:
        return (
            [],
            [],
            [],
            create_location_status(),
            {"width": "100%", "height": "200px"},
            dash.no_update,
        )

    seq_ids = [record["id"] for record in records]
    rows = init_location_rows(seq_ids, existing_rows)
    grid_style = {"width": "100%", "height": grid_height(len(rows))}
    grid_options = {
        "suppressPropertyNamesCheck": True,
        "rowHeight": 40,
        "headerHeight": 44,
        "enableCellTextSelection": True,
        "ensureDomOrder": True,
    }
    return (
        rows,
        rows,
        seq_ids,
        create_location_status(seq_ids),
        grid_style,
        grid_options,
    )


clientside_callback(
    """
    function(ev, rowData) {
        if (!ev || !rowData) {
            return [window.dash_clientside.no_update, window.dash_clientside.no_update];
        }
        var shipHeader = ev.data && ev.data.ship_header;
        var colId = ev.colId;
        if (!shipHeader || !colId) {
            return [window.dash_clientside.no_update, window.dash_clientside.no_update];
        }
        var newValue = (ev.data && ev.data[colId] !== undefined) ? ev.data[colId] : ev.value;
        var updatedRows = rowData.map(function(row) {
            if (row.ship_header !== shipHeader) {
                return row;
            }
            var copy = Object.assign({}, row);
            copy[colId] = newValue;
            return copy;
        });
        return [updatedRows, updatedRows];
    }
    """,
    Output("submit-location-rows", "data", allow_duplicate=True),
    Output(LOCATION_GRID_ID, "rowData", allow_duplicate=True),
    Input(LOCATION_GRID_ID, "cellValueChanged"),
    State(LOCATION_GRID_ID, "rowData"),
    prevent_initial_call=True,
)


@callback(
    [
        Output("submit-location-rows", "data", allow_duplicate=True),
        Output(LOCATION_GRID_ID, "rowData", allow_duplicate=True),
        Output("submit-location-tsv-feedback", "children"),
    ],
    Input("submit-location-tsv-upload", "contents"),
    [
        State("submit-location-tsv-upload", "filename"),
        State("submit-primary-seq-ids", "data"),
        State("submit-location-rows", "data"),
    ],
    prevent_initial_call=True,
)
def import_location_tsv(contents, filename, seq_ids, existing_rows):
    """Import location metadata from TSV/CSV matched by ship header."""
    if not contents:
        raise PreventUpdate

    if not seq_ids:
        return (
            dash.no_update,
            dash.no_update,
            dmc.Alert(
                "Upload a FASTA file before importing metadata.",
                title="No sequences loaded",
                color="var(--mantine-color-orange-6)",
                variant="light",
            ),
        )

    try:
        decoded = decode_upload_contents(contents)
        rows, imported_count = parse_location_metadata_tsv(
            decoded, seq_ids, existing_rows
        )
        return (
            rows,
            rows,
            dmc.Alert(
                f"Imported metadata for {imported_count} row(s) from {filename or 'uploaded file'}. "
                "Empty cells were skipped; ships not listed in the file were unchanged.",
                title="Metadata imported",
                color="var(--mantine-color-green-6)",
                variant="light",
            ),
        )
    except WebValidationError as exc:
        return (
            dash.no_update,
            dash.no_update,
            dmc.Alert(
                str(exc.message),
                title="Metadata import failed",
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
        )
    except Exception as exc:
        logger.error(f"Failed to import location metadata: {exc}", exc_info=True)
        return (
            dash.no_update,
            dash.no_update,
            dmc.Alert(
                "Could not read the metadata file. Check the format and try again.",
                title="Metadata import failed",
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
        )


@callback(
    Output("email-validation-message", "children"),
    Input("uploader", "value"),
    prevent_initial_call=True,
)
def validate_email_input(email):
    """Provide real-time email validation feedback."""
    if not email or email.strip() == "":
        return ""

    if validate_email(email):
        return dmc.Text("✓ Valid email", size="xs", c="var(--mantine-color-green-6)")
    else:
        return dmc.Text(
            "Enter a valid email (e.g., name@example.com)",
            size="xs",
            c="var(--mantine-color-red-6)",
        )


@callback(
    [
        Output("submit-fasta-sequence-upload", "children"),
        Output("loading-output-1", "children"),
    ],
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-prefill", "data"),
    ],
)
def update_fasta_details(seq_contents, seq_filename, fasta_prefill):
    placeholder = "Choose a FASTA file (.fa, .fasta, .fna)"

    # Manually uploaded file always takes priority
    if seq_contents is not None:
        try:
            filename = seq_filename or "upload.fa"
            records = parse_fasta_records(seq_contents, filename)
            return (
                upload_box_status(
                    "✓ File selected", color="var(--mantine-color-green-7)"
                ),
                create_fasta_display(records, filename),
            )
        except WebValidationError as e:
            return (
                upload_box_status("Invalid file", color="var(--mantine-color-red-6)"),
                dmc.Alert(
                    str(e.message),
                    title="File error",
                    color="var(--mantine-color-red-6)",
                    variant="light",
                ),
            )
        except Exception as e:
            logger.error(e)
            return (
                upload_box_status("Invalid file", color="var(--mantine-color-red-6)"),
                dmc.Alert(
                    "We couldn't parse this FASTA file. Check that it's valid and try again.",
                    title="File error",
                    color="var(--mantine-color-red-6)",
                    variant="light",
                ),
            )

    # Show prefill details when no file has been manually uploaded
    if fasta_prefill:
        prefill_contents = fasta_prefill.get("contents")
        prefill_filename = fasta_prefill.get("filename", "sequence_from_blast.fa")
        if prefill_contents:
            try:
                records = parse_fasta_records(prefill_contents, prefill_filename)
                return (
                    upload_box_status(
                        "✓ Pre-filled from BLAST",
                        color="var(--mantine-color-green-7)",
                    ),
                    dmc.Alert(
                        dmc.Stack(
                            [
                                dmc.Text(
                                    f"File name: {prefill_filename}",
                                    size="sm",
                                    fw=500,
                                ),
                                dmc.Text(
                                    f"Number of sequences: {len(records)}",
                                    size="sm",
                                ),
                                dmc.Text(
                                    "Upload a different file above to replace.",
                                    size="xs",
                                    c="dimmed",
                                ),
                            ],
                            gap="xs",
                        ),
                        title="FASTA pre-filled",
                        color="var(--mantine-color-green-6)",
                        variant="light",
                    ),
                )
            except Exception as e:
                logger.error(e)

    return placeholder, html.Div()


@callback(
    [
        Output("submit-output-gff-upload", "children"),
        Output("loading-output-2", "children"),
    ],
    [
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-prefill", "data"),
    ],
)
def update_gff_details(
    anno_contents, anno_filename, seq_contents, seq_filename, fasta_prefill
):
    placeholder = "Choose a GFF file (.gff, .gff3) — optional"

    if anno_contents is None:
        return placeholder, html.Div()
    try:
        decoded = decode_upload_contents(anno_contents)
        row_count = sum(
            1
            for line in decoded.split("\n")
            if line.strip() and not line.startswith("#") and "\t" in line
        )
        fasta_seq_ids = []
        resolved_contents, resolved_filename, _seq_date = _resolve_primary_fasta(
            fasta_prefill, seq_contents, seq_filename, None
        )
        if resolved_contents:
            try:
                fasta_seq_ids = [
                    record["id"]
                    for record in parse_fasta_records(
                        resolved_contents, resolved_filename or "upload.fa"
                    )
                ]
            except WebValidationError:
                fasta_seq_ids = []

        return (
            upload_box_status("✓ File selected", color="var(--mantine-color-green-7)"),
            create_gff_display(
                anno_filename,
                row_count,
                fasta_seq_ids=fasta_seq_ids,
                gff_seqids=get_gff_seqids(anno_contents),
            ),
        )
    except Exception as e:
        logger.error(e)
        return (
            upload_box_status("Invalid file", color="var(--mantine-color-red-6)"),
            dmc.Alert(
                "We couldn't parse this GFF file. Ensure it's valid GFF3 format and try again.",
                title="File error",
                color="var(--mantine-color-red-6)",
                variant="light",
            ),
        )

import base64
import dash
from dash_iconify import DashIconify
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State, ALL
from dash import dcc, html, callback, ctx
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
from src.config.celery_config import run_task
from src.tasks import process_submission_task
from src.utils.web_submission_adapter import (
    WebValidationError,
    parse_fasta_records,
    get_gff_seqids,
    build_submission_entries,
    decode_upload_contents,
)

import uuid
import json
from src.utils.email_notifications import (
    send_curator_notification,
    send_submission_confirmation,
)
from src.components.submission_queue import create_submission_queue


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


def create_location_placeholder():
    """Placeholder before a FASTA file is uploaded."""
    return dmc.Text(
        "Upload a FASTA file above to enter organism and location details for each sequence.",
        size="sm",
        c="dimmed",
    )


def create_location_table(seq_ids):
    """Table with one row per FASTA header for host genome location metadata."""
    rows = []
    for index, seq_id in enumerate(seq_ids):
        rows.append(
            html.Tr(
                [
                    html.Td(
                        dmc.Text(
                            seq_id,
                            size="sm",
                            fw=500,
                            style={"wordBreak": "break-word", "maxWidth": "180px"},
                        )
                    ),
                    html.Td(
                        dmc.TextInput(
                            id={"type": "primary-loc-genus", "index": index},
                            placeholder="e.g., Alternaria",
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.TextInput(
                            id={"type": "primary-loc-species", "index": index},
                            placeholder="e.g., alternata",
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.TextInput(
                            id={"type": "primary-loc-strain", "index": index},
                            placeholder="optional",
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.TextInput(
                            id={"type": "primary-loc-assembly", "index": index},
                            placeholder="e.g., GCA_000001305.1",
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.TextInput(
                            id={"type": "primary-loc-hostchr", "index": index},
                            placeholder="e.g., chr1",
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.NumberInput(
                            id={"type": "primary-loc-shipstart", "index": index},
                            placeholder="Start",
                            min=1,
                            step=1,
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.NumberInput(
                            id={"type": "primary-loc-shipend", "index": index},
                            placeholder="End",
                            min=1,
                            step=1,
                            size="sm",
                        )
                    ),
                    html.Td(
                        dmc.SegmentedControl(
                            id={"type": "primary-loc-strand", "index": index},
                            data=[
                                {"label": "+", "value": "1"},
                                {"label": "-", "value": "2"},
                            ],
                            value="1",
                            size="xs",
                        )
                    ),
                ]
            )
        )

    children = [
        dmc.Text(
            "Starship submissions",
            size="sm",
            fw=600,
            c="dimmed",
            tt="uppercase",
        ),
        dmc.Text(
            f"{len(seq_ids)} sequence{'s' if len(seq_ids) != 1 else ''} — enter organism and location details for each row.",
            size="sm",
            c="dimmed",
        ),
    ]
    if len(seq_ids) > 1:
        children.append(
            dmc.Alert(
                "Annotations in your GFF file are matched to FASTA headers automatically.",
                color="var(--mantine-color-blue-6)",
                variant="light",
                title="Multiple sequences",
            )
        )

    children.append(
        html.Div(
            dmc.Table(
                [
                    html.Thead(
                        html.Tr(
                            [
                                html.Th("Sequence"),
                                html.Th("Genus"),
                                html.Th("Species"),
                                html.Th("Strain"),
                                html.Th("Genome Accession"),
                                html.Th("Host contig"),
                                html.Th("Start"),
                                html.Th("End"),
                                html.Th("Strand"),
                            ]
                        )
                    ),
                    html.Tbody(rows),
                ],
                striped=True,
                highlightOnHover=True,
                withTableBorder=True,
                withColumnBorders=True,
            ),
            style={"overflowX": "auto"},
        )
    )

    return dmc.Stack(children, gap="md")


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

# Section: Location (table — one row per FASTA header)
location_section = dmc.Stack(
    [
        html.Div(
            id="submit-location-fields",
            children=create_location_placeholder(),
        ),
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

sidebar = dmc.Box(
    dmc.Stack(
        [data_policy_card, create_submission_queue(max_items=15)],
        gap="md",
    ),
    style={
        "position": "sticky",
        "top": "var(--mantine-spacing-md)",
        "maxHeight": "calc(100vh - 2rem)",
        "overflowY": "auto",
    },
)

layout = dmc.Container(
    size="lg",
    children=[
        dcc.Location(id="submit-url", refresh=False),
        dcc.Store(id="submit-fasta-prefill"),
        dcc.Store(id="submit-primary-seq-ids", data=[]),
        submission_header,
        dmc.Grid(
            [
                dmc.GridCol(
                    dmc.Stack(
                        [submission_info_card, submission_body],
                        gap="xl",
                    ),
                    span={"base": "auto", "md": 8},
                ),
                dmc.GridCol(
                    sidebar,
                    span={"base": "auto", "md": 4},
                ),
            ],
        ),
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

        # Enrich classification with family/navis/haplotype from main DB when closest_match is set
        if classification:
            classification = _enrich_classification_from_db(classification)

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
    primary_loc_genus,
    primary_loc_species,
    primary_loc_strain,
    primary_loc_assembly,
    primary_loc_hostchr,
    primary_loc_shipstart,
    primary_loc_shipend,
    primary_loc_strand,
    uploader,
    evidence,
    comment,
):
    """Collect and validate all Starship entries from the FASTA table."""
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

    if not primary_loc_hostchr or len(primary_loc_hostchr) != expected_count:
        raise WebValidationError(
            "Provide organism and location details for each sequence in the table below.",
            "location",
        )

    def _cell(values, index):
        if not values or index >= len(values):
            return None
        value = values[index]
        if isinstance(value, str):
            value = value.strip()
        return value or None

    locations = []
    for index in range(expected_count):
        strand_val = (
            primary_loc_strand[index]
            if primary_loc_strand and index < len(primary_loc_strand)
            else 1
        )
        if isinstance(strand_val, str):
            strand_val = int(strand_val) if strand_val else 1

        locations.append(
            {
                "genus": _cell(primary_loc_genus, index),
                "species": _cell(primary_loc_species, index),
                "strain": _cell(primary_loc_strain, index),
                "assembly_accession": _cell(primary_loc_assembly, index),
                "hostchr": primary_loc_hostchr[index],
                "shipstart": primary_loc_shipstart[index],
                "shipend": primary_loc_shipend[index],
                "strand_radio": strand_val or 1,
            }
        )

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


def process_validated_submission_entry(validated_data):
    """Queue one validated submission entry for async processing."""
    submission_id = str(uuid.uuid4())
    shipstrand = "+" if (validated_data.get("strand_radio") or 1) == 1 else "-"
    seq_contents = validated_data.get("seq_contents")
    seq_filename = validated_data.get("seq_filename") or validated_data.get("filename")
    seq_date = validated_data.get("seq_date") or datetime.datetime.now().timestamp()

    try:
        db_submission_id = insert_submission(
            seq_contents=seq_contents,
            seq_filename=seq_filename,
            seq_date=seq_date,
            anno_contents=validated_data.get("anno_contents"),
            anno_filename=validated_data.get("anno_filename"),
            anno_date=validated_data.get("anno_date"),
            uploader=validated_data["uploader"],
            evidence=validated_data["evidence"],
            genus=validated_data["genus"],
            species=validated_data["species"],
            strain=validated_data.get("strain"),
            hostchr=validated_data["hostchr"],
            shipstart=validated_data["shipstart"],
            shipend=validated_data["shipend"],
            shipstrand=shipstrand,
            comment=validated_data.get("comment") or "",
            accession=None,
            needs_review=True,
            classification=validated_data.get("classification"),
            assembly_accession=validated_data.get("assembly_accession"),
        )
        validated_data["db_submission_id"] = db_submission_id
    except Exception as e:
        logger.warning(f"Failed to insert submission into queue (continuing): {e}")
        validated_data["db_submission_id"] = None

    create_submission_status(submission_id, "queued")
    task_result = run_task(process_submission_task, validated_data, submission_id)

    accession_assigned = None
    if hasattr(task_result, "id"):
        update_submission_status(
            submission_id,
            "processing",
            progress=10,
            message="Your submission is in the queue. Processing usually takes a few minutes.",
        )
        cache.set(f"task:{submission_id}", task_result.id, timeout=3600)
    else:
        update_submission_status(
            submission_id,
            "completed" if task_result.get("success") else "failed",
            progress=100,
            message=task_result.get("message", "Processing complete"),
            result=task_result,
        )
        if task_result.get("success"):
            accession_assigned = task_result.get("accession")

    return submission_id, task_result, accession_assigned


def build_submission_success_message(entries, results):
    """Build modal feedback for one or more queued submissions."""
    _submission_ids = [item["submission_id"] for item in results]
    async_pending = any(item["async"] for item in results)

    if len(entries) == 1:
        result = results[0]
        submission_id = result["submission_id"]
        task_result = result["task_result"]
        if async_pending:
            return dmc.Alert(
                children=[
                    dmc.Text(
                        "Your submission is in the queue and will be processed shortly.",
                        fw=600,
                        size="lg",
                    ),
                    html.Br(),
                    dmc.Text(f"Submission ID: {submission_id}", size="sm", c="dimmed"),
                    html.Br(),
                    dmc.Text(
                        "We'll email you when processing is complete. Our curation team will then review your submission.",
                        size="sm",
                    ),
                ],
                title="Submission received",
                color="var(--mantine-color-green-6)",
                variant="light",
            )
        if task_result.get("success"):
            return dmc.Alert(
                children=[
                    dmc.Text(
                        f"Accession: {task_result['accession']}",
                        size="md",
                        fw=500,
                        c="var(--mantine-color-green-6)",
                    ),
                    dmc.Text(
                        f"StarshipID: {task_result['filename']}",
                        size="sm",
                        c="dimmed",
                    ),
                    html.Br(),
                    dmc.Text(
                        "The submission still has to be curated in order to be visible on starbase. You may be contacted via email if we need some more information.",
                        size="sm",
                    ),
                ],
                title="Submission complete",
                color="var(--mantine-color-green-6)",
                variant="light",
            )
        return dmc.Alert(
            children=[
                dmc.Text(
                    task_result.get(
                        "user_message",
                        "We couldn't complete your submission. Please check your input and try again.",
                    ),
                    size="sm",
                ),
            ],
            title="Submission failed",
            color="var(--mantine-color-red-6)",
            variant="light",
        )

    list_items = []
    for entry, result in zip(entries, results):
        label = entry.get("filename") or entry.get("seq_id")
        list_items.append(dmc.ListItem(f"{label} — ID {result['submission_id'][:8]}…"))

    if async_pending:
        return dmc.Alert(
            children=[
                dmc.Text(
                    f"Queued {len(entries)} Starship submissions.",
                    fw=600,
                    size="lg",
                ),
                dmc.List(list_items, size="sm", spacing="xs"),
                dmc.Text(
                    "We'll email you when processing is complete. Our curation team will then review each submission.",
                    size="sm",
                ),
            ],
            title="Submissions received",
            color="var(--mantine-color-green-6)",
            variant="light",
        )

    return dmc.Alert(
        children=[
            dmc.Text(
                f"Processed {len(entries)} Starship submissions.",
                fw=600,
                size="lg",
            ),
            dmc.List(list_items, size="sm", spacing="xs"),
        ],
        title="Submissions complete",
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
        Output("submit-location-fields", "children", allow_duplicate=True),
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
        State({"type": "primary-loc-genus", "index": ALL}, "value"),
        State({"type": "primary-loc-species", "index": ALL}, "value"),
        State({"type": "primary-loc-strain", "index": ALL}, "value"),
        State({"type": "primary-loc-assembly", "index": ALL}, "value"),
        State({"type": "primary-loc-hostchr", "index": ALL}, "value"),
        State({"type": "primary-loc-shipstart", "index": ALL}, "value"),
        State({"type": "primary-loc-shipend", "index": ALL}, "value"),
        State({"type": "primary-loc-strand", "index": ALL}, "value"),
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
    primary_loc_genus,
    primary_loc_species,
    primary_loc_strain,
    primary_loc_assembly,
    primary_loc_hostchr,
    primary_loc_shipstart,
    primary_loc_shipend,
    primary_loc_strand,
    uploader,
    evidence,
    comment,
):
    """Validate and queue one or more Starship submissions."""
    triggered_id = ctx.triggered_id if ctx.triggered else None
    no_update_tail = (dash.no_update,) * 5

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
            primary_loc_genus=primary_loc_genus,
            primary_loc_species=primary_loc_species,
            primary_loc_strain=primary_loc_strain,
            primary_loc_assembly=primary_loc_assembly,
            primary_loc_hostchr=primary_loc_hostchr,
            primary_loc_shipstart=primary_loc_shipstart,
            primary_loc_shipend=primary_loc_shipend,
            primary_loc_strand=primary_loc_strand,
            uploader=uploader,
            evidence=evidence,
            comment=comment,
        )

        results = []
        for entry in entries:
            submission_id, task_result, accession_assigned = (
                process_validated_submission_entry(entry)
            )
            results.append(
                {
                    "submission_id": submission_id,
                    "task_result": task_result,
                    "accession_assigned": accession_assigned,
                    "async": hasattr(task_result, "id"),
                }
            )
            try:
                send_curator_notification(submission_id, entry, accession_assigned)
                if uploader:
                    send_submission_confirmation(
                        uploader, submission_id, accession_assigned
                    )
            except Exception as e:
                logger.error(f"Failed to send email notifications: {str(e)}")

        message = build_submission_success_message(entries, results)
        reset_form = all(
            item["async"] or item["task_result"].get("success") for item in results
        )

        if reset_form:
            return (
                True,
                message,
                "",
                "",
                "",
                [],
                create_location_placeholder(),
            )

        return (
            True,
            message,
            *no_update_tail,
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
        Output("submit-location-fields", "children"),
        Output("submit-primary-seq-ids", "data"),
    ],
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-prefill", "data"),
    ],
)
def update_primary_location_fields(seq_contents, seq_filename, fasta_prefill):
    """Render location table with one row per FASTA header."""
    resolved_contents, resolved_filename, _seq_date = _resolve_primary_fasta(
        fasta_prefill, seq_contents, seq_filename, None
    )
    if not resolved_contents:
        return create_location_placeholder(), []

    try:
        records = parse_fasta_records(
            resolved_contents, resolved_filename or "upload.fa"
        )
    except WebValidationError:
        return create_location_placeholder(), []

    seq_ids = [record["id"] for record in records]
    return create_location_table(seq_ids), seq_ids


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

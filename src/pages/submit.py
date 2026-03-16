import base64
import dash
from dash_iconify import DashIconify
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback, ctx
from dash.exceptions import PreventUpdate

import datetime
from typing import Dict, Any, Optional
from src.utils.seq_utils import parse_fasta, parse_gff
import re
from src.database.sql_engine import get_submissions_session
from src.database.models.schema import Submission
from sqlalchemy.exc import SQLAlchemyError
from src.components.ui import create_file_upload

from src.config.logging import get_logger
from src.config.cache import cache
from src.config.celery_config import run_task
from src.tasks import process_submission_task
from src.utils.web_submission_adapter import (
    validate_submission_data,
    WebValidationError,
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


dash.register_page(__name__)

submission_header = dmc.Title(
    [
        "Submission of Candidate ",
        html.Span("Starship", style={"fontStyle": "italic"}),
        " Sequences to ",
        html.Span("starbase", className="logo-text"),
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
                dmc.Text(
                    "Each submission is processed by our automated pipeline, then manually reviewed by our curation team.",
                    # TODO: "You'll receive a confirmation email once your submission is processed."
                    size="md",
                    mb="md",
                ),
                dmc.List(
                    [
                        dmc.ListItem(dmc.Text("Sequence validation and parsing")),
                        dmc.ListItem(dmc.Text("Duplicate checking and classification")),
                        dmc.ListItem(dmc.Text("Accession number assignment")),
                    ],
                    type="ordered",
                    size="sm",
                    spacing="xs",
                    icon=dmc.ThemeIcon(
                        DashIconify(icon="tabler:check", width=16),
                        size="sm",
                        variant="light",
                        color="var(--mantine-color-indigo-6)",
                    ),
                ),
                dmc.Text(
                    "Each submission will also go through a manual review. Submissions will be included in the next database release.",
                    size="sm",
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
    title="Submission",
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
                                children=html.Div(
                                    id="loading-output-1",
                                    style={"minHeight": "60px"},
                                ),
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
                                children=html.Div(
                                    id="loading-output-2",
                                    style={"minHeight": "60px"},
                                ),
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

# Section: Contact & Organism (two-column)
contact_organism_section = dmc.Stack(
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
                    span={"base": 12, "md": 6},
                ),
                dmc.GridCol(
                    dmc.Stack(
                        [
                            dmc.Text(
                                "Organism",
                                size="sm",
                                fw=600,
                                c="dimmed",
                                tt="uppercase",
                            ),
                            dmc.Group(
                                [
                                    dmc.TextInput(
                                        id="genus",
                                        label="Genus",
                                        placeholder="e.g., Alternaria",
                                        required=True,
                                        style={"flex": 1},
                                    ),
                                    dmc.TextInput(
                                        id="species",
                                        label="Species",
                                        placeholder="e.g., alternata",
                                        required=True,
                                        style={"flex": 1},
                                    ),
                                ],
                            ),
                        ],
                        gap="md",
                    ),
                    span={"base": 12, "md": 6},
                ),
            ],
        ),
    ],
    gap="md",
)

# Section: Location
location_section = dmc.Stack(
    [
        dmc.Text(
            "Location in host genome",
            size="sm",
            fw=600,
            c="dimmed",
            tt="uppercase",
        ),
        dmc.TextInput(
            id="hostchr",
            label="Host contig or scaffold ID",
            placeholder="e.g., chr1, scaffold_001, NC_123456",
            description="The identifier from your genome assembly where this Starship was found",
            required=True,
        ),
        dmc.Grid(
            [
                dmc.GridCol(
                    dmc.NumberInput(
                        id="shipstart",
                        label="Start coordinate",
                        placeholder="e.g., 1200",
                        description="5' boundary of the Starship in the host assembly",
                        required=True,
                        min=1,
                        step=1,
                    ),
                    span={"base": 12, "sm": 6},
                ),
                dmc.GridCol(
                    dmc.NumberInput(
                        id="shipend",
                        label="End coordinate",
                        placeholder="e.g., 20500",
                        description="3' boundary of the Starship in the host assembly",
                        required=True,
                        min=1,
                        step=1,
                    ),
                    span={"base": 12, "sm": 6},
                ),
            ],
        ),
        dmc.RadioGroup(
            id="strand-radios",
            label="Strand",
            value=1,
            children=[
                dmc.Radio(
                    label="Positive strand",
                    value=1,
                    color="indigo",
                ),
                dmc.Space(h="sm"),
                dmc.Radio(
                    label="Negative strand",
                    value=2,
                    color="indigo",
                ),
            ],
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
                    dmc.Text("Add optional notes", size="sm", c="dimmed"),
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
                    "Submit Starship",
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
            contact_organism_section,
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

layout = dmc.Container(
    size="lg",
    children=[
        dcc.Location(id="submit-url", refresh=False),
        dcc.Store(id="submit-fasta-prefill"),
        submission_header,
        dmc.Grid(
            [
                dmc.GridCol(
                    submission_info_card,
                    span={"base": 12, "md": 6},
                ),
                dmc.GridCol(
                    dmc.Stack(
                        [data_policy_card, create_submission_queue(max_items=15)],
                        gap="md",
                    ),
                    span={"base": 12, "md": 6},
                ),
            ],
        ),
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

    if not fasta_contents:
        info_banner = dmc.Alert(
            "We pre-filled metadata from your BLAST search, but the sequence file wasn't found. Upload your FASTA file below.",
            title="Partial pre-fill from BLAST",
            color="var(--mantine-color-orange-6)",
            variant="light",
            mb="sm",
        )
        return None, comment, "BLAST", info_banner

    try:
        _ct, content_string = fasta_contents.split(",", 1)
        fasta_text = base64.b64decode(content_string).decode("utf-8")
        records = parse_fasta(fasta_text, fasta_filename)
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

    return (
        {"contents": fasta_contents, "filename": fasta_filename},
        comment,
        "BLAST",
        info_banner,
    )


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
):
    try:
        _ct, content_string = seq_contents.split(",", 1)
        seq_decoded = base64.b64decode(content_string).decode("utf-8")
        seq_datetime_str = datetime.datetime.fromtimestamp(seq_date).strftime(
            "%Y-%m-%d %H:%M:%S"
        )

        anno_decoded = None
        anno_filename_val = None
        anno_datetime_str = None
        if anno_contents:
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
                hostchr=hostchr,
                shipstart=shipstart,
                shipend=shipend,
                shipstrand=shipstrand,
                comment=comment,
                accession_tag=accession,
                needs_review=needs_review,
            )
            session.add(submission)
            session.commit()
            logger.debug(f"Successfully inserted submission for {seq_filename}")
            return True

    except SQLAlchemyError as e:
        logger.error(f"Database error during submission: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error processing submission: {str(e)}")
        raise


def create_fasta_display(records, filename):
    """Create HTML components to display FASTA file info."""
    return dmc.Alert(
        children=[
            dmc.Text(f"✓ File: {filename}", size="sm", fw=500),
            dmc.Text(f"Sequences found: {len(records)}", size="sm", c="dimmed"),
        ],
        title="FASTA loaded",
        color="var(--mantine-color-green-6)",
        variant="light",
    )


@callback(
    Output("submit-ship", "disabled"),
    Input("submit-consent-checkbox", "checked"),
)
def toggle_submit_button(checked):
    return not bool(checked)


@callback(
    [
        Output("submit-modal", "opened"),
        Output("output-data-upload", "children"),
        Output("submit-ship", "loading"),
        Output("uploader", "value"),
        Output("evidence", "value"),
        Output("genus", "value"),
        Output("species", "value"),
        Output("hostchr", "value"),
        Output("shipstart", "value"),
        Output("shipend", "value"),
        Output("comment", "value"),
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
        State("uploader", "value"),
        State("evidence", "value"),
        State("genus", "value"),
        State("species", "value"),
        State("hostchr", "value"),
        State("shipstart", "value"),
        State("shipend", "value"),
        State("strand-radios", "value"),
        State("comment", "value"),
    ],
    prevent_initial_call=True,
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
    uploader,
    evidence,
    genus,
    species,
    hostchr,
    shipstart,
    shipend,
    strand_radio,
    comment,
):
    """
    Main submission callback with asynchronous processing.

    This function validates input, starts an async task, and provides immediate feedback.
    Users can check progress via submission status.
    """
    # Determine which button was clicked
    triggered_id = ctx.triggered_id if ctx.triggered else None

    # Handle modal close button
    if triggered_id == "close":
        return (
            False,
            "",
            False,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
        )

    # If submit button wasn't clicked, prevent update
    if triggered_id != "submit-ship":
        raise PreventUpdate

    loading = True

    # Use prefilled FASTA from BLAST session if no file was uploaded directly
    if not seq_contents and fasta_prefill:
        seq_contents = fasta_prefill.get("contents")
        seq_filename = seq_filename or fasta_prefill.get(
            "filename", "sequence_from_blast.fa"
        )
        seq_date = seq_date or datetime.datetime.now().timestamp()
    # Early validation - email format
    if uploader and not validate_email(uploader):
        return (
            True,
            dmc.Alert(
                "Use a valid email address (e.g., name@example.com) so we can send your confirmation.",
                title="Invalid email",
                color="var(--mantine-color-orange-6)",
                variant="light",
            ),
            False,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
        )

    try:
        # Step 1: Validate input data (quick validation)
        # This must pass before we proceed with submission
        # If validation fails, validate_submission_data will raise ValidationError
        logger.debug("Validating submission data")
        validated_data = validate_submission_data(
            seq_contents,
            seq_filename,
            uploader,
            evidence,
            genus,
            species,
            hostchr,
            shipstart,
            shipend,
        )

        # Step 2: Create submission data package, using `validated_data` output, which is already a dict

        validated_data["seq_date"] = seq_date
        validated_data["anno_contents"] = anno_contents
        validated_data["anno_filename"] = anno_filename
        validated_data["anno_date"] = anno_date
        validated_data["strand_radio"] = strand_radio
        validated_data["comment"] = comment or ""

        # Step 3: Generate unique submission ID
        submission_id = str(uuid.uuid4())
        logger.info(f"Created submission {submission_id} for file {seq_filename}")

        # Step 4: Create initial status
        create_submission_status(submission_id, "queued")

        # Step 5: Start async processing task
        logger.debug(f"Starting async submission processing for {submission_id}")
        task_result = run_task(process_submission_task, validated_data, submission_id)

        # Store task ID for status checking
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
            # Task ran synchronously (no Celery)
            update_submission_status(
                submission_id,
                "completed" if task_result.get("success") else "failed",
                progress=100,
                message=task_result.get("message", "Processing complete"),
                result=task_result,
            )
            if task_result.get("success"):
                accession_assigned = task_result.get("accession")

        # Step 6: Send email notifications
        try:
            # Send curator notification
            send_curator_notification(submission_id, validated_data, accession_assigned)
            # Send confirmation to submitter
            if uploader:
                send_submission_confirmation(
                    uploader, submission_id, accession_assigned
                )
        except Exception as e:
            logger.error(f"Failed to send email notifications: {str(e)}")

        # Step 7: Show immediate feedback
        loading = False
        modal_open = True

        if hasattr(task_result, "id"):
            # Async processing
            message = dmc.Alert(
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
            # Reset form
            return (
                modal_open,
                message,
                loading,
                "",
                "",
                "",
                "",
                "",
                None,
                None,
                "",  # Reset form fields
            )
        else:
            # Sync processing completed
            if task_result.get("success"):
                message = dmc.Alert(
                    children=[
                        dmc.Text("Submission complete.", fw=600, size="lg"),
                        html.Br(),
                        dmc.Text(
                            f"Accession: {task_result['accession']}",
                            size="md",
                            fw=500,
                            c="var(--mantine-color-green-6)",
                        ),
                        dmc.Text(
                            f"Filename: {task_result['filename']}",
                            size="sm",
                            c="dimmed",
                        ),
                        html.Br(),
                        dmc.Text(
                            "Our curation team will review it. You'll receive an email confirmation.",
                            size="sm",
                        ),
                    ],
                    title="Submission complete",
                    color="var(--mantine-color-green-6)",
                    variant="light",
                )
                # Reset form on success
                return (
                    modal_open,
                    message,
                    loading,
                    "",
                    "",
                    "",
                    "",
                    "",
                    None,
                    None,
                    "",  # Reset form fields
                )
            else:
                message = dmc.Alert(
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
                # Keep form filled on error
                return (
                    modal_open,
                    message,
                    loading,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                    dash.no_update,
                )

    except WebValidationError as e:
        # Handle validation errors from web_submission_adapter
        # Convert to our ValidationError for consistent handling
        validation_error = ValidationError(
            str(e.message) if hasattr(e, "message") else str(e),
            getattr(e, "field", None),
        )
        error_response = handle_submission_error(validation_error)
        loading = False
        return (
            error_response["modal_open"],
            error_response["message"],
            error_response["loading"],
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
        )

    except SubmissionError as e:
        # Handle our custom submission errors
        error_response = handle_submission_error(e)
        return (
            error_response["modal_open"],
            error_response["message"],
            error_response["loading"],
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
        )

    except Exception as e:
        # Handle unexpected errors
        logger.error(f"Unexpected error in submit_ship: {str(e)}", exc_info=True)
        error_response = handle_submission_error(e)
        return (
            error_response["modal_open"],
            error_response["message"],
            error_response["loading"],
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
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
    Output("submit-fasta-sequence-upload", "children"),
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-prefill", "data"),
    ],
)
def update_fasta_details(seq_contents, seq_filename, fasta_prefill):
    # Manually uploaded file always takes priority
    if seq_contents is not None:
        try:
            content_type, content_string = seq_contents.split(",")
            query_string = base64.b64decode(content_string).decode("utf-8")
            records = parse_fasta(query_string, seq_filename)
            return create_fasta_display(records, seq_filename)
        except Exception as e:
            logger.error(e)
            return dmc.Alert(
                "We couldn't parse this FASTA file. Check that it's valid and try again.",
                title="File error",
                color="var(--mantine-color-red-6)",
                variant="light",
            )

    # Show prefill details when no file has been manually uploaded
    if fasta_prefill:
        prefill_contents = fasta_prefill.get("contents")
        prefill_filename = fasta_prefill.get("filename", "sequence_from_blast.fa")
        if prefill_contents:
            try:
                _ct, content_string = prefill_contents.split(",", 1)
                query_string = base64.b64decode(content_string).decode("utf-8")
                records = parse_fasta(query_string, prefill_filename)
                return dmc.Stack(
                    [
                        dmc.Text(f"File name: {prefill_filename}", size="sm", fw=500),
                        dmc.Text(f"Number of sequences: {len(records)}", size="sm"),
                        dmc.Text(
                            "Pre-filled from BLAST. Upload a different file to replace.",
                            size="xs",
                            c="dimmed",
                        ),
                    ],
                    gap="xs",
                )
            except Exception as e:
                logger.error(e)

    return dmc.Text(
        "No file selected. Choose a FASTA file above to continue.",
        size="sm",
        c="dimmed",
    )


@callback(
    Output("submit-output-gff-upload", "children"),
    [
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
    ],
)
def update_gff_details(anno_contents, anno_filename):
    if anno_contents is None:
        return dmc.Text(
            "No file selected. GFF annotations are optional.",
            size="sm",
            c="dimmed",
        )
    try:
        children = parse_gff(anno_contents, anno_filename)
        if isinstance(children, (list, tuple)) and len(children) > 0:
            return dmc.Alert(
                children=[
                    dmc.Text(f"✓ File: {anno_filename}", size="sm", fw=500),
                    *(
                        [children]
                        if not isinstance(children, (list, tuple))
                        else children
                    ),
                ],
                title="GFF loaded",
                color="var(--mantine-color-green-6)",
                variant="light",
            )
        return children
    except Exception as e:
        logger.error(e)
        return dmc.Alert(
            "We couldn't parse this GFF file. Ensure it's valid GFF3 format and try again.",
            title="File error",
            color="var(--mantine-color-red-6)",
            variant="light",
        )

import base64
import dash
from dash_iconify import DashIconify
import dash_bootstrap_components as dbc
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
from src.components.leaderboard import create_leaderboard
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
                title="Validation Error",
                children=str(error.user_message),
                color="var(--mantine-color-yellow-6)",
                variant="light",
            ),
            "loading": False,
        }
    elif isinstance(error, ProcessingError):
        return {
            "modal_open": False,
            "message": dmc.Alert(
                title="Processing Error",
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
                title="Database Error",
                children=dmc.Stack(
                    [
                        dmc.Text(str(error.user_message)),
                        dmc.Text(
                            "If this persists, please contact support.",
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
                title="Unexpected Error",
                children=dmc.Stack(
                    [
                        dmc.Text("An unexpected error occurred. Please try again."),
                        dmc.Text(
                            "If this persists, please contact support.",
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


layout = dmc.Container(
    size="lg",
    children=[
        dcc.Location(id="submit-url", refresh=False),
        dcc.Store(id="submit-fasta-prefill"),
        # Header Section
        dmc.Title(
            [
                "Submission of Candidate ",
                html.Span("Starship", style={"fontStyle": "italic"}),
                " Sequences to ",
                html.Span("starbase", className="logo-text"),
            ],
            order=1,
            mb="md",
        ),
        dmc.Paper(
            children=[
                dmc.Text(
                    "Comparative genomics projects are a collaborative effort. Submit your Starship discoveries to the community to help us build the most comprehensive database of Starship elements.",
                    size="md",
                    c="dimmed",
                    mb="md",
                ),
                dmc.Text(
                    "Each submission is processed by our automated pipeline, then manually reviewed by our curation team.",
                    # TODO: "You'll receive a confirmation email once your submission is processed."
                    size="md",
                    c="dimmed",
                    mb="md",
                ),
                dmc.List(
                    [
                        dmc.ListItem(dmc.Text("Sequence validation and parsing")),
                        dmc.ListItem(dmc.Text("Duplicate checking and classification")),
                        dmc.ListItem(dmc.Text("Accession number assignment")),
                        dmc.ListItem(dmc.Text("Curator review")),
                        dmc.ListItem(
                            dmc.Text("Inclusion in next database release (if approved)")
                        ),
                    ],
                    type="ordered",
                    size="sm",
                    spacing="xs",
                    icon=dmc.ThemeIcon(
                        DashIconify(icon="tabler:check", width=16),
                        size="sm",
                        variant="light",
                        color="blue",
                    ),
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
        ),
        # Leaderboard Section
        create_leaderboard(limit=10, title="🏆 Top Contributors"),
        dmc.Space(h=20),
        # Pending Submissions Queue
        create_submission_queue(max_items=15),
        dmc.Space(h=20),
        # Modal
        dbc.Modal(
            [
                dbc.ModalHeader(
                    dbc.ModalTitle("Submission Received", className="text-info"),
                    close_button=True,
                ),
                dbc.ModalBody(
                    [
                        html.Div(
                            [
                                html.Div(id="output-data-upload", className="mt-3"),
                            ],
                            className="text-center",
                        )
                    ]
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close",
                        id="close",
                        className="ms-auto",
                        color="primary",
                        n_clicks=0,
                    )
                ),
            ],
            id="submit-modal",
            is_open=False,
            centered=True,
        ),
        dmc.Alert(
            "Fields marked with * are required",
            color="red",
            variant="light",
            title="Submission Requirements",
        ),
        dmc.Space(h=10),
        # Form Content
        dmc.Paper(
            style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
            children=dbc.Form(
                [
                    dmc.Stack(
                        [
                            # File Upload Section
                            dmc.Stack(
                                [
                                    dmc.Title("Upload Files", order=2, mb="md"),
                                    # FASTA Upload
                                    html.Div(id="submit-prefill-info"),
                                    dmc.Paper(
                                        children=[
                                            dmc.Text(
                                                [
                                                    "Starship Sequence ",
                                                    html.Span(
                                                        "*",
                                                        style={
                                                            "color": "var(--mantine-color-red-6)"
                                                        },
                                                    ),
                                                ],
                                                fw=500,
                                                mb="sm",
                                            ),
                                            create_file_upload(
                                                upload_id="submit-fasta-upload",
                                                output_id="submit-fasta-sequence-upload",
                                                accept_types=[
                                                    ".fa",
                                                    ".fas",
                                                    ".fasta",
                                                    ".fna",
                                                ],
                                                placeholder_text="Select a FASTA file to upload",
                                            ),
                                            dcc.Loading(
                                                id="loading-1",
                                                type="circle",
                                                children=html.Div(
                                                    id="loading-output-1"
                                                ),
                                            ),
                                        ],
                                        p="md",
                                        radius="md",
                                        withBorder=False,
                                    ),
                                    # GFF Upload
                                    dmc.Paper(
                                        children=[
                                            dmc.Text(
                                                "Gene Annotations (GFF3)",
                                                fw=500,
                                                mb="sm",
                                            ),
                                            create_file_upload(
                                                upload_id="submit-upload-gff",
                                                output_id="submit-output-gff-upload",
                                                accept_types=[".gff", ".gff3", ".tsv"],
                                                placeholder_text="Select a GFF file to upload",
                                            ),
                                            dcc.Loading(
                                                id="loading-2",
                                                type="circle",
                                                children=html.Div(
                                                    id="loading-output-2"
                                                ),
                                            ),
                                        ],
                                        p="md",
                                        radius="md",
                                        withBorder=False,
                                    ),
                                ],
                                gap="lg",
                            ),
                            # Metadata Section
                            dmc.Stack(
                                [
                                    dmc.Title(
                                        "Metadata", order=2, mb="md"
                                    ),
                                    # Curator Info
                                    dmc.TextInput(
                                        id="uploader",
                                        label="Email of curator",
                                        placeholder="Enter email",
                                        required=True,
                                        leftSection=DashIconify(icon="fas fa-envelope"),
                                    ),
                                    html.Div(
                                        id="email-validation-message",
                                        style={"minHeight": "20px"},
                                    ),
                                    dmc.TextInput(
                                        id="evidence",
                                        label="How were Starships annotated?",
                                        placeholder="e.g., starfish, manual curation, etc.",
                                        required=True,
                                        description="Please specify the tool or method used",
                                    ),
                                    # Taxonomy Info
                                    dmc.Group(
                                        [
                                            dmc.TextInput(
                                                id="genus",
                                                label="Genus",
                                                placeholder="Alternaria",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                            dmc.TextInput(
                                                id="species",
                                                label="Species",
                                                placeholder="alternata",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                        ]
                                    ),
                                    # Location Info
                                    dmc.TextInput(
                                        id="hostchr",
                                        label="Host genome contig/scaffold/chromosome ID",
                                        placeholder="e.g., chr1, scaffold_001, or GenBank Accession",
                                        description="Identifier from the host genome assembly",
                                        required=True,
                                    ),
                                    dmc.Group(
                                        [
                                            dmc.NumberInput(
                                                id="shipstart",
                                                label="Start coordinate",
                                                placeholder="1200",
                                                description="5' boundary position",
                                                required=True,
                                                min=1,
                                                step=1,
                                                style={"flex": 1},
                                            ),
                                            dmc.NumberInput(
                                                id="shipend",
                                                label="End coordinate",
                                                placeholder="20500",
                                                description="3' boundary position",
                                                required=True,
                                                min=1,
                                                step=1,
                                                style={"flex": 1},
                                            ),
                                        ]
                                    ),
                                    dmc.RadioGroup(
                                        id="strand-radios",
                                        label="Strand",
                                        value=1,
                                        children=[
                                            dmc.Radio(label="Positive strand", value=1, color="indigo"),
                                            dmc.Space(h="sm"),
                                            dmc.Radio(label="Negative strand", value=2, color="indigo"),
                                        ],
                                    ),
                                    # Comments
                                    dmc.Textarea(
                                        id="comment",
                                        label="Additional information (optional)",
                                        placeholder="Any comments about the Starship features, annotations, or host genome?",
                                        description="Share any relevant details about this submission",
                                        minRows=3,
                                        autosize=True,
                                        maxRows=6,
                                    ),
                                ],
                                gap="md",
                            ),
                            # Submit Button
                            dmc.Center(
                                dmc.Button(
                                    "Submit Starship",
                                    id="submit-ship",
                                    size="lg",
                                    variant="filled",
                                    color="indigo",
                                    loading=False,
                                    leftSection=html.I(className="fas fa-rocket")
                                    if False
                                    else None,  # Icon placeholder
                                    fullWidth=False,
                                    style={"marginTop": "1rem"},
                                ),
                            ),
                        ],
                        gap="xl",
                    ),
                ]
            ),
            p="xl",
            radius="md",
            withBorder=True,
        ),
        dmc.Space(h="xl"),
        # Submission Status Section
        dmc.Paper(
            style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
            children=[
                dmc.Title(
                    "Submission Status", order=3, mb="md"
                ),
                dmc.Stack(
                    [
                        dmc.TextInput(
                            id="status-submission-id",
                            label="Check Submission Status",
                            placeholder="Enter submission ID",
                            description="Enter the submission ID to check processing status",
                        ),
                        dmc.Button(
                            "Check Status",
                            id="check-status-btn",
                            variant="light",
                            size="sm",
                            color="indigo",
                            fullWidth=False,
                            radius="md",
                        ),
                        html.Div(id="status-display", className="mt-3"),
                    ],
                    gap="sm",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
        ),
        # Modal
        dbc.Modal(
            [
                dbc.ModalHeader(
                    dbc.ModalTitle("Submission Status"),
                    close_button=True,
                ),
                dbc.ModalBody(
                    [
                        html.Div(
                            [
                                html.Div(id="output-data-upload", className="mt-3"),
                            ],
                            className="text-center",
                        )
                    ]
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close",
                        id="close",
                        className="ms-auto",
                        color="primary",
                        n_clicks=0,
                    )
                ),
            ],
            id="submit-modal",
            is_open=False,
            centered=True,
        ),
    ],
    style={
        "margin": "0 auto",
        "padding": "var(--mantine-spacing-xl)",
    },
)


@callback(
    [
        Output("submit-fasta-prefill", "data"),
        Output("comment", "value"),
        Output("evidence", "value"),
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
            "Classification metadata pre-filled from BLAST, but no FASTA file was found. "
            "Please upload the sequence manually.",
            title="Partial Pre-fill from BLAST",
            color="yellow",
            variant="light",
            mb="sm",
        )
        return None, comment, "Starbase BLAST", info_banner

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
                    f"Sequence pre-filled from BLAST: {fasta_filename}"
                    + (f" ({n_seqs} sequence{'s' if n_seqs != 1 else ''})" if n_seqs else ""),
                    fw=500,
                ),
                dmc.Text(
                    "You can upload a different file below to override it.",
                    size="sm",
                    c="dimmed",
                ),
            ],
            gap="xs",
        ),
        title="Pre-filled from BLAST",
        color="green",
        variant="light",
        mb="sm",
    )

    return (
        {"contents": fasta_contents, "filename": fasta_filename},
        comment,
        "Starbase BLAST",
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
        title="FASTA File Loaded Successfully",
        color="green",
        variant="light",
    )


@callback(
    [
        Output("submit-modal", "is_open"),
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
        State("submit-modal", "is_open"),
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
        State("submit-modal", "is_open"),
    ],
    prevent_initial_call=True,
)
def submit_ship(
    submit_clicks,
    close_clicks,
    is_open,
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
        seq_filename = seq_filename or fasta_prefill.get("filename", "sequence_from_blast.fa")
        seq_date = seq_date or datetime.datetime.now().timestamp()
    # Early validation - email format
    if uploader and not validate_email(uploader):
        return (
            True,
            html.Div(
                [
                    html.H4("Validation Error", className="text-warning"),
                    html.P(
                        "Please enter a valid email address.", className="text-muted"
                    ),
                ]
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
                message="Submission is being processed...",
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

        # Step 6: Show immediate feedback with submission ID
        modal = not is_open if not close_modal else False
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
            message = html.Div(
                [
                    dmc.Alert(
                        children=[
                            dmc.Text(
                                "✓ Your submission has been received and is being processed!",
                                fw=600,
                                size="lg",
                            ),
                            html.Br(),
                            dmc.Text(
                                f"Submission ID: {submission_id}", size="sm", c="dimmed"
                            ),
                            html.Br(),
                            dmc.Text(
                                "You'll receive a confirmation email shortly. Our curation team will review your submission.",
                                size="sm",
                            ),
                        ],
                        title="Submission Received",
                        color="green",
                        variant="light",
                    ),
                ]
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
                        dmc.Text("✓ Successfully submitted!", fw=600, size="lg"),
                        html.Br(),
                        dmc.Text(
                            f"Accession: {task_result['accession']}",
                            size="md",
                            fw=500,
                            c="green",
                        ),
                        dmc.Text(
                            f"Filename: {task_result['filename']}",
                            size="sm",
                            c="dimmed",
                        ),
                        html.Br(),
                        dmc.Text(
                            "Your submission will be reviewed by our curation team. You'll receive an email confirmation.",
                            size="sm",
                        ),
                    ],
                    title="Submission Complete",
                    color="green",
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
                                "user_message", "An error occurred during submission."
                            ),
                            size="sm",
                        ),
                    ],
                    title="Submission Failed",
                    color="red",
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
        return dmc.Text("✓ Valid email", size="xs", c="green")
    else:
        return dmc.Text("✗ Please enter a valid email address", size="xs", c="red")


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
            return html.Div(["There was an error processing this file."])

    # Show prefill details when no file has been manually uploaded
    if fasta_prefill:
        prefill_contents = fasta_prefill.get("contents")
        prefill_filename = fasta_prefill.get("filename", "sequence_from_blast.fa")
        if prefill_contents:
            try:
                _ct, content_string = prefill_contents.split(",", 1)
                query_string = base64.b64decode(content_string).decode("utf-8")
                records = parse_fasta(query_string, prefill_filename)
                return html.Div(
                    [
                        html.H6(f"File name: {prefill_filename}"),
                        html.H6(f"Number of sequences: {len(records)}"),
                        dmc.Text(
                            "Pre-filled from BLAST — upload a new file to override",
                            size="xs",
                            c="dimmed",
                            mt="xs",
                        ),
                    ]
                )
            except Exception as e:
                logger.error(e)

    try:
        content_type, content_string = seq_contents.split(",")
        query_string = base64.b64decode(content_string).decode("utf-8")
        records = parse_fasta(query_string, seq_filename)
        return create_fasta_display(records, seq_filename)
    except Exception as e:
        logger.error(e)
        return dmc.Alert(
            "There was an error processing this file. Please check the FASTA format.",
            title="Error",
            color="red",
            variant="light",
        )
    return html.Div(html.P(["Select a FASTA file to upload"]))


@callback(
    Output("submit-output-gff-upload", "children"),
    [
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
    ],
)
def update_gff_details(anno_contents, anno_filename):
    if anno_contents is None:

        except Exception as e:
            logger.error(e)
            return html.Div(["There was an error processing this file."])
    return html.Div(["Select a GFF file to upload"])


@callback(
    Output("status-display", "children"),
    Input("check-status-btn", "n_clicks"),
    State("status-submission-id", "value"),
    prevent_initial_call=True,
)
def check_submission_status(n_clicks, submission_id):
    """Check and display submission status."""
    if not n_clicks or not submission_id or not submission_id.strip():
        return html.Div("Please enter a submission ID", className="text-muted")

    submission_id = submission_id.strip()
    status_data = get_submission_status(submission_id)

    if not status_data:
        return dmc.Alert(
            "Submission not found or expired. Submission status is only available for 1 hour after submission.",
            title="Status Not Found",
            color="var(--mantine-color-orange-6)",
            variant="light",
        )

    status = status_data.get("status", "unknown")
    progress = status_data.get("progress", 0)
    message = status_data.get("message", "")
    created_at = status_data.get("created_at", "")
    updated_at = status_data.get("updated_at", "")

    # Status color mapping (design system: indigo for info, semantic for states)
    status_colors = {
        "queued": "indigo",
        "processing": "var(--mantine-color-yellow-6)",
        "completed": "var(--mantine-color-green-6)",
        "failed": "var(--mantine-color-red-6)",
        "unknown": "gray",
    }

    color = status_colors.get(status, "gray")

    status_display = [
        dmc.Group(
            [
                dmc.Text("Status:", fw=700),
                dmc.Badge(status.upper(), color=color, variant="light"),
            ]
        ),
        dmc.Progress(value=progress, color="var(--mantine-primary-color-6)", size="lg"),
        dmc.Text(f"Progress: {progress}%", size="sm", c="dimmed"),
        dmc.Text(message, size="sm"),
    ]

    if created_at:
        status_display.append(
            dmc.Text(f"Submitted: {created_at}", size="xs", c="dimmed")
        )

    if updated_at and updated_at != created_at:
        status_display.append(
            dmc.Text(f"Last updated: {updated_at}", size="xs", c="dimmed")
        )

    # Show result if completed
    result = status_data.get("result")
    if result and status == "completed" and result.get("success"):
        status_display.extend(
            [
                dmc.Divider(),
                dmc.Text("Result:", fw=700),
                dmc.Text(f"Accession: {result.get('accession', 'N/A')}", size="sm"),
                dmc.Text(
                    f"Review needed: {'Yes' if result.get('needs_review') else 'No'}",
                    size="sm",
                ),
                dmc.Text(f"Filename: {result.get('filename', 'N/A')}", size="sm"),
            ]
        )
    elif result and status == "failed":
        status_display.extend(
            [
                dmc.Divider(),
                dmc.Alert(
                    result.get("user_message", "Processing failed"),
                    title="Error",
                    color="var(--mantine-color-red-6)",
                    variant="light",
                ),
            ]
        )

    try:
        children = parse_gff(anno_contents, anno_filename)
        # Wrap in success alert if it's not already formatted
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
                title="GFF File Loaded Successfully",
                color="green",
                variant="light",
            )
        return children
    except Exception as e:
        logger.error(e)
        return dmc.Alert(
            "There was an error processing this file. Please check the format.",
            title="Error",
            color="red",
            variant="light",
        )
    return dmc.Stack(status_display, gap="xs")

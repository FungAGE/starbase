import base64
import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback

import datetime
from typing import Dict, Any, Optional
from src.utils.seq_utils import parse_fasta, parse_gff

from src.database.sql_engine import get_submissions_session
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy import text
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

logger = get_logger(__name__)


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
            "message": html.Div(
                [
                    html.H4("Validation Error", className="text-warning"),
                    html.P(str(error.user_message), className="text-muted"),
                ]
            ),
            "loading": False,
        }
    elif isinstance(error, ProcessingError):
        return {
            "modal_open": False,
            "message": html.Div(
                [
                    html.H4("Processing Error", className="text-danger"),
                    html.P(str(error.user_message), className="text-muted"),
                ]
            ),
            "loading": False,
        }
    elif isinstance(error, DatabaseError):
        return {
            "modal_open": False,
            "message": html.Div(
                [
                    html.H4("Database Error", className="text-danger"),
                    html.P(str(error.user_message), className="text-muted"),
                    html.P(
                        "If this persists, please contact support.",
                        className="text-muted small",
                    ),
                ]
            ),
            "loading": False,
        }
    else:
        return {
            "modal_open": False,
            "message": html.Div(
                [
                    html.H4("Unexpected Error", className="text-danger"),
                    html.P(
                        "An unexpected error occurred. Please try again.",
                        className="text-muted",
                    ),
                    html.P(
                        "If this persists, please contact support.",
                        className="text-muted small",
                    ),
                ]
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
    size="md",
    children=[
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title(
                    [
                        "Submit Starships to ",
                        html.Span("starbase", className="logo-text"),
                    ],
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    "Fields marked with * are required",
                    c="red",
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        # Form Content
        dmc.Paper(
            children=dbc.Form(
                [
                    dmc.Stack(
                        [
                            # File Upload Section
                            dmc.Stack(
                                [
                                    dmc.Title("Upload Files", order=2, mb="md"),
                                    # FASTA Upload
                                    dmc.Paper(
                                        children=[
                                            dmc.Text(
                                                [
                                                    "Starship Sequence ",
                                                    html.Span(
                                                        "*", style={"color": "red"}
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
                                    dmc.Title("Metadata", order=2, mb="md"),
                                    # Curator Info
                                    dmc.TextInput(
                                        id="uploader",
                                        label="Email of curator",
                                        placeholder="Enter email",
                                        required=True,
                                    ),
                                    dmc.TextInput(
                                        id="evidence",
                                        label="How were Starships annotated?",
                                        placeholder="i.e. starfish",
                                        required=True,
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
                                        placeholder="'chr1', or GenBank Accession",
                                        required=True,
                                    ),
                                    dmc.Group(
                                        [
                                            dmc.NumberInput(
                                                id="shipstart",
                                                label="Start coordinate",
                                                placeholder="1200",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                            dmc.NumberInput(
                                                id="shipend",
                                                label="End coordinate",
                                                placeholder="20500",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                        ]
                                    ),
                                    dmc.RadioGroup(
                                        id="strand-radios",
                                        label="Strand",
                                        value=1,
                                        children=[
                                            dmc.Radio(label="Positive strand", value=1),
                                            dmc.Space(h=10),
                                            dmc.Radio(label="Negative strand", value=2),
                                        ],
                                    ),
                                    # Comments
                                    dmc.Textarea(
                                        id="comment",
                                        label="Additional information",
                                        placeholder="Any comments about the Starship features, annotations, or host genome?",
                                        minRows=3,
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
                                    variant="gradient",
                                    gradient={"from": "indigo", "to": "cyan"},
                                    loading=False,
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
        # Submission Status Section
        dmc.Paper(
            children=[
                dmc.Title("Submission Status", order=3, mb="md"),
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
                    dbc.ModalTitle("Submission Status", className="text-info"),
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
        "padding": "2rem",
    },
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
        content_type, content_string = seq_contents.split(",")
        seq_decoded = base64.b64decode(content_string).decode("utf-8")
        seq_datetime_obj = datetime.datetime.fromtimestamp(seq_date).strftime(
            "%Y-%m-%d %H:%M:%S"
        )

        anno_contents = ""
        anno_datetime_obj = ""

        if anno_contents:
            content_type, content_string = anno_contents.split(",")
            anno_contents = base64.b64decode(content_string).decode("utf-8")
            anno_datetime_obj = datetime.datetime.fromtimestamp(anno_date).strftime(
                "%Y-%m-%d %H:%M:%S"
            )

        with get_submissions_session() as session:
            query = """
                INSERT INTO submissions (
                    seq_contents, seq_filename, seq_date, anno_contents,
                    anno_filename, anno_date, uploader, evidence,
                    genus, species, hostchr, shipstart, shipend,
                    shipstrand, comment, accession_tag, needs_review
                )
                VALUES (
                    :seq_contents, :seq_filename, :seq_date, :anno_contents, 
                    :anno_filename, :anno_date, :uploader, :evidence, 
                    :genus, :species, :hostchr, :shipstart, :shipend, 
                    :shipstrand, :comment, :accession_tag, :needs_review
                )
            """
            session.execute(
                text(query),
                {
                    "seq_contents": seq_decoded,
                    "seq_filename": seq_filename,
                    "seq_date": seq_datetime_obj,
                    "anno_contents": anno_contents,
                    "anno_filename": anno_filename,
                    "anno_date": anno_datetime_obj,
                    "uploader": uploader,
                    "evidence": evidence,
                    "genus": genus,
                    "species": species,
                    "hostchr": hostchr,
                    "shipstart": shipstart,
                    "shipend": shipend,
                    "shipstrand": shipstrand,
                    "comment": comment,
                    "accession_tag": accession,
                    "needs_review": needs_review,
                },
            )
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
    return html.Div(
        [
            html.H6(f"File name: {filename}"),
            html.H6(f"Number of sequences: {len(records)}"),
        ]
    )


@callback(
    [
        Output("submit-modal", "is_open"),
        Output("output-data-upload", "children"),
        Output("submit-ship", "loading"),
    ],
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-upload", "last_modified"),
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
        Input("submit-upload-gff", "last_modified"),
        Input("uploader", "value"),
        Input("evidence", "value"),
        Input("genus", "value"),
        Input("species", "value"),
        Input("hostchr", "value"),
        Input("shipstart", "value"),
        Input("shipend", "value"),
        Input("strand-radios", "value"),
        Input("comment", "value"),
        Input("submit-ship", "n_clicks"),
        Input("close", "n_clicks"),
    ],
    [State("submit-modal", "is_open")],
)
def submit_ship(
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
    n_clicks,
    close_modal,
    is_open,
):
    """
    Main submission callback with asynchronous processing.

    This function validates input, starts an async task, and provides immediate feedback.
    Users can check progress via submission status.
    """
    modal = is_open
    message = ""
    loading = False

    # Only process if submit button was clicked
    if not n_clicks or n_clicks <= 0:
        return modal, message, loading

    loading = True

    try:
        # Step 1: Validate input data (quick validation)
        # This must pass before we proceed with submission
        # If validation fails, validate_submission_data will raise ValidationError
        logger.debug("Validating submission data")
        _ = validate_submission_data(
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

        # Step 2: Create submission data package
        submission_data = {
            "seq_contents": seq_contents,
            "seq_filename": seq_filename,
            "seq_date": seq_date,
            "anno_contents": anno_contents,
            "anno_filename": anno_filename,
            "anno_date": anno_date,
            "uploader": uploader,
            "evidence": evidence,
            "genus": genus,
            "species": species,
            "hostchr": hostchr,
            "shipstart": shipstart,
            "shipend": shipend,
            "strand_radio": strand_radio,
            "comment": comment or "",
        }

        # Step 3: Generate unique submission ID
        submission_id = str(uuid.uuid4())
        logger.info(f"Created submission {submission_id} for file {seq_filename}")

        # Step 4: Create initial status
        create_submission_status(submission_id, "queued")

        # Step 5: Start async processing task
        logger.debug(f"Starting async submission processing for {submission_id}")
        task_result = run_task(process_submission_task, submission_data, submission_id)

        # Store task ID for status checking
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

        # Step 6: Show immediate feedback with submission ID
        modal = not is_open if not close_modal else False
        if hasattr(task_result, "id"):
            # Async processing
            message = html.Div(
                [
                    html.H4("Submission Queued!", className="mb-3"),
                    dmc.Text(
                        "Your submission is being processed in the background.",
                        className="text-muted",
                    ),
                    dmc.Text(f"Submission ID: {submission_id}", className="text-muted"),
                    html.Br(),
                    dmc.Text(
                        "You can close this dialog. Processing will continue.",
                        className="text-muted small",
                    ),
                    html.Br(),
                    dmc.Text(
                        "Check the status below to see when processing is complete.",
                        className="text-muted small",
                    ),
                ]
            )
        else:
            # Sync processing completed
            if task_result.get("success"):
                message = html.Div(
                    [
                        html.H4("Successfully submitted!", className="mb-3"),
                        dmc.Text(
                            f"Assigned accession: {task_result['accession']}",
                            className="text-muted",
                        ),
                        dmc.Text(
                            f"Review status: {'Needs review' if task_result['needs_review'] else 'Auto-approved'}",
                            className="text-muted",
                        ),
                        dmc.Text(
                            f"Filename: {task_result['filename']}",
                            className="text-muted",
                        ),
                        dmc.Text(
                            f"Uploaded by: {task_result['uploader']}",
                            className="text-muted",
                        ),
                    ]
                )
            else:
                message = html.Div(
                    [
                        html.H4("Submission Failed", className="text-danger mb-3"),
                        dmc.Text(
                            task_result.get("user_message", "An error occurred"),
                            className="text-muted",
                        ),
                    ]
                )

        loading = False
        return modal, message, loading

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
        )

    except Exception as e:
        # Handle unexpected errors
        logger.error(f"Unexpected error in submit_ship: {str(e)}", exc_info=True)
        error_response = handle_submission_error(e)
        return (
            error_response["modal_open"],
            error_response["message"],
            error_response["loading"],
        )


@callback(
    Output("submit-fasta-sequence-upload", "children"),
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
    ],
)
def update_fasta_details(seq_contents, seq_filename):
    if seq_contents is None:
        return html.Div(html.P(["Select a FASTA file to upload"]))

    try:
        content_type, content_string = seq_contents.split(",")
        query_string = base64.b64decode(content_string).decode("utf-8")
        records = parse_fasta(query_string, seq_filename)
        return create_fasta_display(records, seq_filename)
    except Exception as e:
        logger.error(e)
        return html.Div(["There was an error processing this file."])


@callback(
    Output("submit-output-gff-upload", "children"),
    [
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
    ],
)
def update_gff_details(anno_contents, anno_filename):
    if anno_contents is None:
        return [
            html.Div(["Select a GFF file to upload"]),
        ]
    else:
        try:
            children = parse_gff(anno_contents, anno_filename)
            return children

        except Exception as e:
            logger.error(e)
            return html.Div(["There was an error processing this file."])


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
            color="orange",
            variant="light",
        )

    status = status_data.get("status", "unknown")
    progress = status_data.get("progress", 0)
    message = status_data.get("message", "")
    created_at = status_data.get("created_at", "")
    updated_at = status_data.get("updated_at", "")

    # Status color mapping
    status_colors = {
        "queued": "blue",
        "processing": "yellow",
        "completed": "green",
        "failed": "red",
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
        dmc.Progress(value=progress, color=color, size="lg"),
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
                    color="red",
                    variant="light",
                ),
            ]
        )

    return dmc.Stack(status_display, gap="xs")

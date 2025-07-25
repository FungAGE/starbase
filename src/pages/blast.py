import dash
from dash import dcc, html, callback, clientside_callback
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import os
import base64
import pandas as pd
import time
import tempfile

from src.config.cache import cache
from src.utils.seq_utils import (
    check_input,
    write_temp_fasta,
    seq_processing_error_alert,
)
from src.utils.blast_utils import (
    create_no_matches_alert,
    parse_blast_xml,
    select_ship_family,
)

from src.components.callbacks import (
    curated_switch,
    create_file_upload,
)
from src.database.sql_manager import fetch_meta_data
from src.utils.classification_utils import WORKFLOW_STAGES
from src.tasks import (
    run_blast_search_task,
    run_hmmer_search_task,
    run_classification_workflow_task,
)

from src.config.logging import get_logger

from src.utils.blast_data import BlastData, ClassificationData

dash.register_page(__name__)

logger = get_logger(__name__)

progress_section = dmc.Stack(
    [
        dmc.Group(
            [
                dbc.Progress(
                    id="classification-progress",
                    value=0,
                    color="blue",
                    animated=True,
                    striped=True,
                    style={
                        "width": "100%",
                        "marginBottom": "5px",
                    },
                ),
            ]
        ),
        dmc.Group(
            [
                dmc.Text(
                    "Classification Status:",
                    size="lg",
                    fw=500,
                ),
                dmc.Text(
                    id="classification-stage-display",
                    size="lg",
                    c="blue",
                ),
            ]
        ),
    ],
    gap="md",
    id="classification-progress-section",
    style={"display": "none"},
)

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        # Upload stores
        dcc.Store(id="blast-sequences-store"),
        dcc.Store(id="upload-error-store"),
        # Submission id
        dcc.Store(id="submission-id-store"),
        # Processed data stores
        dcc.Store(id="processed-metadata-store"),
        dcc.Store(id="processed-blast-store"),
        dcc.Store(
            id="blast-processed-tabs", data=[0]
        ),  # Store which tabs have been processed
        # Results stores
        dcc.Store(id="blast-results-store"),
        dcc.Store(id="classification-result-store", data=None),
        # Single data store with all workflow state
        dcc.Store(id="classification-stage", data="Upload a sequence"),
        dcc.Store(
            id="classification-workflow-state",
            data={"current_stage": None, "complete": False},
        ),
        # Store for each workflow stage
        *[
            dcc.Store(id={"type": "classification-stage-data", "index": stage["id"]})
            for stage in WORKFLOW_STAGES
        ],
        # Tab stores
        dcc.Store(id="blast-active-tab", data=0),  # Store active tab index
        # Interval for polling workflow state is no longer needed
        # but kept disabled for backward compatibility
        dcc.Interval(
            id="classification-interval",
            interval=500,
            disabled=True,  # Always disabled
            max_intervals=0,  # Never trigger
        ),
        dcc.Store(id="workflow-state-store", data=None),
        dmc.Space(h=20),
        dmc.Grid(
            children=[
                dmc.GridCol(
                    span={"sm": 12, "lg": 4},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack(
                                [
                                    dmc.Title("BLAST Search", order=1),
                                    dmc.Text(
                                        "Search protein/nucleotide sequences for Starships and Starship-associated genes",
                                        c="dimmed",
                                        size="lg",
                                    ),
                                    # Input Section
                                    dmc.Stack(
                                        [
                                            dmc.Title("Input Sequence(s)", order=3),
                                            dmc.Textarea(
                                                id="query-text",
                                                placeholder="Paste a maximum of 10 FASTA sequences here...",
                                                autosize=True,
                                                minRows=5,
                                                maxRows=15,
                                                style={"width": "100%"},
                                            ),
                                            # Upload Section
                                            dmc.Stack(
                                                [
                                                    dmc.Center(
                                                        dmc.Text("Or", size="lg"),
                                                    ),
                                                    dmc.Paper(
                                                        children=create_file_upload(
                                                            upload_id="blast-fasta-upload",
                                                            output_id="upload-details",
                                                            accept_types=[
                                                                ".fa",
                                                                ".fas",
                                                                ".fasta",
                                                                ".fna",
                                                            ],
                                                            placeholder_text=dmc.Stack(
                                                                [
                                                                    dmc.Text(
                                                                        "Drag and drop or click to select a FASTA file",
                                                                        size="sm",
                                                                    ),
                                                                    dmc.Text(
                                                                        "Maximum of 10 sequences",
                                                                        size="sm",
                                                                        c="dimmed",
                                                                    ),
                                                                ],
                                                                gap="xs",
                                                            ),
                                                        ),
                                                        withBorder=False,
                                                        shadow="sm",
                                                        radius="md",
                                                        style={"cursor": "pointer"},
                                                    ),
                                                    html.Div(
                                                        id="upload-error-message",
                                                        style={"color": "red"},
                                                    ),
                                                ],
                                                gap="md",
                                            ),
                                            # Options Section
                                            dmc.Stack(
                                                [
                                                    dmc.Title(
                                                        "Search Options", order=3
                                                    ),
                                                    dmc.Grid(
                                                        children=[
                                                            dmc.GridCol(
                                                                span={
                                                                    "xs": 12,
                                                                    "sm": 6,
                                                                },
                                                                children=[
                                                                    dmc.Stack(
                                                                        [
                                                                            curated_switch(
                                                                                text="Only search curated Starships",
                                                                                size="sm",
                                                                            ),
                                                                        ],
                                                                        gap="xs",
                                                                    ),
                                                                ],
                                                            ),
                                                            dmc.GridCol(
                                                                span={
                                                                    "xs": 12,
                                                                    "sm": 6,
                                                                },
                                                                children=[
                                                                    dmc.NumberInput(
                                                                        id="evalue-threshold",
                                                                        label="E-value Threshold",
                                                                        value=0.001,
                                                                        min=0,
                                                                        max=1,
                                                                        step=0.005,
                                                                    ),
                                                                ],
                                                            ),
                                                        ],
                                                        gutter="md",
                                                        align="center",
                                                    ),
                                                ],
                                                gap="xs",
                                            ),
                                        ],
                                        gap="md",
                                    ),
                                    dmc.Space(h="md"),
                                    # Submit Section
                                    dmc.Stack(
                                        [
                                            dmc.Center(
                                                dmc.Button(
                                                    "Submit BLAST",
                                                    id="submit-button",
                                                    variant="gradient",
                                                    gradient={
                                                        "from": "indigo",
                                                        "to": "cyan",
                                                    },
                                                    size="lg",
                                                    loading=False,
                                                    loaderProps={
                                                        "variant": "dots",
                                                        "color": "white",
                                                    },
                                                ),
                                            ),
                                        ],
                                        gap="xs",
                                    ),
                                ],
                                gap="md",
                            ),
                            p="xl",
                            radius="md",
                            shadow="sm",
                            withBorder=True,
                            # style={"height": "100%"},
                        ),
                    ],
                ),
                # Right Column - Results Panel
                dmc.GridCol(
                    span={"sm": 12, "lg": 8},
                    children=[
                        dmc.Stack(
                            pos="relative",
                            children=[
                                dmc.LoadingOverlay(
                                    visible=False,
                                    id="results-loading-overlay",
                                    overlayProps={"radius": "sm", "blur": 2},
                                    loaderProps={
                                        "variant": "oval",
                                        "size": "xl",
                                        "color": "blue",
                                    },
                                    zIndex=10,
                                ),
                                html.Div(
                                    id="right-column-content",
                                    children=[
                                        # This div will be replaced with tabs when more than one sequence is in query
                                        dmc.Stack(
                                            children=[
                                                html.Div(
                                                    id="classification-output",
                                                    className="mt-4",
                                                ),
                                                progress_section,
                                                # BLAST results section
                                                dmc.Stack(
                                                    [
                                                        html.Div(
                                                            id="blast-container",
                                                            className="blast-container",
                                                            style={
                                                                "width": "100%",
                                                                "display": "flex",
                                                                "flexDirection": "column",
                                                                "minHeight": "300px",
                                                                "alignItems": "flex-start",
                                                                "textAlign": "left",
                                                            },
                                                        ),
                                                    ]
                                                ),
                                            ],
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                    style={"justify-content": "flex-start"},
                ),
            ],
            gutter="xl",
        ),
    ],
)

########################################################
# Handle Sequence Input
########################################################


@callback(
    [
        Output("submit-button", "disabled", allow_duplicate=True),
        Output("submit-button", "children", allow_duplicate=True),
        Output("blast-sequences-store", "data", allow_duplicate=True),
        Output("upload-details", "children", allow_duplicate=True),
        Output("upload-error-message", "children", allow_duplicate=True),
    ],
    [
        Input("blast-fasta-upload", "contents"),
        Input("blast-fasta-upload", "filename"),
    ],
    prevent_initial_call=True,
)
def update_fasta_details(seq_content, seq_filename):
    """
    Handle file uploads only - immediately process the file to show a summary
    but don't display errors in the main interface - just in the upload area.
    """
    button_disabled = False
    button_text = "Submit BLAST"
    seq_list = None
    error_alert = None

    # Default upload details - maintain this if there's an error
    upload_details = html.Div(
        html.P(
            ["Select a FASTA file to upload"],
        )
    )

    if seq_content is None:
        return button_disabled, button_text, seq_list, upload_details, error_alert

    try:
        # Check for file size and prevent upload if too large
        try:
            max_size = 10 * 1024 * 1024  # 10 MB
            content_type, content_string = seq_content.split(",")
            decoded = base64.b64decode(content_string)
            file_size = len(decoded)
        except Exception as e:
            logger.error(f"Error decoding file contents: {e}")
            error_alert = dmc.Alert(
                title="Invalid File Format",
                children="The file appears to be corrupted or in an unsupported format.",
                color="red",
                variant="filled",
            )
            return True, "Error", None, upload_details, error_alert

        # Check file size constraint
        if file_size > max_size:
            error_msg = f"The file '{seq_filename}' exceeds the 10 MB limit."
            error_alert = dmc.Alert(
                title="File Too Large",
                children=error_msg,
                color="red",
                variant="filled",
            )
            return True, "Error", None, upload_details, error_alert

        # Parse the file contents using check_input
        input_type, seq_list, n_seqs, error = check_input(
            query_text_input=None, query_file_contents=seq_content
        )

        # Save success/error status in a log
        if seq_list:
            logger.debug(
                f"Successfully parsed file upload: {seq_filename}, type={input_type}, sequences={n_seqs}"
            )
        else:
            logger.warning(
                f"Failed to parse file upload: {seq_filename}, error={error}"
            )

        # Handle parsing errors
        if error:
            error_alert = (
                error
                if not isinstance(error, str)
                else dmc.Alert(
                    title="Error Processing File",
                    children=error,
                    color="red",
                    variant="filled",
                )
            )
            return True, "Error", None, upload_details, error_alert

        # Check if seq_list is None
        if seq_list is None:
            error_alert = dmc.Alert(
                title="Parsing Error",
                children="Failed to parse sequences from the file.",
                color="red",
                variant="filled",
            )
            return True, "Error", None, upload_details, error_alert

        # Process valid sequences - create a summary for display
        if n_seqs > 1:
            upload_details = html.Div(
                [
                    html.P(f"File name: {seq_filename}"),
                    html.P(
                        [
                            f"Found {n_seqs} sequences in the file",
                            html.Span(
                                " (maximum of 10 will be processed)",
                                style={"color": "orange", "fontStyle": "italic"}
                                if n_seqs > 10
                                else {"display": "none"},
                            ),
                        ]
                    ),
                    html.Ul(
                        [
                            html.Li(f"{seq['header']} ({len(seq['sequence'])} bp)")
                            for seq in seq_list[:3]
                        ]
                        + ([html.Li("...")] if n_seqs > 3 else [])
                    ),
                    html.P(
                        "Only the first 10 sequences will be processed."
                        if n_seqs > 10
                        else "All sequences will be processed when you submit."
                    ),
                ]
            )
        else:
            upload_details = html.Div(
                [
                    html.H6(f"File name: {seq_filename}"),
                    html.H6(f"Number of sequences: {n_seqs}"),
                    html.P(
                        f"Sequence: {seq_list[0]['header']} ({len(seq_list[0]['sequence'])} bp)",
                        style={"fontSize": "14px", "color": "#666"},
                    ),
                ]
            )

        # Store sequence information in an additional data attribute to help with debugging
        upload_details.n_sequences = n_seqs  # Add custom attribute

        return button_disabled, button_text, seq_list, upload_details, error_alert

    except Exception as e:
        logger.error(f"Error processing file: {e}")
        error_alert = dmc.Alert(
            title="Error Processing File",
            children=f"An unexpected error occurred: {str(e)}",
            color="red",
            variant="filled",
        )
        return True, "Error", None, upload_details, error_alert


@callback(
    Output("query-text", "value", allow_duplicate=True),
    Input("blast-fasta-upload", "contents"),
    prevent_initial_call=True,
    id="clear-text-on-file-upload",
)
def clear_text_on_file_upload(file_contents):
    """Clear the text area when a file is uploaded"""
    if file_contents:
        return ""
    raise PreventUpdate


@callback(
    [
        Output("blast-fasta-upload", "contents", allow_duplicate=True),
        Output("upload-details", "children", allow_duplicate=True),
    ],
    Input("query-text", "value"),
    State("blast-fasta-upload", "contents"),
    prevent_initial_call=True,
    id="clear-file-on-text-input",
)
def clear_file_on_text_input(text_value, current_file_contents):
    """Clear the file upload when text is entered"""
    if text_value and len(text_value.strip()) > 10 and current_file_contents:
        return None, html.Div(html.P(["Select a FASTA file to upload"]))
    raise PreventUpdate


@callback(
    [
        Output("blast-sequences-store", "data", allow_duplicate=True),
        Output("query-text", "helperText"),
        Output("upload-error-message", "children", allow_duplicate=True),
    ],
    [
        Input("query-text", "value"),
    ],
    prevent_initial_call=True,
)
def preprocess_text_input(query_text_input):
    """
    Process text input as it's entered to show a summary and limit sequences.
    This doesn't store the sequences for processing yet - that happens when submit is clicked.
    """
    if not query_text_input:
        return None, None, None

    try:
        # Parse the text input to check for multiple sequences
        input_type, seq_list, n_seqs, error = check_input(
            query_text_input=query_text_input,
            query_file_contents=None,
            max_sequences=10,
        )

        # If there's a real error, show it but don't update sequences store
        if error and isinstance(error, dmc.Alert) and error.color == "red":
            return None, None, error

        # If we have sequences, create a summary message
        if seq_list and n_seqs > 0:
            if n_seqs == 1:
                helper_text = f"Found 1 sequence: {seq_list[0]['header']} ({len(seq_list[0]['sequence'])} bp)"
            else:
                helper_text = f"Found {n_seqs} sequences"
                if n_seqs > 3:
                    seq_summary = (
                        ", ".join([seq["header"] for seq in seq_list[:3]]) + "..."
                    )
                else:
                    seq_summary = ", ".join([seq["header"] for seq in seq_list])
                helper_text += f" ({seq_summary})"

            # Add warning if we limited the sequences
            if n_seqs > 10:
                helper_text += " (only the first 10 will be processed)"

            # Don't store sequences yet - just return helper text
            return None, helper_text, None

        return None, "No valid sequences found", None

    except Exception as e:
        logger.error(f"Error preprocessing text input: {e}")
        return (
            None,
            None,
            dmc.Alert(
                title="Error",
                children=f"Error preprocessing sequences: {str(e)}",
                color="red",
                variant="filled",
            ),
        )


# Submission ID is now generated by the preprocess callback to ensure proper ordering


# Add the clientside callback to immediately set button loading state
clientside_callback(
    """
    function updateLoadingState(n_clicks) {
        return true;
    }
    """,
    Output("submit-button", "loading", allow_duplicate=True),
    Input("submit-button", "n_clicks"),
    prevent_initial_call=True,
)

# Add a clientside callback to show the loading overlay
clientside_callback(
    """
    function showLoadingOverlay(n_clicks) {
        return true;
    }
    """,
    Output("results-loading-overlay", "visible", allow_duplicate=True),
    Input("submit-button", "n_clicks"),
    prevent_initial_call=True,
)


########################################################
# Setup Layout
########################################################
# Helper functions to create layouts
def create_single_layout():
    return dmc.Stack(
        children=[
            html.Div(id="classification-output", className="mt-4"),
            # BLAST results section - use create_blast_container for consistency
            dmc.Stack(
                [
                    dbc.Spinner(
                        children=html.Div(
                            id="blast-container",
                            className="blast-container",
                            style={
                                "width": "100%",
                                "display": "flex",
                                "flexDirection": "column",
                                "minHeight": "300px",
                                "alignItems": "flex-start",
                                "textAlign": "left",
                            },
                        ),
                        color="primary",
                        type="border",
                        fullscreen=False,
                        spinner_style={"width": "3rem", "height": "3rem"},
                    )
                ]
            ),
        ],
    )


def create_tab_layout(seq_list, n_seqs=10):
    """Create tabs for multiple sequences"""
    tabs = []
    for idx, seq in enumerate(seq_list[:n_seqs]):
        tabs.append(
            dbc.Tab(
                label=f"Sequence {idx + 1}: {seq['header'][:20]}{'...' if len(seq['header']) > 20 else ''}",
                tab_id=f"tab-{idx}",
            )
        )

    return dbc.Card(
        [
            dbc.CardHeader(
                dbc.Tabs(id="blast-tabs", active_tab="tab-0", children=tabs)
            ),
            dbc.CardBody([html.Div(id="tab-content")]),
        ]
    )


def create_blast_container(sequence_results, tab_id=None):
    """Create the BLAST container with a spinner wrapper

    Args:
        sequence_results: The BLAST results data
        tab_id: Optional ID suffix for creating unique container IDs in a tabbed interface
    """
    blast_content = sequence_results.get("blast_content")
    if not blast_content:
        return html.Div(
            [
                html.H2(
                    "BLAST Results", style={"marginTop": "15px", "marginBottom": "20px"}
                ),
                create_no_matches_alert(),
            ]
        )

    # Create a unique container ID if tab_id is provided
    container_id = (
        f"blast-container-{tab_id}" if tab_id is not None else "blast-container"
    )

    # For single sequences, ensure we're using the standard ID
    if tab_id is None:
        logger.debug(
            "Creating blast container for single sequence with ID: blast-container"
        )
    else:
        logger.debug(
            f"Creating blast container for tab {tab_id} with ID: {container_id}"
        )

    # Create with a regular string ID and wrap in a spinner
    return html.Div(
        [
            dbc.Spinner(
                children=html.Div(
                    id=container_id,
                    className="blast-container",  # Add a common class for all containers
                    style={
                        "width": "100%",
                        "display": "flex",
                        "flexDirection": "column",
                        "minHeight": "300px",  # Increased height to make spinner more visible
                        "alignItems": "flex-start",  # Left align the content
                        "textAlign": "left",  # Ensure text is left-aligned
                    },
                ),
                color="primary",
                type="border",  # Options: "border", "grow"
                fullscreen=False,
                spinner_style={"width": "3rem", "height": "3rem"},
            )
        ]
    )


########################################################
# Processing Data
########################################################
@callback(
    Output("processed-metadata-store", "data"),
    Input("curated-input", "value"),
)
def process_metadata(curated):
    try:
        initial_df = cache.get("meta_data")
        if initial_df is None:
            initial_df = fetch_meta_data(curated=curated)

        return (
            initial_df[["accession_tag", "familyName"]]
            .drop_duplicates()
            .to_dict("records")
        )
    except Exception as e:
        logger.error(f"Error processing metadata: {str(e)}")
        return None


@callback(
    Output("processed-blast-store", "data"),
    [Input("blast-results-store", "data"), Input("blast-active-tab", "data")],
)
def process_blast_results(blast_results_store, active_tab_idx):
    if not blast_results_store:
        logger.warning("No blast_results_store data")
        return None

    # Convert active_tab_idx to string (since keys in sequence_results are strings)
    tab_idx = str(active_tab_idx or 0)
    logger.debug(f"Processing BLAST results for tab index: {tab_idx}")

    # Get the correct sequence results based on active tab
    sequence_results = blast_results_store.get("sequence_results", {}).get(tab_idx)
    if not sequence_results:
        logger.warning(f"No sequence results for tab index {tab_idx}")
        return None

    # Check for blast_file or blast_content directly in sequence_results
    blast_results_file = sequence_results.get("blast_file")
    blast_content = sequence_results.get("blast_content")

    # If we have direct blast_content, use it without reading file
    if blast_content:
        logger.debug(f"Using direct blast_content for tab {tab_idx}")
        return {"blast_text": blast_content}

    # If no blast_file, return empty blast text with warning
    if not blast_results_file:
        logger.warning(f"No blast file in sequence results for tab {tab_idx}")
        return {"blast_text": ""}

    logger.debug(f"Reading BLAST file: {blast_results_file}")
    try:
        # Check if file exists
        if not os.path.exists(blast_results_file):
            logger.error(f"BLAST file does not exist: {blast_results_file}")
            return {"blast_text": ""}

        # Read the BLAST output as text
        with open(blast_results_file, "r") as f:
            blast_results = f.read()

        # Add size limit check
        results_size = len(blast_results)
        logger.debug(f"Read BLAST results, size: {results_size} bytes")
        if results_size > 5 * 1024 * 1024:  # 5MB limit
            logger.warning(f"BLAST results too large: {results_size} bytes")
            return {"blast_text": "BLAST results too large to display"}

        # Format data for BlasterJS
        data = {
            "blast_text": blast_results  # Pass the raw BLAST text directly
        }

        return data
    except Exception as e:
        logger.error(f"Error processing BLAST results: {str(e)}")
        # Always return consistent data format, even on error
        return {"blast_text": ""}  # Return empty string instead of None


########################################################
# Main Processing
########################################################


@callback(
    [
        Output("blast-sequences-store", "data", allow_duplicate=True),
        Output("upload-error-store", "data", allow_duplicate=True),
        Output("upload-error-message", "children", allow_duplicate=True),
        Output("right-column-content", "children", allow_duplicate=True),
        Output("submission-id-store", "data", allow_duplicate=True),
    ],
    [
        Input("submit-button", "n_clicks"),
    ],
    [
        State("query-text", "value"),
        State("blast-sequences-store", "data"),
        State(
            "blast-fasta-upload", "contents"
        ),  # Add file contents as state to handle race conditions
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
    id="preprocess-callback",
)
def preprocess(n_clicks, query_text_input, seq_list, file_contents):
    """
    Process input when the submit button is pressed.
    For file uploads: Use the already parsed sequences from blast-sequences-store
    For text input: Parse the input text now, limiting to max 10 sequences
    """
    if not n_clicks:
        raise PreventUpdate

    error_alert = None
    ui_content = None  # Initialize UI content

    try:
        # Log what we're working with to help with debugging
        logger.debug(
            f"Preprocess called with seq_list={seq_list is not None}, query_text_input={query_text_input is not None}, file_contents={file_contents is not None}"
        )

        # If we have parsed sequences from a file upload, use those
        if seq_list is not None:
            logger.debug(
                f"Using pre-parsed sequences from blast-sequences-store: {len(seq_list)} sequences"
            )
            # Ensure limit to 10 sequences
            if len(seq_list) > 10:
                seq_list = seq_list[:10]

            # Create UI structure for these sequences
            ui_content = (
                create_tab_layout(seq_list)
                if len(seq_list) > 1
                else create_single_layout()
            )

            return seq_list, None, None, ui_content, n_clicks

        # If we have file contents but no parsed sequences, try to parse them now
        # This handles cases where the file upload callback hasn't completed before submit
        if file_contents and not seq_list:
            logger.debug("Parsing file contents directly from state")
            try:
                input_type, seq_list, n_seqs, error = check_input(
                    query_text_input=None, query_file_contents=file_contents
                )

                if error:
                    logger.error(f"Error parsing file contents: {error}")
                    error_alert = (
                        error
                        if not isinstance(error, str)
                        else dmc.Alert(
                            title="Error Processing File",
                            children=error,
                            color="red",
                            variant="filled",
                        )
                    )
                    return (
                        None,
                        str(error) if isinstance(error, str) else "Parsing error",
                        error_alert,
                        None,
                        n_clicks,
                    )

                if seq_list and len(seq_list) > 0:
                    logger.debug(
                        f"Successfully parsed {len(seq_list)} sequences from file contents"
                    )
                    # Ensure limit to 10 sequences
                    if len(seq_list) > 10:
                        seq_list = seq_list[:10]

                    # Create UI structure for these sequences
                    ui_content = (
                        create_tab_layout(seq_list)
                        if len(seq_list) > 1
                        else create_single_layout()
                    )

                    return seq_list, None, None, ui_content, n_clicks
            except Exception as e:
                logger.error(f"Error processing file contents directly: {e}")
                # Fall through to text input handling

        # Otherwise, parse the text input with max_sequences=10
        if query_text_input:
            logger.debug("Processing text input")
            input_type, seq_list, n_seqs, error = check_input(
                query_text_input=query_text_input,
                query_file_contents=None,
                max_sequences=10,
            )

            if error:
                logger.error(f"Error parsing text input: {error}")
                error_alert = (
                    error
                    if not isinstance(error, str)
                    else seq_processing_error_alert(error)
                )
                return (
                    None,
                    str(error) if isinstance(error, str) else "Parsing error",
                    error_alert,
                    None,
                    n_clicks,
                )

            logger.debug(f"Text input processed: {input_type}, {n_seqs} sequences")

            # Create UI structure based on number of sequences from text input
            ui_content = (
                create_tab_layout(seq_list)
                if len(seq_list) > 1
                else create_single_layout()
            )

            return seq_list, None, None, ui_content, n_clicks

        # If we have neither text input nor uploaded sequences, show an error
        error_alert = dmc.Alert(
            title="No Input Provided",
            children="Please enter a sequence or upload a FASTA file.",
            color="yellow",
            variant="filled",
        )
        return None, "No input provided", error_alert, None, n_clicks

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        error_alert = dmc.Alert(
            title="Error Processing Input",
            children=f"An unexpected error occurred: {str(e)}",
            color="red",
            variant="filled",
        )
        return None, str(e), error_alert, None, n_clicks


def process_single_sequence(seq_data, evalue_threshold):
    """Process a single sequence and return structured results"""

    # extract from seq_data
    query_header = seq_data.get("header", "query")
    query_seq = seq_data.get("sequence", "")
    query_type = seq_data.get("type", "nucl")
    query_length = len(query_seq)

    class_dict = BlastData(
        seq_type=query_type,
        sequence=query_seq,
    )

    classification_data = ClassificationData()

    # Write sequence to temporary FASTA file
    tmp_query_fasta = write_temp_fasta(query_header, query_seq)
    if not tmp_query_fasta:
        error_message = "Failed to create temporary file"
        logger.error(error_message)
        class_dict.error = error_message

    # fetch metadata
    meta_df = fetch_meta_data()

    try:
        # Run BLAST search
        blast_results = run_blast_search_task(
            query_header=query_header,
            query_seq=query_seq,
            query_type=query_type,
            eval_threshold=evalue_threshold,
        )

        # Handle case where blast_results is None or invalid
        if not blast_results:
            logger.warning("BLAST search returned no results")
        elif isinstance(blast_results, dict) and "content" in blast_results:
            # Handle the new format (content-based instead of file path)
            blast_content = blast_results["content"]
            blast_results_file = tempfile.NamedTemporaryFile(
                suffix=".blast", delete=False
            ).name
            with open(blast_results_file, "w") as f:
                f.write(blast_content)
            logger.debug(f"Created temporary BLAST results file: {blast_results_file}")
        elif isinstance(blast_results, str) and os.path.exists(blast_results):
            # Handle the old format (file path)
            blast_results_file = blast_results
            # Also read content from file for direct access
            try:
                with open(blast_results_file, "r") as f:
                    blast_content = f.read()
            except Exception as e:
                logger.error(f"Error reading blast file content: {e}")
            logger.debug(f"Using existing BLAST results file: {blast_results_file}")
        else:
            # Handle unexpected format
            logger.error(f"Unexpected BLAST results format: {type(blast_results)}")
            if isinstance(blast_results, str):
                # If it's a string but not a file, it might be the actual content
                blast_content = blast_results
                blast_results_file = None
                logger.debug("Using blast results string as content")
            else:
                blast_results_file = None
                blast_content = None

        # Update class_dict with blast results
        class_dict.blast_file = blast_results_file
        class_dict.blast_content = blast_content

        # Process either with file or direct content
        if blast_results_file and os.path.exists(blast_results_file):
            blast_tsv = parse_blast_xml(blast_results_file)
            if blast_tsv and os.path.exists(blast_tsv):
                blast_df = pd.read_csv(blast_tsv, sep="\t")
                if len(blast_df) == 0:
                    error_message = "No BLAST hits found"
                    logger.debug(error_message)

        # append to class_dict
        class_dict.fasta_file = tmp_query_fasta
        class_dict.blast_df = (
            blast_df.to_dict("records")
            if isinstance(blast_df, pd.DataFrame)
            else blast_df
        )

        try:
            workflow_result = None

            # Check for classification result data first
            if workflow_result and workflow_result.get("complete", False):
                match_stage = workflow_result.get("match_stage")
                match_result = workflow_result.get("match_result")
                found_match = workflow_result.get("found_match", False)

                classification_data.source = match_stage
                classification_data.closest_match = match_result

                # Handle exact/contained/similar matches
                if found_match and match_stage in [
                    "exact",
                    "contained",
                    "similar",
                ]:
                    meta_df_sub = meta_df[meta_df["accession_tag"] == match_result]

                    if not meta_df_sub.empty:
                        family = (
                            meta_df_sub["familyName"].iloc[0]
                            if "familyName" in meta_df_sub.columns
                            else None
                        )
                        classification_data.family = family

                        navis = (
                            meta_df_sub["navis_name"].iloc[0]
                            if "navis_name" in meta_df_sub.columns
                            else None
                        )
                        classification_data.navis = navis

                        haplotype = (
                            meta_df_sub["haplotype_name"].iloc[0]
                            if "haplotype_name" in meta_df_sub.columns
                            else None
                        )
                        classification_data.haplotype = haplotype

                        # Set confidence based on match type
                        if match_stage in ["exact", "contained", "similar"]:
                            classification_data.confidence = "High"
                        else:
                            classification_data.confidence = "Low"

                # Check for family/navis/haplotype classification from workflow
                if match_stage == "family" and match_result:
                    if isinstance(match_result, dict) and "family" in match_result:
                        family_name = match_result["family"]
                    elif (
                        isinstance(match_result, tuple)
                        and len(match_result) > 0
                        and isinstance(match_result[0], dict)
                        and "family" in match_result[0]
                    ):
                        family_name = match_result[0]["family"]
                    else:
                        family_name = match_result

                    if family_name:
                        # Extract aln_length and evalue from result if available, otherwise use defaults
                        top_aln_length = 0
                        top_evalue = 1.0

                        if isinstance(match_result, dict):
                            top_aln_length = match_result.get("aln_length", 0)
                            top_evalue = match_result.get("evalue", 1.0)
                        elif (
                            isinstance(match_result, tuple)
                            and len(match_result) > 0
                            and isinstance(match_result[0], dict)
                        ):
                            top_aln_length = match_result[0].get("aln_length", 0)
                            top_evalue = match_result[0].get("evalue", 1.0)

                        classification_data.source = "hmmsearch"
                        classification_data.family = family_name
                        classification_data.match_details = f"Captain gene match with length {top_aln_length}, Evalue: {top_evalue}"
                        classification_data.confidence = (
                            "High"
                            if top_aln_length > 80 and top_evalue < 0.001
                            else "Low"
                        )

                if match_stage == "navis" and match_result:
                    navis_name = match_result
                    # Make sure family_name is defined
                    family_name = workflow_result.get("family", "Unknown")
                    classification_data.source = "captain clustering"
                    classification_data.family = family_name
                    classification_data.navis = navis_name
                    classification_data.confidence = "Medium"

                if match_stage == "haplotype" and match_result:
                    haplotype_name = match_result
                    # Make sure family_name and navis_name are defined
                    family_name = workflow_result.get("family", "Unknown")
                    navis_name = workflow_result.get("navis", "Unknown")
                    classification_data.source = "k-mer similarity"
                    classification_data.family = family_name
                    classification_data.navis = navis_name
                    classification_data.haplotype = haplotype_name
                    classification_data.confidence = "Medium"

        except Exception as e:
            logger.error(f"Error processing classification workflow for sequence: {e}")

        # Don't process BLAST results here since classification workflow will handle everything
        # The classification workflow will run in the update_classification_workflow_state callback
        # and will use BLAST results as part of its comprehensive classification process

        # Only use BLAST fallback if sequence is too short for classification workflow
        sequence_length = len(query_seq)
        skip_classification = sequence_length < 5000

        if (
            skip_classification
            and classification_data.source is None
            and blast_df is not None
        ):
            # Process the BLAST results
            try:
                # Sort by evalue (ascending) and pident (descending) to get best hits
                blast_df = blast_df.sort_values(
                    ["evalue", "pident"], ascending=[True, False]
                )
                top_hit = blast_df.iloc[0]
                logger.debug(f"Top hit: {top_hit}")

                top_evalue = float(top_hit["evalue"])
                top_aln_length = int(top_hit["aln_length"])
                top_pident = float(top_hit["pident"])

                if top_pident >= 90:
                    # look up family name from accession tag
                    hit_IDs = top_hit["hit_IDs"]
                    # Convert hit_IDs to a list if it's a string
                    hit_IDs_list = [hit_IDs] if isinstance(hit_IDs, str) else hit_IDs
                    meta_df_sub = meta_df[meta_df["accession_tag"].isin(hit_IDs_list)]

                    if not meta_df_sub.empty:
                        top_family = meta_df_sub["familyName"].iloc[0]
                        navis = (
                            meta_df_sub["navis_name"].iloc[0]
                            if "navis_name" in meta_df_sub.columns
                            else None
                        )
                        haplotype = (
                            meta_df_sub["haplotype_name"].iloc[0]
                            if "haplotype_name" in meta_df_sub.columns
                            else None
                        )

                        classification_data.source = "blast_hit"
                        classification_data.family = top_family
                        classification_data.navis = navis
                        classification_data.haplotype = haplotype
                        classification_data.closest_match = hit_IDs
                        classification_data.match_details = f"length {top_aln_length}bp with {top_pident:.1f}% identity. E-value: {top_evalue}"
                        if top_pident >= 90 and top_aln_length > query_length * 0.8:
                            classification_data.confidence = "High"
                        elif top_pident >= 50 and top_aln_length > query_length * 0.5:
                            classification_data.confidence = "Medium"
                        else:
                            classification_data.confidence = "Low"

            except Exception as e:
                logger.error(f"Error processing BLAST results: {e}")

        # If still no classification and sequence is too short for workflow, try HMMER
        if skip_classification and classification_data.source is None:
            try:
                hmmer_results = run_hmmer_search_task(
                    query_header=query_header,
                    query_seq=query_seq,
                    query_type=query_type,
                    eval_threshold=evalue_threshold,
                )

                # Handle the new format with protein_content instead of protein_file
                captain_results = None
                if hmmer_results:
                    captain_results = hmmer_results.get("results")

                    # If we need the protein file for something else, recreate it
                    if (
                        "protein_content" in hmmer_results
                        and hmmer_results["protein_content"]
                    ):
                        protein_temp_file = tempfile.NamedTemporaryFile(
                            suffix=".faa", delete=False
                        ).name
                        with open(protein_temp_file, "w") as f:
                            f.write(hmmer_results["protein_content"])
                        # If needed: hmmer_results["protein_file"] = protein_temp_file

                if captain_results and len(captain_results) > 0:
                    first_result = captain_results[0]  # Get just the first dictionary
                    captain_family, aln_length, evalue = select_ship_family(
                        first_result
                    )

                    if captain_family:
                        if aln_length >= 80 and evalue <= 0.001:
                            classification_data.source = "hmmsearch"
                            classification_data.family = captain_family
                            classification_data.match_details = f"Captain gene match (length = {aln_length}, e-value = {evalue})"
                            classification_data.confidence = "High"
                        elif aln_length >= 50:
                            classification_data.source = "hmmsearch"
                            classification_data.family = captain_family
                            classification_data.match_details = f"Captain gene partial match (length = {aln_length}, e-value = {evalue})"
                            classification_data.confidence = "Medium"
                        else:
                            classification_data.source = "hmmsearch"
                            classification_data.family = captain_family
                            classification_data.match_details = f"Captain gene partial match (length = {aln_length}, e-value = {evalue})"
                            classification_data.confidence = "Low"

            except Exception as e:
                logger.error(f"Error processing HMMER results: {e}")

        # Update classification data
        if classification_data.source:
            class_dict.classification = classification_data

        class_dict.processed = True

        # Return structured results - convert to dict for backward compatibility
        return class_dict.to_dict()

    except Exception as e:
        class_dict.error = str(e)
        logger.error(f"Error in process_single_sequence: {e}")
        return class_dict.to_dict()


@callback(
    [
        Output("blast-results-store", "data"),
        Output("classification-interval", "disabled", allow_duplicate=True),
        Output("workflow-state-store", "data"),
        Output("submit-button", "loading", allow_duplicate=True),
        Output("results-loading-overlay", "visible", allow_duplicate=True),
    ],
    [
        Input("submission-id-store", "data"),
    ],
    [
        State("blast-sequences-store", "data"),
        State("evalue-threshold", "value"),
        State("blast-fasta-upload", "contents"),  # Add file contents as state
    ],
    running=[
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
def process_sequences(submission_id, seq_list, evalue_threshold, file_contents):
    """Process only the first sequence initially, for both text input and file uploads"""
    if not submission_id:
        logger.warning("No submission_id provided to process_sequences")
        raise PreventUpdate

    if not seq_list:
        # If we have no seq_list but do have file contents, parse the file directly
        if file_contents:
            logger.debug("No seq_list but file_contents available - parsing directly")
            try:
                input_type, direct_seq_list, n_seqs, error = check_input(
                    query_text_input=None, query_file_contents=file_contents
                )

                if error or not direct_seq_list or len(direct_seq_list) == 0:
                    logger.warning(f"Failed to parse file contents directly: {error}")
                    error_state = {
                        "complete": True,
                        "error": "Failed to parse sequence data",
                        "status": "failed",
                        "task_id": str(submission_id),
                        "found_match": False,
                        "match_stage": None,
                        "match_result": None,
                        "workflow_started": True,
                        "current_stage": None,
                        "current_stage_idx": 0,
                        "start_time": 0.0,
                        "stages": {},
                        "class_dict": {},
                    }
                    return (
                        None,
                        True,
                        error_state,
                        False,
                        False,
                    )  # Set loading to False when done

                logger.debug(
                    f"Successfully parsed {len(direct_seq_list)} sequences directly from file contents"
                )
                seq_list = direct_seq_list
            except Exception as e:
                logger.error(f"Error parsing file contents: {e}")
                error_state = {
                    "complete": True,
                    "error": "Error parsing sequence data",
                    "status": "failed",
                    "task_id": str(submission_id),
                    "found_match": False,
                    "match_stage": None,
                    "match_result": None,
                    "workflow_started": True,
                    "current_stage": None,
                    "current_stage_idx": 0,
                    "start_time": 0.0,
                    "stages": {},
                    "class_dict": {},
                }
                return (
                    None,
                    True,
                    error_state,
                    False,
                    False,
                )  # Set loading to False on error too
        else:
            logger.warning("No seq_list or file_contents provided to process_sequences")
            error_state = {
                "complete": True,
                "error": "No sequence data available",
                "status": "failed",
                "task_id": str(submission_id) if submission_id else None,
                "found_match": False,
                "match_stage": None,
                "match_result": None,
                "workflow_started": True,
                "current_stage": None,
                "current_stage_idx": 0,
                "start_time": 0.0,
                "stages": {},
                "class_dict": {},
            }
            return (
                None,
                True,
                error_state,
                False,
                False,
            )  # Set loading to False when done

    logger.debug(
        f"Processing sequence submission with ID: {submission_id}, sequences: {len(seq_list)}"
    )

    try:
        blast_results = None
        classification_interval_disabled = True  # Always disable the interval
        workflow_state = None

        # Process only the first sequence to start
        if seq_list and len(seq_list) > 0:
            # Log sequence details for debugging
            first_seq = seq_list[0]
            logger.debug(
                f"Processing first sequence: header={first_seq.get('header', 'unknown')[:30]}..., length={len(first_seq.get('sequence', ''))}"
            )

            # Process sequence 0
            sequence_result = process_single_sequence(seq_list[0], evalue_threshold)

            if sequence_result:
                logger.debug(
                    f"Sequence result: processed={sequence_result.get('processed', False)}, blast_file={sequence_result.get('blast_file') is not None}"
                )

                # Structure the blast results properly
                blast_results = {
                    "processed_sequences": [0],  # Indices of processed sequences
                    "sequence_results": {
                        "0": sequence_result
                    },  # Results keyed by sequence index
                    "total_sequences": len(seq_list),  # Store total number of sequences
                }

                # Determine if we should enable classification interval
                sequence_length = len(sequence_result.get("sequence", ""))
                skip_classification = (
                    sequence_length < 5000 or sequence_result.get("error") is not None
                )

                logger.debug(
                    f"Classification decision: skip={skip_classification}, seq_length={sequence_length}"
                )

                # Default to disabled - always
                classification_interval_disabled = True

                # Set up workflow state if needed
                if not skip_classification and sequence_result.get("blast_content"):
                    # Write the file content to include in class_dict
                    tmp_query_fasta = None

                    try:
                        tmp_query_fasta = write_temp_fasta(
                            seq_list[0].get("header", "query"),
                            seq_list[0].get("sequence", ""),
                        )
                    except Exception as e:
                        logger.error(f"Error writing FASTA to temp file: {e}")
                        # If we can't write the file, skip classification
                        skip_classification = True

                    if not skip_classification:
                        # Initialize a workflow state structure for Dash (without large data)
                        workflow_state = {
                            "task_id": str(submission_id),
                            "status": "initialized",
                            "complete": False,
                            "workflow_started": False,  # Will be set to True when processing starts
                            "current_stage": None,
                            "current_stage_idx": 0,
                            "error": None,
                            "found_match": False,
                            "match_stage": None,
                            "match_result": None,
                            "start_time": time.time(),
                            "stages": {
                                stage["id"]: {"progress": 0, "complete": False}
                                for stage in WORKFLOW_STAGES
                            },
                            "class_dict": {
                                "seq_type": seq_list[0].get("type", "nucl"),
                                "fasta_file": tmp_query_fasta,
                                "has_blast_results": sequence_result.get("blast_df")
                                is not None,
                                "fetch_ship_params": {
                                    "curated": False,
                                    "with_sequence": True,
                                    "dereplicate": True,
                                },
                                "fetch_captain_params": {
                                    "curated": True,
                                    "with_sequence": True,
                                },
                            },
                        }

                        # Since we're no longer using Celery, we don't need the interval
                        # always set to disabled since we run synchronously
                        classification_interval_disabled = True

            else:
                logger.warning(
                    "No sequence result returned from process_single_sequence"
                )
                blast_results = {
                    "processed_sequences": [0],
                    "sequence_results": {
                        "0": {
                            "processed": False,
                            "error": "Failed to process sequence",
                            "sequence": seq_list[0].get("sequence", ""),
                        }
                    },
                    "total_sequences": len(seq_list),
                }

        logger.debug(
            f"Completed process_sequences: has_blast_results={blast_results is not None}, interval_disabled={classification_interval_disabled}"
        )
        return (
            blast_results,
            classification_interval_disabled,
            workflow_state,
            False,
            False,
        )  # Set loading to False when done
    except Exception as e:
        logger.error(f"Error in process_sequences: {str(e)}")
        # Return basic data on error
        error_state = {
            "complete": True,
            "error": str(e),
            "status": "failed",
            "task_id": str(submission_id) if submission_id else None,
        }
        return (
            None,
            True,
            error_state,
            False,
            False,
        )  # Set loading to False on error too


@callback(
    Output("blast-results-store", "data", allow_duplicate=True),
    Input("blast-tabs", "active_tab"),
    [
        State("blast-results-store", "data"),
        State("blast-sequences-store", "data"),
        State("evalue-threshold", "value"),
    ],
    prevent_initial_call=True,
)
def process_additional_sequence(active_tab, results_store, seq_list, evalue_threshold):
    """Process a sequence when its tab is selected if not already processed"""
    try:
        if not active_tab or not results_store or not seq_list:
            logger.warning(
                f"Missing data for tab processing: active_tab={active_tab is not None}, results_store={results_store is not None}, seq_list={seq_list is not None}"
            )
            raise PreventUpdate

        # Extract tab index from tab id (e.g., "tab-2" -> 2)
        try:
            tab_idx = (
                int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0
            )
        except (ValueError, IndexError) as e:
            logger.error(f"Invalid tab format '{active_tab}': {e}")
            raise PreventUpdate

        # Check if this sequence has already been processed
        processed_sequences = results_store.get("processed_sequences", [])
        logger.debug(
            f"Tab {tab_idx} clicked. Processed sequences: {processed_sequences}"
        )

        if tab_idx in processed_sequences:
            # Already processed this sequence
            logger.debug(f"Tab {tab_idx} already processed, skipping")
            raise PreventUpdate

        # Ensure tab_idx is valid
        if tab_idx >= len(seq_list):
            logger.error(
                f"Tab index {tab_idx} out of range (seq_list length: {len(seq_list)})"
            )
            raise PreventUpdate

        # Process this sequence
        logger.info(
            f"Processing sequence for tab {tab_idx}: {seq_list[tab_idx].get('header', 'unknown')[:50]}..."
        )
        sequence_result = process_single_sequence(seq_list[tab_idx], evalue_threshold)

        if not sequence_result:
            logger.error(
                f"Failed to process sequence for tab {tab_idx} - no result returned"
            )
            raise PreventUpdate

        if sequence_result.get("error"):
            logger.error(
                f"Error processing sequence for tab {tab_idx}: {sequence_result.get('error')}"
            )

        # Check if this sequence should get classification workflow (≥5000bp)
        sequence_length = len(seq_list[tab_idx].get("sequence", ""))
        should_classify = sequence_length >= 5000

        if should_classify and sequence_result.get("blast_content"):
            logger.info(
                f"Running classification workflow for tab {tab_idx} (length: {sequence_length})"
            )
            try:
                from src.utils.blast_data import (
                    BlastData,
                    FetchShipParams,
                    FetchCaptainParams,
                )
                from src.tasks import run_classification_workflow_task
                from src.database.sql_manager import fetch_meta_data
                import tempfile

                # Create temporary FASTA file for the workflow
                tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
                with open(tmp_fasta, "w") as f:
                    f.write(
                        f">{seq_list[tab_idx].get('header', 'query')}\n{seq_list[tab_idx].get('sequence', '')}\n"
                    )

                # Set up workflow data
                upload_data = BlastData(
                    seq_type=seq_list[tab_idx].get("type", "nucl"),
                    fetch_ship_params=FetchShipParams(
                        curated=False, with_sequence=True, dereplicate=True
                    ),
                    fetch_captain_params=FetchCaptainParams(
                        curated=True, with_sequence=True
                    ),
                    fasta_file=tmp_fasta,
                    blast_df=sequence_result.get("blast_df"),
                )

                meta_df = fetch_meta_data()
                meta_dict = meta_df.to_dict("records") if meta_df is not None else None

                # Run classification workflow
                workflow_result = run_classification_workflow_task(
                    class_dict=upload_data, meta_dict=meta_dict
                )

                # If workflow found a match, update sequence results with classification
                if (
                    workflow_result
                    and workflow_result.get("complete", False)
                    and workflow_result.get("found_match", False)
                ):
                    match_stage = workflow_result.get("match_stage")
                    match_accession = workflow_result.get("match_result")

                    logger.info(
                        f"Workflow found {match_stage} match for tab {tab_idx}: {match_accession}"
                    )

                    # Create classification data
                    classification_data = {
                        "source": match_stage,
                        "closest_match": match_accession,
                        "confidence": "High"
                        if match_stage in ["exact", "contained"]
                        else "Medium",
                    }

                    # Look up metadata for the matched accession
                    try:
                        if meta_df is not None and not meta_df.empty:
                            meta_match = meta_df[
                                meta_df["accession_tag"] == match_accession
                            ]
                            if not meta_match.empty:
                                # Add family information
                                if "familyName" in meta_match.columns:
                                    family = meta_match["familyName"].iloc[0]
                                    if family and family != "None":
                                        classification_data["family"] = family

                                # Add navis information
                                if "navis_name" in meta_match.columns:
                                    navis = meta_match["navis_name"].iloc[0]
                                    if navis and navis != "None":
                                        classification_data["navis"] = navis

                                # Add haplotype information
                                if "haplotype_name" in meta_match.columns:
                                    haplotype = meta_match["haplotype_name"].iloc[0]
                                    if haplotype and haplotype != "None":
                                        classification_data["haplotype"] = haplotype

                                # Add match details
                                if match_stage == "exact":
                                    classification_data["match_details"] = (
                                        f"Exact sequence match to {match_accession}"
                                    )
                                elif match_stage == "contained":
                                    classification_data["match_details"] = (
                                        f"Query sequence contained within {match_accession}"
                                    )
                                elif match_stage == "similar":
                                    classification_data["match_details"] = (
                                        f"High similarity to {match_accession}"
                                    )
                                else:
                                    classification_data["match_details"] = (
                                        f"{match_stage.replace('_', ' ').title()} match to {match_accession}"
                                    )
                    except Exception as e:
                        logger.error(
                            f"Error looking up metadata for {match_accession}: {e}"
                        )

                    # Update sequence results with classification
                    sequence_result["classification"] = classification_data
                    logger.info(
                        f"Updated tab {tab_idx} with classification: {classification_data}"
                    )

            except Exception as e:
                logger.error(
                    f"Error running classification workflow for tab {tab_idx}: {e}"
                )
        else:
            if not should_classify:
                logger.debug(
                    f"Skipping classification for tab {tab_idx} - sequence too short ({sequence_length}bp)"
                )
            else:
                logger.debug(
                    f"Skipping classification for tab {tab_idx} - no BLAST content"
                )

        # Update the results store
        updated_results = dict(results_store)
        updated_results["processed_sequences"].append(tab_idx)
        updated_results["sequence_results"][str(tab_idx)] = sequence_result

        logger.info(f"Successfully processed and updated results for tab {tab_idx}")
        return updated_results

    except PreventUpdate:
        raise
    except Exception as e:
        logger.error(
            f"Unexpected error in process_additional_sequence for tab {active_tab}: {e}"
        )
        logger.exception("Full traceback:")
        raise PreventUpdate


########################################################
# Tabbed Output
########################################################
@callback(
    Output("tab-content", "children"),
    [
        Input("blast-tabs", "active_tab"),
        Input("blast-results-store", "data"),
    ],
    prevent_initial_call=True,
)
def render_tab_content(active_tab, results_store):
    """Render the content for the current tab"""
    from src.utils.classification_utils import create_classification_output

    try:
        if not active_tab or not results_store:
            logger.warning(
                f"Missing data for tab rendering: active_tab={active_tab}, results_store={results_store is not None}"
            )
            raise PreventUpdate

        # Extract tab index
        try:
            tab_idx = (
                int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0
            )
        except (ValueError, IndexError):
            logger.error(f"Invalid tab format: {active_tab}")
            return dmc.Alert(
                title="Error",
                children="Invalid tab format",
                color="red",
                variant="filled",
            )

        # Check if this sequence has been processed
        processed_sequences = results_store.get("processed_sequences", [])
        logger.debug(
            f"Rendering tab {tab_idx}. Processed sequences: {processed_sequences}"
        )

        if tab_idx not in processed_sequences:
            # Not processed yet - show loading state
            logger.debug(f"Tab {tab_idx} not processed yet, showing loading state")
            return dmc.Stack(
                [
                    dmc.Center(dmc.Loader(size="xl")),
                    dmc.Text("Processing sequence...", size="lg"),
                ]
            )

        # Get results for this sequence
        sequence_results = results_store.get("sequence_results", {}).get(str(tab_idx))
        if not sequence_results:
            logger.error(f"No sequence results found for tab {tab_idx}")
            return dmc.Alert(
                title="Error",
                children=f"No results found for sequence {tab_idx + 1}. Please try clicking on this tab again to reprocess.",
                color="red",
                variant="filled",
            )

        # Check for errors in the results
        if "error" in sequence_results and sequence_results["error"] is not None:
            logger.error(
                f"Error in sequence results for tab {tab_idx}: {sequence_results['error']}"
            )
            return seq_processing_error_alert(sequence_results["error"])

        # Render classification results (workflow state not available for tabbed interface)
        classification_output = create_classification_output(
            sequence_results, workflow_state=None
        )

        # Render BLAST results
        blast_container = create_blast_container(sequence_results, tab_id=tab_idx)

        logger.debug(f"Successfully rendered content for tab {tab_idx}")
        return [
            classification_output,
            progress_section,
            # BLAST results section
            dmc.Stack(
                [
                    blast_container,
                ]
            ),
        ]

    except PreventUpdate:
        raise
    except Exception as e:
        logger.error(
            f"Unexpected error in render_tab_content for tab {active_tab}: {e}"
        )
        logger.exception("Full traceback:")
        return dmc.Alert(
            title="Error Rendering Tab",
            children=f"An unexpected error occurred while rendering this tab: {str(e)}",
            color="red",
            variant="filled",
        )


@callback(
    Output("blast-active-tab", "data"),
    Input("blast-tabs", "active_tab"),
    prevent_initial_call=True,
)
def update_active_tab(active_tab):
    if not active_tab:
        raise PreventUpdate

    # Extract tab index from tab id (e.g., "tab-2" -> 2)
    tab_idx = int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0
    return tab_idx


########################################################
# Clientside callbacks for BlasterJS
########################################################

# Define the clientside JavaScript function directly in the callback
clientside_callback(
    """
    function(data, active_tab_idx) {
        console.log("BlasterJS callback triggered", data, "active tab:", active_tab_idx);
        
        // Check if we have valid data
        if (!data || !data.blast_text) {
            console.log("No valid BLAST data available");
            return window.dash_clientside.no_update;
        }
        
        try {
            console.log("Initializing BlasterJS with text length:", data.blast_text.length);
            
            // Find the correct container based on active tab
            let containerId = active_tab_idx !== null ? 
                `blast-container-${active_tab_idx}` : 'blast-container';
            
            // Create a function that will attempt to initialize BlasterJS
            const initializeBlasterJS = (attempts = 0, maxAttempts = 5) => {
                // Use standard ID selector
                let container = document.getElementById(containerId);
                
                // If no container is found and we haven't exceeded max attempts, retry
                if (!container && attempts < maxAttempts) {
                    console.log(`Container ${containerId} not found, retrying in 100ms (attempt ${attempts + 1}/${maxAttempts})`);
                    setTimeout(() => initializeBlasterJS(attempts + 1, maxAttempts), 100);
                    return;
                }
                
                // If still no container after max attempts, try fallback options
                if (!container) {
                    console.log(`No container found with ID: ${containerId} after ${maxAttempts} attempts`);
                    // Try to find the default container first since this is likely a single sequence view
                    container = document.getElementById('blast-container');
                    
                    if (!container) {
                        console.log("No default container found, trying by class");
                        // If still not found, try to find a container with class blast-container 
                        let containers = document.getElementsByClassName('blast-container');
                        if (containers.length > 0) {
                            console.log("Found container by class instead");
                            container = containers[0];
                        } else {
                            console.error("No blast containers found in the DOM");
                            return;
                        }
                    }
                }
                
                console.log("Found container:", container);
                
                // Clear existing content first
                container.innerHTML = '';
                
                // If we have empty or invalid blast text, show a message
                if (!data.blast_text || !data.blast_text.trim() || data.blast_text === "BLAST results too large to display") {
                    let messageText = data.blast_text === "BLAST results too large to display" ?
                        "The BLAST results are too large to display in the browser." :
                        "No BLAST results to display";
                        
                    container.innerHTML = '<div style="padding: 20px; text-align: center; color: #666;">' + messageText + '</div>';
                    return;
                }
                
                // Create the title element
                const titleElement = document.createElement('h2');
                titleElement.innerHTML = 'BLAST Results';
                titleElement.style.marginTop = '15px';
                titleElement.style.marginBottom = '20px';
                titleElement.style.textAlign = 'left';
                titleElement.style.width = '100%';
                
                // Create unique IDs for this container to avoid conflicts in multi-tab view
                const uniquePrefix = containerId.replace('blast-container', 'blast');
                const alignmentsDivId = `${uniquePrefix}-multiple-alignments`;
                const tablesDivId = `${uniquePrefix}-alignments-table`;
                
                // Create simple divs for BlasterJS with explicit left alignment
                const alignmentsDiv = document.createElement('div');
                alignmentsDiv.id = alignmentsDivId;
                alignmentsDiv.style.textAlign = 'left';
                alignmentsDiv.style.width = '100%';
                
                const tableDiv = document.createElement('div');
                tableDiv.id = tablesDivId;
                tableDiv.style.textAlign = 'left';
                tableDiv.style.width = '100%';
                
                // Add elements to the container
                container.appendChild(titleElement);
                container.appendChild(alignmentsDiv);
                container.appendChild(tableDiv);
                
                // Add a style element to ensure BlasterJS output is properly aligned
                const styleEl = document.createElement('style');
                styleEl.textContent = `
                    #${alignmentsDivId}, #${tablesDivId} {
                        text-align: left !important;
                        margin-left: 0 !important;
                        padding-left: 0 !important;
                    }
                    #${alignmentsDivId} div, #${tablesDivId} div,
                    #${alignmentsDivId} table, #${tablesDivId} table {
                        text-align: left !important;
                        margin-left: 0 !important;
                    }
                    .alignment-viewer {
                        text-align: left !important;
                    }
                `;
                container.appendChild(styleEl);
                
                // Basic BlasterJS initialization
                try {
                    var blasterjs = require("biojs-vis-blasterjs");
                    console.log("BlasterJS loaded successfully");
                    var instance = new blasterjs({
                        string: data.blast_text,
                        multipleAlignments: alignmentsDivId,
                        alignmentsTable: tablesDivId
                    });
                    console.log("BlasterJS initialized successfully for", containerId);
                    
                    // Store instance ID in a data attribute for potential future reference
                    container.dataset.blasterjsInstance = uniquePrefix;
                    container.dataset.initialized = "true";
                    
                    // Additional styling fix after BlasterJS renders
                    setTimeout(function() {
                        const tables = container.querySelectorAll('table');
                        tables.forEach(function(table) {
                            table.style.marginLeft = '0';
                            table.style.textAlign = 'left';
                        });
                    }, 100);
                } catch (blasterError) {
                    console.error("Error initializing BlasterJS library:", blasterError);
                    container.innerHTML += "<div style='color:red;'>Error initializing BLAST viewer: " + blasterError + "</div>";
                    container.dataset.initialized = "error";
                }
            };
            
            // Start the initialization process
            initializeBlasterJS();
            
            return window.dash_clientside.no_update;
        } catch (error) {
            console.error('Overall error in callback:', error);
            // Error handling - display error in a div
            return window.dash_clientside.no_update;
        }
    }
    """,
    Output("blast-container", "children"),
    [Input("processed-blast-store", "data"), Input("blast-active-tab", "data")],
    prevent_initial_call=True,
)

# Add a clientside callback for tab-specific BLAST visualization
clientside_callback(
    """
    function(active_tab, blast_data) {
        console.log("Tab BlasterJS callback triggered", 
                    active_tab, 
                    blast_data ? (blast_data.blast_text ? "has text" : "no text") : "no data");
        
        if (!active_tab) {
            console.warn("No active tab provided");
            return window.dash_clientside.no_update;
        }
        
        if (!blast_data) {
            // console.warn("No blast data available");
            return window.dash_clientside.no_update;
        }
        
        if (!blast_data.blast_text) {
            console.warn("No blast text in data");
            return window.dash_clientside.no_update;
        }

        try {
            // Extract the tab index from the active tab ID
            const tabIdx = parseInt(active_tab.split("-")[1]);
            if (isNaN(tabIdx)) {
                console.error("Invalid tab index:", active_tab);
                return window.dash_clientside.no_update;
            }

            // Find the correct container based on active tab
            const containerId = `blast-container-${tabIdx}`;
            
            // Create a function that will attempt to initialize BlasterJS
            const initializeBlasterJS = (attempts = 0, maxAttempts = 5) => {
                let container = document.getElementById(containerId);
                
                // If no container is found and we haven't exceeded max attempts, retry
                if (!container && attempts < maxAttempts) {
                    console.log(`Container ${containerId} not found, retrying in 100ms (attempt ${attempts + 1}/${maxAttempts})`);
                    setTimeout(() => initializeBlasterJS(attempts + 1, maxAttempts), 100);
                    return;
                }
                
                // If still no container after max attempts, try fallback options
                if (!container) {
                    console.log(`No container found with ID: ${containerId} after ${maxAttempts} attempts`);
                    // Try fallback options
                    const containers = document.getElementsByClassName('blast-container');
                    if (containers.length > 0) {
                        console.log("Found container by class instead");
                        container = containers[0];
                    } else {
                        const mainContainer = document.getElementById('blast-container');
                        if (mainContainer) {
                            console.log("Using main blast-container as fallback");
                            container = mainContainer;
                        } else {
                            console.error("No blast containers found in the DOM");
                            return;
                        }
                    }
                }

                // Only initialize if empty or not initialized yet
                if (container.children.length === 0 || !container.dataset.initialized) {
                    console.log(`Initializing BlasterJS for tab ${tabIdx}`);
                    
                    // Clear existing content
                    container.innerHTML = '';
                    
                    // Create the title element
                    const titleElement = document.createElement('h2');
                    titleElement.innerHTML = 'BLAST Results';
                    titleElement.style.marginTop = '15px';
                    titleElement.style.marginBottom = '20px';
                    titleElement.style.textAlign = 'left';
                    titleElement.style.width = '100%';
                    
                    // Create unique IDs for this container
                    const uniquePrefix = `blast-${tabIdx}`;
                    const alignmentsDivId = `${uniquePrefix}-multiple-alignments`;
                    const tablesDivId = `${uniquePrefix}-alignments-table`;
                    
                    // Create divs for BlasterJS
                    const alignmentsDiv = document.createElement('div');
                    alignmentsDiv.id = alignmentsDivId;
                    alignmentsDiv.style.textAlign = 'left';
                    alignmentsDiv.style.width = '100%';
                    
                    const tableDiv = document.createElement('div');
                    tableDiv.id = tablesDivId;
                    tableDiv.style.textAlign = 'left';
                    tableDiv.style.width = '100%';
                    
                    // Add elements to the container
                    container.appendChild(titleElement);
                    container.appendChild(alignmentsDiv);
                    container.appendChild(tableDiv);
                    
                    // Add styling
                    const styleEl = document.createElement('style');
                    styleEl.textContent = `
                        #${alignmentsDivId}, #${tablesDivId} {
                            text-align: left !important;
                            margin-left: 0 !important;
                            padding-left: 0 !important;
                        }
                        #${alignmentsDivId} div, #${tablesDivId} div,
                        #${alignmentsDivId} table, #${tablesDivId} table {
                            text-align: left !important;
                            margin-left: 0 !important;
                        }
                        .alignment-viewer {
                            text-align: left !important;
                        }
                    `;
                    container.appendChild(styleEl);
                    
                    // Validate blast text
                    if (!blast_data.blast_text.trim() || blast_data.blast_text === "BLAST results too large to display") {
                        let messageText = blast_data.blast_text === "BLAST results too large to display" ?
                            "The BLAST results are too large to display in the browser." :
                            "No BLAST results to display";
                        
                        const emptyDiv = document.createElement('div');
                        emptyDiv.innerHTML = '<div style="padding: 20px; text-align: center; color: #666;">' + messageText + '</div>';
                        container.appendChild(emptyDiv);
                        container.dataset.initialized = "true";
                        return;
                    }
                    
                    // Initialize BlasterJS
                    try {
                        // Check if biojs-vis-blasterjs is available
                        if (typeof require !== 'function') {
                            console.error("require function not available - can't load BlasterJS");
                            container.innerHTML += '<div style="color:orange;padding:10px;">Error: BlasterJS library not available</div>';
                            container.dataset.initialized = "error";
                            return;
                        }
                        
                        let blasterjs = require("biojs-vis-blasterjs");
                        if (!blasterjs) {
                            console.error("Failed to load BlasterJS library");
                            container.innerHTML += '<div style="color:orange;padding:10px;">Error loading BlasterJS library</div>';
                            container.dataset.initialized = "error";
                            return;
                        }
                        
                        let instance = new blasterjs({
                            string: blast_data.blast_text,
                            multipleAlignments: alignmentsDivId,
                            alignmentsTable: tablesDivId
                        });
                        
                        // Mark as initialized
                        container.dataset.initialized = "true";
                        
                        // Apply additional styling
                        setTimeout(function() {
                            const tables = container.querySelectorAll('table');
                            tables.forEach(function(table) {
                                table.style.marginLeft = '0';
                                table.style.textAlign = 'left';
                            });
                        }, 100);
                        
                        console.log(`BlasterJS initialized for tab ${tabIdx}`);
                    } catch (error) {
                        console.error("Error initializing BlasterJS:", error);
                        container.innerHTML += `<div style="color:red;padding:10px;">Error initializing BLAST viewer: ${error.toString()}</div>`;
                        container.dataset.initialized = "error";
                    }
                }
            };
            
            // Start the initialization process
            initializeBlasterJS();
            
        } catch (error) {
            console.error("Error in tab visualization callback:", error);
        }
        
        return window.dash_clientside.no_update;
    }
    """,
    Output("tab-content", "children", allow_duplicate=True),
    [Input("blast-tabs", "active_tab"), Input("processed-blast-store", "data")],
    prevent_initial_call=True,
)


########################################################
# Classification Output
########################################################
@callback(
    Output("classification-output", "children"),
    [
        Input("blast-results-store", "data"),
        Input("workflow-state-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_single_sequence_classification(blast_results_store, workflow_state):
    """Update classification output for single sequence submissions"""
    from src.utils.classification_utils import create_classification_output

    logger.info(
        f"update_single_sequence_classification CALLED: blast_results_store={blast_results_store is not None}, workflow_state={workflow_state is not None}"
    )

    if workflow_state:
        logger.info(
            f"Workflow state: complete={workflow_state.get('complete')}, found_match={workflow_state.get('found_match')}, match_stage={workflow_state.get('match_stage')}, match_result={workflow_state.get('match_result')}"
        )

    if not blast_results_store:
        logger.info("No blast_results_store, returning None")
        return None

    # Get the results for sequence 0
    sequence_results = blast_results_store.get("sequence_results", {}).get("0")
    if not sequence_results:
        logger.info("No sequence results for sequence 0, returning None")
        return None

    logger.info(
        f"Original sequence_results has classification: {'classification' in sequence_results}"
    )

    # Check if we have workflow results that should override the initial classification
    if (
        workflow_state
        and workflow_state.get("complete", False)
        and workflow_state.get("found_match", True)
    ):
        logger.info("Workflow completed with match found, updating classification")
        # Update sequence_results with workflow classification data
        match_result = workflow_state.get("match_result")
        if isinstance(match_result, dict):
            # Use workflow results
            sequence_results = sequence_results.copy()
            sequence_results["classification"] = match_result
            logger.info("Used workflow match_result dict as classification")
        elif workflow_state.get("match_stage") and workflow_state.get("match_result"):
            # Convert simple workflow results to classification format and look up metadata
            sequence_results = sequence_results.copy()
            match_stage = workflow_state.get("match_stage")
            match_accession = workflow_state.get("match_result")

            logger.info(
                f"Processing workflow result: {match_stage} -> {match_accession}"
            )

            # Create base classification
            classification_data = {
                "source": match_stage,
                "closest_match": match_accession,
                "confidence": "High"
                if match_stage in ["exact", "contained"]
                else "Medium",
            }

            # Look up additional metadata for the matched accession
            try:
                meta_df = fetch_meta_data()
                if meta_df is not None and not meta_df.empty:
                    meta_match = meta_df[meta_df["accession_tag"] == match_accession]
                    if not meta_match.empty:
                        logger.debug(f"Found metadata for {match_accession}")
                        # Add family information
                        if "familyName" in meta_match.columns:
                            family = meta_match["familyName"].iloc[0]
                            if family and family != "None":
                                classification_data["family"] = family

                        # Add navis information
                        if "navis_name" in meta_match.columns:
                            navis = meta_match["navis_name"].iloc[0]
                            if navis and navis != "None":
                                classification_data["navis"] = navis

                        # Add haplotype information
                        if "haplotype_name" in meta_match.columns:
                            haplotype = meta_match["haplotype_name"].iloc[0]
                            if haplotype and haplotype != "None":
                                classification_data["haplotype"] = haplotype

                        # Add match details
                        if match_stage == "exact":
                            classification_data["match_details"] = (
                                f"Exact sequence match to {match_accession}"
                            )
                        elif match_stage == "contained":
                            classification_data["match_details"] = (
                                f"Query sequence contained within {match_accession}"
                            )
                        elif match_stage == "similar":
                            classification_data["match_details"] = (
                                f"High similarity to {match_accession}"
                            )
                        else:
                            classification_data["match_details"] = (
                                f"{match_stage.replace('_', ' ').title()} match to {match_accession}"
                            )
                    else:
                        logger.warning(f"No metadata found for {match_accession}")

            except Exception as e:
                logger.error(f"Error looking up metadata for {match_accession}: {e}")

            sequence_results["classification"] = classification_data
            logger.info(f"Created classification_data: {classification_data}")
    else:
        logger.info(
            f"Workflow not complete or no match found. Complete: {workflow_state.get('complete') if workflow_state else 'N/A'}, Found match: {workflow_state.get('found_match') if workflow_state else 'N/A'}"
        )

    # Create the classification output using the existing function
    logger.info(
        f"Final sequence_results has classification: {'classification' in sequence_results}"
    )
    result = create_classification_output(sequence_results, workflow_state)
    logger.info(f"create_classification_output returned: {type(result)}")
    return result


@callback(
    [
        Output("workflow-state-store", "data", allow_duplicate=True),
        Output("results-loading-overlay", "visible", allow_duplicate=True),
        Output("blast-results-store", "data", allow_duplicate=True),
    ],
    [
        # Remove dependency on interval n_intervals
        Input("workflow-state-store", "data"),
    ],
    [
        State("blast-results-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_classification_workflow_state(workflow_state, blast_results_store):
    """Initialize workflow state (one-time operation)"""
    # If no workflow state, nothing to do
    if workflow_state is None:
        raise PreventUpdate

    # Return early if workflow is complete, has error, or already started
    if (
        workflow_state.get("complete", False)
        or workflow_state.get("error") is not None
        or workflow_state.get("workflow_started", False)
    ):
        # If workflow is complete, also set loading overlay to false
        if (
            workflow_state.get("complete", False)
            or workflow_state.get("error") is not None
        ):
            return workflow_state, False, blast_results_store
        raise PreventUpdate

    # Check if task is running
    task_id = workflow_state.get("task_id")
    if not task_id:
        logger.warning("No task ID found in workflow state")
        workflow_state = workflow_state.copy()
        workflow_state["error"] = "Missing task ID"
        workflow_state["complete"] = True
        workflow_state["status"] = "failed"
        return workflow_state, False, blast_results_store

    try:
        # Start the classification workflow - ONE TIME ONLY
        if workflow_state.get("task_id") and not workflow_state.get("workflow_started"):
            class_dict = workflow_state.get("class_dict", {})

            from src.utils.blast_data import (
                BlastData,
                FetchShipParams,
                FetchCaptainParams,
            )

            # Extract BLAST data from the blast_results_store
            blast_df = None
            if blast_results_store and "sequence_results" in blast_results_store:
                sequence_results = blast_results_store["sequence_results"].get("0")
                if sequence_results and "blast_df" in sequence_results:
                    blast_df = sequence_results["blast_df"]

            upload_data = BlastData(
                seq_type=class_dict.get("seq_type", "nucl"),
                fetch_ship_params=FetchShipParams(
                    **class_dict.get("fetch_ship_params", {})
                ),
                fetch_captain_params=FetchCaptainParams(
                    **class_dict.get("fetch_captain_params", {})
                ),
                fasta_file=class_dict.get("fasta_file"),
                blast_df=blast_df,  # Pass BLAST results to the workflow
            )

            meta_df = fetch_meta_data()
            meta_dict = meta_df.to_dict("records") if meta_df is not None else None

            logger.debug("Running classification workflow directly (no Celery)")
            # Run the workflow function directly
            result = run_classification_workflow_task(
                class_dict=upload_data, meta_dict=meta_dict
            )

            logger.debug(f"Workflow result type: {type(result)}, content: {result}")

            # Use a copy to modify
            workflow_state = workflow_state.copy()

            # Mark as started so we don't run it again
            workflow_state["workflow_started"] = True

            # If we got a result back, merge it with the workflow state
            if isinstance(result, dict):
                logger.debug(
                    f"Merging workflow result into state. Keys: {list(result.keys())}"
                )
                for key, value in result.items():
                    if value is not None or key not in workflow_state:
                        workflow_state[key] = value
                        logger.debug(f"Set workflow_state['{key}'] = {value}")

            # Ensure the status is set properly
            workflow_state["status"] = "complete"
            workflow_state["complete"] = True

            logger.debug(
                f"Final workflow state: complete={workflow_state.get('complete')}, found_match={workflow_state.get('found_match')}, match_stage={workflow_state.get('match_stage')}, match_result={workflow_state.get('match_result')}"
            )

            # IMPORTANT: Also update the sequence results with classification data for tabbed interface
            if (
                workflow_state.get("complete", False)
                and workflow_state.get("found_match", False)
                and blast_results_store
                and "sequence_results" in blast_results_store
            ):
                logger.info(
                    "Updating sequence results with workflow classification for tabbed interface"
                )

                # Make a copy of blast_results_store to avoid modifying the original
                blast_results_store = blast_results_store.copy()
                blast_results_store["sequence_results"] = blast_results_store[
                    "sequence_results"
                ].copy()

                # Get the sequence results for sequence 0
                sequence_results = blast_results_store["sequence_results"].get("0")
                if sequence_results:
                    # Make a copy of sequence_results to avoid modifying the original
                    sequence_results = sequence_results.copy()

                    match_stage = workflow_state.get("match_stage")
                    match_accession = workflow_state.get("match_result")

                    # Create classification data
                    classification_data = {
                        "source": match_stage,
                        "closest_match": match_accession,
                        "confidence": "High"
                        if match_stage in ["exact", "contained"]
                        else "Medium",
                    }

                    # Look up metadata for the matched accession
                    try:
                        meta_df = fetch_meta_data()
                        if meta_df is not None and not meta_df.empty:
                            meta_match = meta_df[
                                meta_df["accession_tag"] == match_accession
                            ]
                            if not meta_match.empty:
                                logger.debug(f"Found metadata for {match_accession}")
                                # Add family information
                                if "familyName" in meta_match.columns:
                                    family = meta_match["familyName"].iloc[0]
                                    if family and family != "None":
                                        classification_data["family"] = family

                                # Add navis information
                                if "navis_name" in meta_match.columns:
                                    navis = meta_match["navis_name"].iloc[0]
                                    if navis and navis != "None":
                                        classification_data["navis"] = navis

                                # Add haplotype information
                                if "haplotype_name" in meta_match.columns:
                                    haplotype = meta_match["haplotype_name"].iloc[0]
                                    if haplotype and haplotype != "None":
                                        classification_data["haplotype"] = haplotype

                                # Add match details
                                if match_stage == "exact":
                                    classification_data["match_details"] = (
                                        f"Exact sequence match to {match_accession}"
                                    )
                                elif match_stage == "contained":
                                    classification_data["match_details"] = (
                                        f"Query sequence contained within {match_accession}"
                                    )
                                elif match_stage == "similar":
                                    classification_data["match_details"] = (
                                        f"High similarity to {match_accession}"
                                    )
                                else:
                                    classification_data["match_details"] = (
                                        f"{match_stage.replace('_', ' ').title()} match to {match_accession}"
                                    )
                    except Exception as e:
                        logger.error(
                            f"Error looking up metadata for {match_accession}: {e}"
                        )

                    # Update the sequence results with classification data
                    sequence_results["classification"] = classification_data
                    blast_results_store["sequence_results"]["0"] = sequence_results
                    logger.info(
                        f"Updated sequence results with classification: {classification_data}"
                    )

            # Initialize stages if not present
            if "stages" not in workflow_state:
                workflow_state["stages"] = {
                    stage["id"]: {"progress": 0, "complete": False}
                    for stage in WORKFLOW_STAGES
                }

            return workflow_state, False, blast_results_store

        # Safety check - should never actually get here due to the early returns above
        raise PreventUpdate

    except Exception as e:
        logger.error(f"Error updating workflow state: {e}")
        # Create a copy to ensure we return a new object
        if not isinstance(workflow_state, dict):
            workflow_state = {}
        else:
            workflow_state = workflow_state.copy()
        workflow_state["error"] = str(e)
        workflow_state["status"] = "failed"
        workflow_state["complete"] = True
        return workflow_state, False, blast_results_store


@callback(
    [
        Output("classification-progress", "value"),
        Output("classification-stage-display", "children"),
        Output("classification-progress-section", "style"),
    ],
    Input("workflow-state-store", "data"),
    prevent_initial_call=True,
    id="update-classification-progress-callback",
)
def update_classification_progress(workflow_state):
    """Update the classification progress UI based on workflow state"""
    if not workflow_state or not isinstance(workflow_state, dict):
        logger.warning("Invalid workflow state type or empty")
        return 0, "No workflow data", {"display": "none"}

    # Only show if we have a valid state with a status
    status = workflow_state.get("status", "")
    if not status or status == "pending" or "task_id" not in workflow_state:
        return 0, "", {"display": "none"}

    # Calculate progress percentage
    progress = 0
    if workflow_state.get("complete", False):
        progress = 100
    elif (
        "current_stage_idx" in workflow_state
        and workflow_state["current_stage_idx"] is not None
    ):
        try:
            stage_idx = int(workflow_state.get("current_stage_idx", 0))
            total_stages = len(WORKFLOW_STAGES)

            # Safely get the stage progress
            current_stage = workflow_state.get("current_stage")
            stages_dict = workflow_state.get("stages", {})

            if current_stage and current_stage in stages_dict:
                stage_data = stages_dict[current_stage]
                stage_progress = (
                    stage_data.get("progress", 0) if isinstance(stage_data, dict) else 0
                )
            else:
                stage_progress = 0

            # Calculate overall progress: stage contribution + progress within stage
            progress = int(
                (stage_idx / total_stages) * 100 + (stage_progress / total_stages)
            )
            # Ensure progress is within valid range
            progress = max(0, min(100, progress))
        except (ValueError, ZeroDivisionError, TypeError) as e:
            logger.error(f"Error calculating progress: {e}")
            progress = 0

    # Get current stage label
    if "error" in workflow_state and workflow_state["error"]:
        stage_text = f"Error: {workflow_state['error']}"
    elif workflow_state.get("complete", False):
        if workflow_state.get("found_match", False):
            match_stage = workflow_state.get("match_stage", "unknown")
            stage_text = f"Complete - {match_stage.capitalize()} match found"
        else:
            stage_text = "Classification complete"
    else:
        current_stage = workflow_state.get("current_stage")
        if current_stage:
            # Find the stage label
            stage_label = None
            for stage in WORKFLOW_STAGES:
                if stage["id"] == current_stage:
                    stage_label = stage["label"]
                    break

            stage_text = stage_label if stage_label else f"Processing {current_stage}"
        else:
            stage_text = "Starting classification..."

    # Hide the progress section if classification is complete
    if workflow_state.get("complete", False):
        style = {"display": "none"}
    else:
        # Show section if classification is in progress
        style = (
            {"display": "block"}
            if status and status != "pending"
            else {"display": "none"}
        )

    logger.debug(
        f"Progress update: progress={progress}, text='{stage_text}', style={style}"
    )
    return progress, stage_text, style


@callback(
    Output("classification-interval", "disabled", allow_duplicate=True),
    [
        Input("workflow-state-store", "data"),
        Input("blast-results-store", "data"),
    ],
    prevent_initial_call=True,
)
def disable_interval_when_complete(workflow_state, blast_results):
    """Disable the interval when classification is complete"""
    # Always disable the interval since we're no longer using polling
    return True

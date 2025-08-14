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
    run_classification_workflow_task,
)

from src.config.logging import get_logger

from src.utils.blast_data import BlastData, WorkflowState, FetchShipParams, FetchCaptainParams, ClassificationData

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
        dcc.Store(id="upload-sequences-store"),
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
        dcc.Store(id="blast-data-store"),
        dcc.Store(id="classification-data-store", data=None),
        dcc.Store(id="workflow-state-store", data=None),
        # Single data store with all workflow state
        dcc.Store(id="classification-stage", data="Upload a sequence"),
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
        Output("upload-sequences-store", "data", allow_duplicate=True),
        Output("upload-details", "children", allow_duplicate=True),
        Output("upload-error-message", "children", allow_duplicate=True),
    ],
    [
        Input("blast-fasta-upload", "contents"),
        Input("blast-fasta-upload", "filename"),
    ],
    prevent_initial_call=True,
)
def update_file_details(seq_content, seq_filename):
    """
    Immediately process the file to show a summary.
    Will display an error related to file problems in the upload area.
    - show and error and block submission if:
        - file is too large (over 10 MB)
        - file is not a valid FASTA format
        - file contains no sequences
        - file parsing fails for any reason
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

        if file_size > max_size:
            error_msg = f"The file '{seq_filename}' exceeds the 10 MB limit."
            error_alert = dmc.Alert(
                title="File Too Large",
                children=error_msg,
                color="red",
                variant="filled",
            )
            return True, "Error", None, upload_details, error_alert

        input_type, seq_list, n_seqs, error = check_input(
            query_text_input=None, query_file_contents=seq_content
        )

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
    [Input("blast-data-store", "data"), Input("blast-active-tab", "data")],
)
def process_blast_results(blast_results_dict, active_tab_idx):
    if not blast_results_dict:
        logger.warning("No blast_results_dict data")
        return None

    blast_results = BlastData.from_dict(blast_results_dict) if blast_results_dict else None

    # Convert active_tab_idx to string (since keys in sequence_results are strings)
    tab_idx = str(active_tab_idx or 0)
    logger.debug(f"Processing BLAST results for tab index: {tab_idx}")

    # Get the correct sequence results based on active tab
    sequence_results = blast_results.sequence_results.get(tab_idx)
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
        Output("upload-sequences-store", "data", allow_duplicate=True),
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
        State("upload-sequences-store", "data"),
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
    Processes input when the submit button is pressed.
    - For file uploads: Use the already parsed sequences from `upload-sequences-store`
    - For text input: Parse the input text now, limiting to max 10 sequences

    If we have parsed sequences from a file upload, use those
    Ensure limit to 10 sequences
    If we have file contents but no parsed sequences, try to parse them now
    This handles cases where the file upload callback hasn't completed before submit


    Input: 
        - n_clicks: Number of times the submit button has been clicked
        - query_text_input: Text input from the user
        - seq_list: Pre-parsed sequences from upload-sequences-store
        - file_contents: Contents of the uploaded FASTA file
    Output:
        - upload-sequences-store: Updated sequences data
        - upload-error-store: Error message if any
        - upload-error-message: Error message to display in the UI
        - right-column-content: Updated UI content for results display
        - submission-id-store: Unique ID for this submission to track processing

    """
    if not n_clicks:
        raise PreventUpdate

    error_alert = None
    ui_content = None 

    try:
        logger.debug(
            f"Preprocess called with seq_list={seq_list is not None}, query_text_input={query_text_input is not None}, file_contents={file_contents is not None}"
        )

        if seq_list is not None:
            logger.debug(
                f"Using pre-parsed sequences from upload-sequences-store: {len(seq_list)} sequences"
            )
            if len(seq_list) > 10:
                seq_list = seq_list[:10]

            ui_content = (
                create_tab_layout(seq_list)
                if len(seq_list) > 1
                else create_single_layout()
            )

            return seq_list, None, None, ui_content, n_clicks
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


def process_single_sequence(seq_data, evalue_threshold, curated=None):
    """Process a single sequence and return structured results"""

    # extract from seq_data
    query_header = seq_data.get("header", "query")
    query_seq = seq_data.get("sequence", "")
    query_type = seq_data.get("type", "nucl")
    query_length = len(query_seq)

    blast_data = BlastData(
        seq_type=query_type,
        sequence=query_seq,
    )

    tmp_query_fasta = write_temp_fasta(query_header, query_seq)
    if not tmp_query_fasta:
        error_message = "Failed to create temporary file"
        logger.error(error_message)
        blast_data.error = error_message

    meta_df = fetch_meta_data()

    try:
        blast_results = run_blast_search_task(
            query_header=query_header,
            query_seq=query_seq,
            query_type=query_type,
            eval_threshold=evalue_threshold,
            curated=curated,
        )

        # Handle cases where blast_results is:
        # -  None or invalid
        # - raw content
        # - a file path
        # - anything unexpected, i.e. if it's a string, assume it's content

        if not blast_results:
            logger.warning("BLAST search returned no results")
        elif isinstance(blast_results, dict) and "content" in blast_results:
            blast_content = blast_results["content"]
            blast_results_file = tempfile.NamedTemporaryFile(
                suffix=".blast", delete=False
            ).name
            with open(blast_results_file, "w") as f:
                f.write(blast_content)
            logger.debug(f"Created temporary BLAST results file: {blast_results_file}")
        elif isinstance(blast_results, str) and os.path.exists(blast_results):
            blast_results_file = blast_results
            try:
                with open(blast_results_file, "r") as f:
                    blast_content = f.read()
            except Exception as e:
                logger.error(f"Error reading blast file content: {e}")
            logger.debug(f"Using existing BLAST results file: {blast_results_file}")
        else:
            
            logger.error(f"Unexpected BLAST results format: {type(blast_results)}")
            if isinstance(blast_results, str):
                blast_content = blast_results
                blast_results_file = None
                logger.debug("Using blast results string as content")
            else:
                blast_results_file = None
                blast_content = None

        if blast_results_file and os.path.exists(blast_results_file):
            blast_tsv = parse_blast_xml(blast_results_file)
            if blast_tsv and os.path.exists(blast_tsv):
                blast_df = pd.read_csv(blast_tsv, sep="\t")
                if len(blast_df) == 0:
                    error_message = "No BLAST hits found"
                    logger.debug(error_message)

        # append to blast_data
        blast_data.blast_file = blast_results_file
        blast_data.blast_content = blast_content
        blast_data.fasta_file = tmp_query_fasta
        blast_data.blast_df = (
            blast_df.to_dict("records")
            if isinstance(blast_df, pd.DataFrame)
            else blast_df
        )

        # Classification will be handled by the workflow system

        blast_data.processed = True

        # Return structured results - convert to dict for backward compatibility
        return blast_data.to_dict()

    except Exception as e:
        blast_data.error = str(e)
        logger.error(f"Error in process_single_sequence: {e}")
        return blast_data.to_dict()


@callback(
    [
        Output("workflow-state-store", "data"),
        Output("classification-data-store", "data"),
        Output("classification-interval", "disabled", allow_duplicate=True),
        Output("blast-data-store", "data", allow_duplicate=True),
        Output("submit-button", "loading", allow_duplicate=True),
        Output("results-loading-overlay", "visible", allow_duplicate=True),
    ],
    [
        Input("submission-id-store", "data"),
    ],
    [
        State("upload-sequences-store", "data"),
        State("evalue-threshold", "value"),
        State("blast-fasta-upload", "contents"),
        State("curated-input", "value"),
    ],
    running=[
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
def process_multiple_sequences(
    submission_id, seq_list, evalue_threshold, file_contents, curated
):
    """Process only the first sequence initially, for both text input and file uploads"""

    workflow_state = WorkflowState(
        stages = {
            stage["id"]: {"progress": 0, "complete": False}
            for stage in WORKFLOW_STAGES
        },
        workflow_started = False
    )

    if not submission_id:
        logger.warning("No submission_id provided to process_multiple_sequences")
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
                    workflow_state.error = "Failed to parse sequence data"
                    workflow_state.complete = True
                    workflow_state.status = "failed"
                    workflow_state.found_match = False
                
                    return (
                        workflow_state.to_dict(),
                        None,
                        True,
                        None,
                        False,
                        False,
                    )

                logger.debug(
                    f"Successfully parsed {len(direct_seq_list)} sequences directly from file contents"
                )
                seq_list = direct_seq_list
            except Exception as e:
                logger.error(f"Error parsing file contents: {e}")
                workflow_state.error = "Failed to parse sequence data"
                workflow_state.complete = True
                workflow_state.status = "failed"
                workflow_state.found_match = False

                return (
                    workflow_state.to_dict(),
                    None,
                    True,
                    None,
                    False,
                    False,
                )
        else:
            logger.warning("No seq_list or file_contents provided to process_multiple_sequences")

            workflow_state.error = "No sequence data available"
            workflow_state.complete = True
            workflow_state.status = "failed"
            workflow_state.found_match = False

            return (
                workflow_state.to_dict(),
                None,
                True,
                None,
                False,
                False,
            )

    logger.debug(
        f"Processing sequence submission with ID: {submission_id}, sequences: {len(seq_list)}"
    )

    try:
        classification_interval_disabled = True  # Always disable the interval

        # Process only the first sequence to start
        if seq_list and len(seq_list) > 0:
            # Log sequence details for debugging
            first_seq = seq_list[0]
            logger.debug(
                f"Processing first sequence: header={first_seq.get('header', 'unknown')[:30]}..., length={len(first_seq.get('sequence', ''))}"
            )

            # Process sequence 0
            sequence_result = process_single_sequence(
                seq_list[0], evalue_threshold, curated
            )

            if sequence_result:
                logger.debug(
                    f"Sequence result: processed={sequence_result.get('processed', False)}, blast_file={sequence_result.get('blast_file') is not None}"
                )

                # Structure the blast results properly
                blast_data = BlastData(processed_sequences = [0],
                    sequence_results = {
                        "0": sequence_result,
                    },
                    total_sequences = len(seq_list)
                )

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
                    # Write the file content to include in classification_data
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

                    if not skip_classification and tmp_query_fasta:
                        # Set up the workflow state with proper classification_data
                        workflow_state.task_id = str(submission_id)
                        workflow_state.fetch_ship_params = FetchShipParams(
                            curated = False,
                            with_sequence = True, 
                            dereplicate = True
                        )
                        workflow_state.fetch_captain_params = FetchCaptainParams(
                            curated = True,
                            with_sequence = True
                        )
                        classification_data = ClassificationData(
                            seq_type = seq_list[0].get("type", "nucl"),
                            fasta_file = tmp_query_fasta
                        )
                        blast_data.seq_type = seq_list[0].get("type", "nucl")
                        blast_data.fasta_file = tmp_query_fasta

                    # Since we're no longer using Celery, we don't need the interval
                    # always set to disabled since we run synchronously
                    classification_interval_disabled = True

            else:
                logger.warning(
                    "No sequence result returned from process_single_sequence"
                )
                blast_data.processed_sequences = [0]
                blast_data.sequence_results = {
                    "0": {
                        "processed": False,
                        "error": "Failed to process sequence",
                        "sequence": seq_list[0].get("sequence", ""),
                    }
                }
                blast_data.total_sequences = len(seq_list)

        logger.debug(
            f"Completed process_multiple_sequences: has_blast_results={blast_data.blast_content is not None}, interval_disabled={classification_interval_disabled}"
        )
        return (
            workflow_state.to_dict(),
            classification_data.to_dict() if 'classification_data' in locals() and classification_data else None,
            classification_interval_disabled,
            blast_data.to_dict(),
            False,
            False,
        )  # Set loading to False when done
    except Exception as e:
        logger.error(f"Error in process_multiple_sequences: {str(e)}")
        # Return basic data on error
        workflow_state.complete = True
        workflow_state.error = str(e)
        workflow_state.status = "failed"
        workflow_state.task_id = str(submission_id) if submission_id else None
        
        return (
            workflow_state.to_dict(),
            None,
            True,
            None,
            False,
            False,
        )  # Set loading to False on error too


@callback(
    Output("blast-data-store", "data", allow_duplicate=True),
    Input("blast-tabs", "active_tab"),
    [
        State("blast-data-store", "data"),
        State("upload-sequences-store", "data"),
        State("evalue-threshold", "value"),
        State("curated-input", "value"),
    ],
    prevent_initial_call=True,
)
def process_additional_sequence(
    active_tab, results_store, seq_list, evalue_threshold, curated
):
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
        sequence_result = process_single_sequence(
            seq_list[tab_idx], evalue_threshold, curated
        )

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
                # Create temporary FASTA file for the workflow
                tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
                with open(tmp_fasta, "w") as f:
                    f.write(
                        f">{seq_list[tab_idx].get('header', 'query')}\n{seq_list[tab_idx].get('sequence', '')}\n"
                    )

                # Set up workflow data
                workflow_state = WorkflowState(
                    stages={
                        stage["id"]: {"progress": 0, "status": "pending"}
                        for stage in WORKFLOW_STAGES
                    }
                )
                # set up workflow state parameters
                workflow_state.fetch_ship_params = FetchShipParams(
                    curated=False, with_sequence=True, dereplicate=True
                )
                workflow_state.fetch_captain_params = FetchCaptainParams(
                    curated=True, with_sequence=True
                )
                
                # set up blast data
                blast_data = BlastData(
                    seq_type=seq_list[tab_idx].get("type", "nucl"),
                    fasta_file=tmp_fasta,
                    blast_df=sequence_result.get("blast_df"),
                )

                meta_df = fetch_meta_data()
                meta_dict = meta_df.to_dict("records") if meta_df is not None else None

                # Run classification workflow
                workflow_result = run_classification_workflow_task(
                    workflow_state_dict=workflow_state, blast_data_dict=blast_data, meta_dict=meta_dict
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
        Input("blast-data-store", "data"),
    ],
    prevent_initial_call=True,
)
def render_tab_content(active_tab, blast_data_dict):
    """Render the content for the current tab"""
    from src.utils.classification_utils import create_classification_output

    blast_data = BlastData.from_dict(blast_data_dict) if blast_data_dict else None

    try:
        if not active_tab or not blast_data:
            logger.warning(
                f"Missing data for tab rendering: active_tab={active_tab}, blast_data={blast_data is not None}"
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
        processed_sequences = blast_data.processed_sequences
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
        sequence_results = blast_data.sequence_results.get(str(tab_idx))
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
        classification_data = sequence_results.get("classification") if sequence_results else None
        classification_output = create_classification_output(workflow_state=None, classification_data=classification_data)

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
                        
                        // BlasterJS should now display clean accession numbers
                        const buttons = container.querySelectorAll('.alignment-table-description');
                        console.log(`Found ${buttons.length} BlasterJS buttons in ${containerId}`);
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
                              
                              // BlasterJS should now display clean accession numbers
                              const buttons = container.querySelectorAll('.alignment-table-description');
                              console.log(`Found ${buttons.length} BlasterJS buttons in tab ${tabIdx}`);
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
        Input("blast-data-store", "data"),
        Input("workflow-state-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_single_sequence_classification(blast_results_dict, workflow_state_dict):
    """Update classification output for single sequence submissions"""
    from src.utils.classification_utils import create_classification_output

    blast_results = BlastData.from_dict(blast_results_dict) if blast_results_dict else None
    workflow_state = WorkflowState.from_dict(workflow_state_dict) if workflow_state_dict else None

    logger.info(
        f"update_single_sequence_classification CALLED: blast_results_dict={blast_results is not None}, workflow_state={workflow_state is not None}"
    )

    if workflow_state:
        logger.info(
            f"Workflow state: complete={workflow_state.complete}, found_match={workflow_state.found_match}, match_stage={workflow_state.match_stage}, match_result={workflow_state.match_result}"
        )

    if not blast_results:
        logger.info("No blast_results, returning None")
        return None

    # Get the results for sequence 0
    sequence_results = blast_results.sequence_results.get("0")
    if not sequence_results:
        logger.info("No sequence results for sequence 0, returning None")
        return None

    logger.info(
        f"Original sequence_results has classification: {'classification' in sequence_results}"
    )

    # Check if we have workflow results that should override the initial classification
    if (
        workflow_state
        and workflow_state.complete
        and workflow_state.found_match
    ):
        logger.info("Workflow completed with match found, updating classification")
        # Update sequence_results with workflow classification data
        classification_data = workflow_state.classification_data
        if classification_data:
            # Use workflow classification data directly
            sequence_results = sequence_results.copy()
            sequence_results["classification"] = classification_data
            logger.info("Used workflow classification_data as classification")
            logger.debug(f"Workflow classification_data: {classification_data}")
        elif workflow_state.match_stage and workflow_state.match_result:
            # Convert simple workflow results to classification format and look up metadata
            sequence_results = sequence_results.copy()
            match_stage = workflow_state.match_stage
            match_accession = workflow_state.match_result

            logger.info(
                f"Processing workflow result: {match_stage} -> {match_accession}"
            )

            # Create base classification
            classification_data = ClassificationData(
                source = match_stage,
                closest_match = match_accession,
                confidence = "High"
                if match_stage in ["exact", "contained"]
                else "Medium"
            )

            # Look up additional metadata for the matched accession
            try:
                meta_df = fetch_meta_data()
                if meta_df is not None and not meta_df.empty:
                    meta_match = meta_df[
                        meta_df["accession_display"] == match_accession
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
                    else:
                        logger.warning(f"No metadata found for {match_accession}")

            except Exception as e:
                logger.error(f"Error looking up metadata for {match_accession}: {e}")

            sequence_results["classification"] = classification_data.to_dict()
            logger.info(f"Created classification_data: {classification_data}")
    else:
        logger.info(
            f"Workflow not complete or no match found. Complete: {workflow_state.complete if workflow_state else 'N/A'}, Found match: {workflow_state.found_match if workflow_state else 'N/A'}"
        )

    # Create the classification output using the existing function
    logger.info(
        f"Final sequence_results has classification: {'classification' in sequence_results}"
    )
    
    # Get classification data from workflow state or sequence results
    # Prefer sequence results since they contain metadata lookup results
    classification_data = None
    if sequence_results and "classification" in sequence_results:
        classification_data = sequence_results["classification"]
        logger.debug(f"Using classification data from sequence_results: {classification_data}")
    elif workflow_state and workflow_state.classification_data:
        classification_data = workflow_state.classification_data
        logger.debug(f"Using classification data from workflow_state: {classification_data}")
    else:
        logger.debug("No classification data found in either source")
    
    result = create_classification_output(workflow_state=workflow_state, classification_data=classification_data)
    logger.info(f"create_classification_output returned: {type(result)}")
    return result


@callback(
    [
        Output("workflow-state-store", "data", allow_duplicate=True),
        Output("results-loading-overlay", "visible", allow_duplicate=True),
        Output("blast-data-store", "data", allow_duplicate=True),
    ],
    [
        # Remove dependency on interval n_intervals
        Input("workflow-state-store", "data"),
        Input("classification-data-store", "data"),
    ],
    [
        State("blast-data-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_classification_workflow_state(workflow_state_dict, classification_data_dict, blast_results_dict):
    """Initialize workflow state (one-time operation)"""
    # If no workflow state, nothing to do
    if workflow_state_dict is None or classification_data_dict is None:
        raise PreventUpdate

    workflow_state = WorkflowState.from_dict(workflow_state_dict) if workflow_state_dict else None
    classification_data = ClassificationData.from_dict(classification_data_dict) if classification_data_dict else None
    blast_results = BlastData.from_dict(blast_results_dict) if blast_results_dict else None

    # Return early if workflow is complete, has error, or already started
    if (
        workflow_state.complete
        or workflow_state.error is not None
        or workflow_state.workflow_started
    ):
        # If workflow is complete, also set loading overlay to false
        if (
            workflow_state.complete
            or workflow_state.error is not None
        ):
            return workflow_state.to_dict(), False, blast_results.to_dict()
        raise PreventUpdate

    # Check if task is running
    task_id = workflow_state.task_id
    if not task_id:
        logger.warning("No task ID found in workflow state")
        workflow_state.error = "Missing task ID"
        workflow_state.complete = True
        workflow_state.status = "failed"
        return workflow_state.to_dict(), False, blast_results.to_dict()

    try:
        # Start the classification workflow - ONE TIME ONLY
        if workflow_state.task_id and not workflow_state.workflow_started:
            # Extract BLAST data from the blast_results_store
            blast_df = None
            if blast_results and blast_results.sequence_results:
                sequence_results = blast_results.sequence_results.get("0")
                if sequence_results and "blast_df" in sequence_results:
                    blast_df = sequence_results["blast_df"]

            blast_data = BlastData(
                seq_type=classification_data.seq_type,
                fasta_file=classification_data.fasta_file,
                blast_df=blast_df,  # Pass BLAST results to the workflow
            )

            meta_df = fetch_meta_data()
            meta_dict = meta_df.to_dict("records") if meta_df is not None else None

            logger.debug("Running classification workflow directly (no Celery)")
            # Run the workflow function directly
            result = run_classification_workflow_task(
                workflow_state=workflow_state.to_dict(), 
                blast_data=blast_data.to_dict(), 
                classification_data=classification_data.to_dict(),
                meta_dict=meta_dict
            )

            logger.debug(f"Workflow result type: {type(result)}, content: {result}")

            # Convert result back to workflow state object
            if isinstance(result, dict):
                workflow_state = WorkflowState.from_dict(result)
            else:
                workflow_state = result

            # Mark as started so we don't run it again
            workflow_state.workflow_started = True

            # The workflow_state should already be updated by the workflow function
            # Just ensure the status is set properly
            workflow_state.status = "complete"
            workflow_state.complete = True

            logger.debug(
                f"Final workflow state: complete={workflow_state.complete}, found_match={workflow_state.found_match}, match_stage={workflow_state.match_stage}, match_result={workflow_state.match_result}"
            )

            # IMPORTANT: Also update the sequence results with classification data for tabbed interface
            if (
                workflow_state.complete
                and workflow_state.found_match
                and blast_results
                and blast_results.sequence_results
            ):
                logger.info(
                    "Updating sequence results with workflow classification for tabbed interface"
                )

                # Get the sequence results for sequence 0 - we'll modify this directly
                sequence_results = blast_results.sequence_results.get("0")
                if sequence_results:
                    # Make a copy of sequence_results to avoid modifying the original
                    sequence_results = sequence_results.copy()

                    match_stage = workflow_state.match_stage
                    match_accession = workflow_state.match_result

                    # Create classification data
                    classification_data = ClassificationData(
                        source=match_stage,
                        closest_match=match_accession,
                        confidence="High"
                        if match_stage in ["exact", "contained"]
                        else "Medium"
                    )

                    # Look up metadata for the matched accession
                    try:
                        meta_df = fetch_meta_data()
                        if meta_df is not None and not meta_df.empty:
                            meta_match = meta_df[
                                meta_df["accession_display"] == match_accession
                            ]
                            if not meta_match.empty:
                                logger.debug(f"Found metadata for {match_accession}")
                                # Add family information
                                if "familyName" in meta_match.columns:
                                    family = meta_match["familyName"].iloc[0]
                                    if family and family != "None":
                                        classification_data.family = family

                                # Add navis information
                                if "navis_name" in meta_match.columns:
                                    navis = meta_match["navis_name"].iloc[0]
                                    if navis and navis != "None":
                                        classification_data.navis = navis

                                # Add haplotype information
                                if "haplotype_name" in meta_match.columns:
                                    haplotype = meta_match["haplotype_name"].iloc[0]
                                    if haplotype and haplotype != "None":
                                        classification_data.haplotype = haplotype

                                # Add match details
                                if match_stage == "exact":
                                    classification_data.match_details = (
                                        f"Exact sequence match to {match_accession}"
                                    )
                                elif match_stage == "contained":
                                    classification_data.match_details = (
                                        f"Query sequence contained within {match_accession}"
                                    )
                                elif match_stage == "similar":
                                    classification_data.match_details = (
                                        f"High similarity to {match_accession}"
                                    )
                                else:
                                    classification_data.match_details = (
                                        f"{match_stage.replace('_', ' ').title()} match to {match_accession}"
                                    )
                    except Exception as e:
                        logger.error(
                            f"Error looking up metadata for {match_accession}: {e}"
                        )

                    # Update the sequence results with classification data
                    sequence_results["classification"] = classification_data.to_dict()
                    blast_results.sequence_results["0"] = sequence_results
                    
                    # Also update the workflow state's classification_data with the enriched data
                    workflow_state.classification_data = classification_data.to_dict()
                    
                    logger.info(
                        f"Updated sequence results with classification: {classification_data}"
                    )

            # Initialize stages if not present
            if not workflow_state.stages:
                workflow_state.stages = {
                    stage["id"]: {"progress": 0, "complete": False}
                    for stage in WORKFLOW_STAGES
                }

            return workflow_state.to_dict(), False, blast_results.to_dict()

        # Safety check - should never actually get here due to the early returns above
        raise PreventUpdate

    except Exception as e:
        logger.error(f"Error updating workflow state: {e}")
        # Update the workflow state with error information
        workflow_state.error = str(e)
        workflow_state.status = "failed"
        workflow_state.complete = True
        return workflow_state.to_dict(), False, blast_results.to_dict()


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
def update_classification_progress(workflow_state_dict):
    """Update the classification progress UI based on workflow state"""
    if not workflow_state_dict or not isinstance(workflow_state_dict, dict):
        logger.warning("Invalid workflow state type or empty")
        return 0, "No workflow data", {"display": "none"}

    workflow_state = WorkflowState.from_dict(workflow_state_dict)

    # Only show if we have a valid state with a status
    status = workflow_state.status
    if not status or status == "pending" or not workflow_state.task_id:
        return 0, "", {"display": "none"}

    # Calculate progress percentage
    progress = 0
    if workflow_state.complete:
        progress = 100
    elif (
        workflow_state.current_stage_idx is not None
        and workflow_state.current_stage_idx is not None
    ):
        try:
            stage_idx = int(workflow_state.current_stage_idx)
            total_stages = len(WORKFLOW_STAGES)

            # Safely get the stage progress
            current_stage = workflow_state.current_stage
            stages_dict = workflow_state.stages

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
    if workflow_state.error:
        stage_text = f"Error: {workflow_state.error}"
    elif workflow_state.complete:
        if workflow_state.found_match:
            match_stage = workflow_state.match_stage or "unknown"
            stage_text = f"Complete - {match_stage.capitalize()} match found"
        else:
            stage_text = "Classification complete"
    else:
        current_stage = workflow_state.current_stage
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
    if workflow_state.complete:
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
        Input("blast-data-store", "data"),
    ],
    prevent_initial_call=True,
)
def disable_interval_when_complete(workflow_state, blast_results):
    """Disable the interval when classification is complete"""
    # Always disable the interval since we're no longer using polling
    return True

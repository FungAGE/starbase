import dash
from dash import dcc, html, callback, clientside_callback
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import os
import base64
import pandas as pd

from src.config.cache import cache
from src.utils.seq_utils import (
    check_input,
    write_temp_fasta,
    seq_processing_error_alert,
)
from src.utils.blast_utils import create_no_matches_alert


from src.components.ui import (
    curated_switch,
    create_file_upload,
)
from src.database.sql_manager import fetch_meta_data
from src.utils.classification_utils import WORKFLOW_STAGES
from src.tasks import (
    run_blast_search_task,
    run_classification_workflow_sync,
)

from src.config.logging import get_logger

from src.utils.blast_data import (
    get_dash_adapter,
    BlastData,
    WorkflowState,
    FetchShipParams,
    FetchCaptainParams,
    ClassificationData,
    safe_convert_sequence_analysis_to_legacy,
)

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
                                        [
                                            "Search protein/nucleotide sequences for ",
                                            html.Span(
                                                "Starships",
                                                style={"fontStyle": "italic"},
                                            ),
                                            " and ",
                                            html.Span(
                                                "Starship-associated genes",
                                                style={"fontStyle": "italic"},
                                            ),
                                        ],
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
                                                                                text=html.Div(
                                                                                    [
                                                                                        "Only search curated ",
                                                                                        html.Span(
                                                                                            "Starships",
                                                                                            style={
                                                                                                "fontStyle": "italic"
                                                                                            },
                                                                                        ),
                                                                                    ]
                                                                                ),
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
                dmc.Title(
                    "BLAST Results",
                    order=2,
                    style={"marginTop": "15px", "marginBottom": "20px"},
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
    """
    Process the BLAST results for the active tab.
    - If we have direct blast_content, use it without reading file
    - If no blast_file, return empty blast text with warning
    - If we have blast_file, read it and pass the raw BLAST text directly for BlasterJS
    - If there's an error, return empty blast text
    Note: there is a 5MB limit on the size of the BLAST output.
    """
    if not blast_results_dict:
        logger.warning("No blast_results_dict data")
        return None

    blast_results = (
        BlastData.from_dict(blast_results_dict) if blast_results_dict else None
    )

    if (
        len(blast_results.sequence_results) == 1
        and "0" in blast_results.sequence_results
    ):
        tab_idx = "0"
        logger.debug(f"Single sequence detected, using tab index: {tab_idx}")
    else:
        tab_idx = str(active_tab_idx or 0)
        logger.debug(f"Processing BLAST results for tab index: {tab_idx}")

    sequence_results = blast_results.sequence_results.get(tab_idx)
    if not sequence_results:
        # Fallback: try to get the first available sequence result
        if blast_results.sequence_results:
            first_key = next(iter(blast_results.sequence_results.keys()))
            sequence_results = blast_results.sequence_results[first_key]
            logger.warning(
                f"No sequence results for tab index {tab_idx}, using first available: {first_key}"
            )
        else:
            logger.warning("No sequence results available at all")
            return None

    blast_results_file = sequence_results.get("blast_file")
    blast_content = sequence_results.get("blast_content")

    if blast_content:
        logger.debug(f"Using direct blast_content for tab {tab_idx}")
        return {"blast_text": blast_content}

    if not blast_results_file:
        logger.warning(f"No blast file in sequence results for tab {tab_idx}")
        return {"blast_text": ""}

    logger.debug(f"Reading BLAST file: {blast_results_file}")
    try:
        if not os.path.exists(blast_results_file):
            logger.error(f"BLAST file does not exist: {blast_results_file}")
            return {"blast_text": ""}

        with open(blast_results_file, "r") as f:
            blast_results = f.read()

        results_size = len(blast_results)
        logger.debug(f"Read BLAST results, size: {results_size} bytes")
        if results_size > 5 * 1024 * 1024:
            logger.warning(f"BLAST results too large: {results_size} bytes")
            return {"blast_text": "BLAST results too large to display"}

        data = {"blast_text": blast_results}

        return data
    except Exception as e:
        logger.error(f"Error processing BLAST results: {str(e)}")
        return {"blast_text": ""}


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
        State("blast-fasta-upload", "contents"),
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
                    if len(seq_list) > 10:
                        seq_list = seq_list[:10]

                    ui_content = (
                        create_tab_layout(seq_list)
                        if len(seq_list) > 1
                        else create_single_layout()
                    )

                    return seq_list, None, None, ui_content, n_clicks
            except Exception as e:
                logger.error(f"Error processing file contents directly: {e}")

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

            ui_content = (
                create_tab_layout(seq_list)
                if len(seq_list) > 1
                else create_single_layout()
            )

            return seq_list, None, None, ui_content, n_clicks

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


def process_single_sequence(seq_data, evalue_threshold, curated=None, sequence_id=None):
    """Process a single sequence and return structured results
    Handle cases where blast_results is:
        - None or invalid
        - raw content (i.e. if it's a string, assume it's content)
        - a file path
        - anything unexpected
    Classification will be handled by the workflow system.
    Returns:
        - structured results
        - converted to dict for backward compatibility

    - Parse BLAST content if available
    - Parse the XML file to TSV format
    - Read the parsed TSV file
    - Attach BLAST result to analysis
    - Mark as complete if no errors
    """

    from src.utils.blast_data import (
        SequenceAnalysis,
        BlastResult,
        SequenceType,
        WorkflowConfig,
        WorkflowStatus,
    )
    from src.utils.blast_utils import parse_blast_xml

    # Extract basic sequence information
    query_header = seq_data.get("header", "query")
    query_seq = seq_data.get("sequence", "")
    query_type = seq_data.get("type", "nucl")

    # Generate sequence ID if not provided
    if not sequence_id:
        import time

        sequence_id = f"seq_{int(time.time() * 1000)}"

    # Create the unified SequenceAnalysis object
    analysis = SequenceAnalysis(
        sequence_id=sequence_id,
        sequence=query_seq,
        sequence_header=query_header,
        sequence_type=SequenceType.NUCLEOTIDE
        if query_type == "nucl"
        else SequenceType.PROTEIN,
        status=WorkflowStatus.RUNNING,
    )

    # Set configuration based on parameters
    analysis.config = WorkflowConfig(
        ship_curated_only=curated or False,
        captain_curated_only=curated or False,
    )

    logger.info(f"Processing sequence {sequence_id}: {query_header[:50]}...")

    try:
        # Create temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        if not tmp_query_fasta:
            error_message = "Failed to create temporary file"
            logger.error(error_message)
            analysis.set_error(error_message)
            return analysis

        blast_results = run_blast_search_task(
            query_header=query_header,
            query_seq=query_seq,
            query_type=query_type,
            eval_threshold=evalue_threshold,
            curated=curated,
        )

        blast_result = BlastResult(
            sequence=query_seq,
            sequence_type=analysis.sequence_type,
            fasta_file=tmp_query_fasta,
        )

        if not blast_results:
            logger.warning("BLAST search returned no results")
            blast_result.error = "No BLAST results returned"
        elif isinstance(blast_results, dict) and "content" in blast_results:
            blast_result.blast_content = blast_results["content"]
            blast_result.blast_file = blast_results.get("file")
        elif isinstance(blast_results, str):
            import os

            if os.path.exists(blast_results):
                blast_result.blast_file = blast_results
                try:
                    with open(blast_results, "r") as f:
                        blast_result.blast_content = f.read()
                except Exception as e:
                    logger.error(f"Failed to read BLAST file {blast_results}: {e}")
                    blast_result.error = f"Failed to read BLAST file: {e}"
            else:
                blast_result.blast_content = blast_results

        if (
            blast_result.blast_file
            and os.path.exists(blast_result.blast_file)
            and not blast_result.error
        ):
            try:
                blast_tsv = parse_blast_xml(blast_result.blast_file)

                if blast_tsv and os.path.exists(blast_tsv):
                    blast_df = pd.read_csv(blast_tsv, sep="\t")

                    if len(blast_df) == 0:
                        logger.warning(f"No BLAST hits found for {sequence_id}")
                        blast_result.blast_hits = []
                    else:
                        blast_result.blast_hits = blast_df.to_dict("records")
                        logger.info(
                            f"Successfully processed {len(blast_df)} BLAST hits for {sequence_id}"
                        )

                    blast_result.processed = True
                else:
                    error_message = "Failed to parse BLAST XML file"
                    logger.error(error_message)
                    blast_result.error = error_message

            except Exception as e:
                error_message = f"Failed to parse BLAST output: {e}"
                logger.error(error_message)
                blast_result.error = error_message

        analysis.blast_result = blast_result

        if not blast_result.error:
            analysis.set_complete()
        else:
            analysis.set_error(blast_result.error)

        logger.info(f"Completed processing sequence {sequence_id}")
        return analysis

    except Exception as e:
        error_message = f"Error in process_single_sequence: {e}"
        logger.error(error_message)
        analysis.set_error(error_message)
        return analysis


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
    """
    Process multiple sequences using centralized state.
    Handle missing seq_list by parsing file_contents if available.
    Performs the following steps:
        - Add sequence to centralized state
        - Process the first sequence using existing logic
        - Import the existing process_single_sequence function
        - Update centralized state with BLAST data
        - Set up workflow state
        - Determine if we should run classification
        - Update centralized state
    Returns:
        - data for Dash stores (for backward compatibility)
    """

    if not submission_id:
        logger.warning("No submission_id provided to process_multiple_sequences")
        raise PreventUpdate

    adapter = get_dash_adapter()
    sequence_id = str(submission_id)

    try:
        if not seq_list and file_contents:
            logger.debug("No seq_list but file_contents available - parsing directly")
            try:
                input_type, direct_seq_list, n_seqs, error = check_input(
                    query_text_input=None, query_file_contents=file_contents
                )
                if error or not direct_seq_list or len(direct_seq_list) == 0:
                    logger.warning(f"Failed to parse file contents directly: {error}")
                    return _return_error_state(
                        adapter, sequence_id, "Failed to parse sequence data"
                    )
                seq_list = direct_seq_list
            except Exception as e:
                logger.error(f"Error parsing file contents: {e}")
                return _return_error_state(
                    adapter, sequence_id, "Failed to parse sequence data"
                )

        if not seq_list:
            logger.warning("No seq_list or file_contents provided")
            return _return_error_state(
                adapter, sequence_id, "No sequence data available"
            )

        logger.debug(
            f"Processing sequence submission with ID: {sequence_id}, sequences: {len(seq_list)}"
        )

        pipeline_state = adapter.pipeline_state
        # Start a new submission - this clears any old state
        sequence_state = pipeline_state.start_new_submission(sequence_id)

        first_seq = seq_list[0]
        logger.debug(
            f"Processing first sequence: header={first_seq.get('header', 'unknown')[:30]}..., length={len(first_seq.get('sequence', ''))}"
        )

        # Process the sequence with the correct sequence ID
        sequence_analysis = process_single_sequence(
            first_seq, evalue_threshold, curated, sequence_id
        )

        if not sequence_analysis:
            logger.warning("No sequence analysis returned from process_single_sequence")
            return _return_error_state(
                adapter, sequence_id, "Failed to process sequence"
            )

        # Convert SequenceAnalysis to legacy format for backward compatibility
        sequence_result = safe_convert_sequence_analysis_to_legacy(
            sequence_analysis, tab_idx=0
        )

        if not sequence_result:
            logger.warning("Failed to convert sequence analysis to legacy format")
            return _return_error_state(
                adapter, sequence_id, "Failed to convert sequence data"
            )

        # Create BlastData with the sequence result
        blast_data = BlastData(
            processed_sequences=[0],
            sequence_results={"0": sequence_result},
            total_sequences=len(seq_list),
        )

        # Update the centralized state with BLAST data using the submission ID
        pipeline_state.update_blast_data(sequence_id, blast_data)

        # Also store the sequence analysis directly in the centralized state
        # This ensures the data is available for the workflow
        if sequence_analysis.blast_result:
            # Update the sequence state with the blast result
            sequence_state = pipeline_state.get_sequence(sequence_id)
            if sequence_state:
                # Store the blast content and file in the blast data
                blast_data.blast_content = sequence_analysis.blast_result.blast_content
                blast_data.blast_file = sequence_analysis.blast_result.blast_file
                blast_data.fasta_file = sequence_analysis.blast_result.fasta_file
                blast_data.seq_type = sequence_analysis.sequence_type.value
                blast_data.processed = sequence_analysis.is_complete()

                # Update the blast data again with the complete information
                pipeline_state.update_blast_data(sequence_id, blast_data)

        # Create workflow state using the submission ID
        workflow_state = WorkflowState(
            stages={
                stage["id"]: {"progress": 0, "complete": False}
                for stage in WORKFLOW_STAGES
            },
            workflow_started=False,
            task_id=sequence_id,
        )

        # Ensure the sequence ID is properly set
        logger.debug(f"Created workflow state with task_id: {workflow_state.task_id}")

        sequence_length = len(sequence_analysis.sequence or "")
        skip_classification = (
            sequence_length < 5000
            or sequence_analysis.has_error()
            or not sequence_analysis.blast_result
            or not sequence_analysis.blast_result.blast_content
        )

        logger.debug(
            f"Classification decision: skip={skip_classification}, seq_length={sequence_length}"
        )

        classification_data = None
        if not skip_classification:
            workflow_state.fetch_ship_params = FetchShipParams(
                curated=False, with_sequence=True, dereplicate=True
            )
            workflow_state.fetch_captain_params = FetchCaptainParams(
                curated=True, with_sequence=True
            )

            # Use the FASTA file from the SequenceAnalysis
            tmp_query_fasta = sequence_analysis.blast_result.fasta_file
            if tmp_query_fasta:
                classification_data = ClassificationData(
                    seq_type=sequence_analysis.sequence_type.value,
                    fasta_file=tmp_query_fasta,
                )
                blast_data.seq_type = sequence_analysis.sequence_type.value
                blast_data.fasta_file = tmp_query_fasta
            else:
                skip_classification = True
                logger.error("No FASTA file available from sequence analysis")

        pipeline_state.update_workflow_state(sequence_id, workflow_state)
        if classification_data:
            pipeline_state.update_classification_data(sequence_id, classification_data)

        store_data = adapter.sync_all_stores(sequence_id)

        logger.debug(f"Completed unified sequence processing for {sequence_id}")
        return (
            store_data["workflow_state"],
            store_data["classification_data"],
            True,  # Always disable interval (no longer using polling)
            store_data["blast_data"],
            False,  # Stop loading
            False,  # Hide loading overlay
        )

    except Exception as e:
        logger.error(f"Error in process_multiple_sequences: {str(e)}")
        return _return_error_state(adapter, sequence_id, str(e))


def _return_error_state(adapter, sequence_id, error_message):
    """Helper function to return consistent error state"""
    pipeline_state = adapter.pipeline_state

    # Ensure sequence exists before setting error
    if sequence_id not in pipeline_state._sequences:
        pipeline_state.add_sequence(sequence_id)

    # Set error in centralized state
    pipeline_state.set_sequence_error(sequence_id, error_message)

    # Create error workflow state
    workflow_state = WorkflowState(
        stages={
            stage["id"]: {"progress": 0, "complete": False} for stage in WORKFLOW_STAGES
        },
        workflow_started=False,
        complete=True,
        status="failed",
        error=error_message,
        task_id=sequence_id,
    )
    pipeline_state.update_workflow_state(sequence_id, workflow_state)

    # Return error state for all stores
    store_data = adapter.sync_all_stores(sequence_id)
    return (
        store_data["workflow_state"],
        store_data["classification_data"],
        True,  # Disable interval
        store_data["blast_data"],
        False,  # Stop loading
        False,  # Hide loading overlay
    )


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
    """
    Process additional sequence using centralized state.

    Returns:
        - updated results store for backward compatibility

    Performs:
        - Creates unique sequence ID for each tab
        - Checks if each sequence has already been processed
        - Processes each sequence using existing logic
        - Uses the new unified approach for each tab
        - Checks if each sequence should get classification workflow (≥5000bp)
        - Adds each sequence to centralized state for classification
        - Sets up BLAST data in centralized state
        - Sets up workflow state
        - Sets up classification data
        - Runs classification workflow
        - If workflow found a match, updates centralized state
        - Creates enriched classification data using the helper function
        - Updates centralized state
        - Updates the results store for backward compatibility
    """
    try:
        if not active_tab or not results_store or not seq_list:
            logger.warning(
                f"Missing data for tab processing: active_tab={active_tab is not None}, results_store={results_store is not None}, seq_list={seq_list is not None}"
            )
            raise PreventUpdate

        try:
            tab_idx = (
                int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0
            )
        except (ValueError, IndexError) as e:
            logger.error(f"Invalid tab format '{active_tab}': {e}")
            raise PreventUpdate

        adapter = get_dash_adapter()
        pipeline_state = adapter.pipeline_state

        main_sequence_id = pipeline_state._active_sequence_id
        if not main_sequence_id:
            # Try to find the main sequence ID from existing sequences
            # Look for sequences that don't have "_tab_" in their ID
            main_sequences = [seq_id for seq_id in pipeline_state._sequences.keys() 
                            if "_tab_" not in seq_id]
            if main_sequences:
                main_sequence_id = main_sequences[0]
                logger.info(f"Recovered main sequence ID: {main_sequence_id}")
                pipeline_state.set_active_sequence(main_sequence_id)
            else:
                logger.error("No active sequence ID found for tab processing and no main sequences available")
                raise PreventUpdate
        tab_sequence_id = f"{main_sequence_id}_tab_{tab_idx}"

        processed_sequences = results_store.get("processed_sequences", [])
        logger.debug(
            f"Tab {tab_idx} clicked. Processed sequences: {processed_sequences}"
        )

        if tab_idx in processed_sequences:
            logger.debug(f"Tab {tab_idx} already processed, skipping")
            raise PreventUpdate

        if tab_idx >= len(seq_list):
            logger.error(
                f"Tab index {tab_idx} out of range (seq_list length: {len(seq_list)})"
            )
            raise PreventUpdate

        logger.info(
            f"Processing sequence for tab {tab_idx}: {seq_list[tab_idx].get('header', 'unknown')[:50]}..."
        )
        sequence_id = f"tab_{tab_idx}"

        # Process using the new unified approach
        logger.info(f"Processing sequence for tab {tab_idx} using unified approach")
        analysis = process_single_sequence(
            seq_list[tab_idx], evalue_threshold, curated, sequence_id
        )

        if not analysis:
            logger.error(
                f"Failed to process sequence for tab {tab_idx} - no analysis returned"
            )
            raise PreventUpdate

        # Convert to legacy format for backward compatibility
        sequence_result = safe_convert_sequence_analysis_to_legacy(
            analysis, tab_idx=tab_idx
        )

        if not sequence_result:
            logger.error(
                f"Failed to convert analysis to legacy format for tab {tab_idx}"
            )
            raise PreventUpdate

        logger.info(
            f"Successfully processed and converted tab {tab_idx} using unified approach"
        )

        if analysis.has_error():
            logger.error(
                f"Error processing sequence for tab {tab_idx}: {analysis.error}"
            )

        sequence_length = len(analysis.sequence or "")
        logger.debug(
            f"Analysis object has blast_result: {hasattr(analysis, 'blast_result')}"
        )
        logger.debug(
            f"BLAST result exists: {analysis.blast_result is not None if hasattr(analysis, 'blast_result') else 'N/A'}"
        )

        # Check for perfect BLAST matches (100% identity, full coverage)
        perfect_blast_match = None
        logger.debug(
            f"Checking for perfect BLAST matches: blast_result={analysis.blast_result is not None if hasattr(analysis, 'blast_result') else False}"
        )
        if analysis.blast_result:
            logger.debug(
                f"BLAST result has hits: {hasattr(analysis.blast_result, 'blast_hits')}"
            )
            if hasattr(analysis.blast_result, "blast_hits"):
                logger.debug(f"BLAST hits: {analysis.blast_result.blast_hits}")

        if (
            analysis.blast_result
            and analysis.blast_result.blast_hits
            and len(analysis.blast_result.blast_hits) > 0
        ):
            logger.debug(
                f"Checking {len(analysis.blast_result.blast_hits)} BLAST hits for perfect matches"
            )
            logger.debug(f"Sequence length: {sequence_length}")

            # Look for a perfect match (100% identity, full query coverage)
            for i, hit in enumerate(analysis.blast_result.blast_hits):
                logger.debug(f"BLAST hit {i}: {hit}")

                pident = hit.get("pident", 0)
                query_start = hit.get("query_start", 0)  # 1-based
                query_end = hit.get("query_end", 0)  # 1-based
                aln_length = hit.get("aln_length", 0)
                hit_ids = hit.get("hit_IDs")

                logger.debug(
                    f"Hit {i}: pident={pident}, query_start={query_start}, query_end={query_end}, aln_length={aln_length}, hit_IDs={hit_ids}"
                )

                # Check if this is a perfect match (100% identity and covers full query)
                covers_full = (
                    query_start <= 2 and query_end >= sequence_length - 1
                )  # Allow small variations
                is_perfect = (
                    pident >= 99.9 and covers_full
                )  # Allow for floating point precision

                logger.debug(
                    f"Hit {i} analysis: pident={pident}, covers_full={covers_full}, is_perfect={is_perfect}"
                )

                if is_perfect:
                    perfect_blast_match = hit_ids
                    logger.info(
                        f"Found perfect BLAST match: {perfect_blast_match} (pident={pident}, coverage={query_start}-{query_end}/{sequence_length})"
                    )
                    break

            if not perfect_blast_match:
                logger.debug(
                    "No perfect BLAST matches found, will run full classification"
                )

        should_classify = (
            sequence_length >= 5000
            and not analysis.has_error()
            and analysis.blast_result
            and analysis.blast_result.blast_content
            and perfect_blast_match
            is None  # Don't classify if we already have a perfect match
        )

        if perfect_blast_match:
            # Create direct exact match classification
            logger.info(
                f"Skipping classification workflow - using perfect BLAST match: {perfect_blast_match}"
            )
            from src.utils.blast_data import ClassificationData

            classification_data = ClassificationData(
                source="exact",
                closest_match=perfect_blast_match,
                confidence="High",
                match_details=f"Exact sequence match to {perfect_blast_match}",
            )

            # Update pipeline state with the classification
            pipeline_state.update_classification_data(
                tab_sequence_id, classification_data
            )

            # Mark as complete in workflow state
            workflow_state = WorkflowState()
            workflow_state.complete = True
            workflow_state.found_match = True
            workflow_state.match_stage = "exact"
            workflow_state.match_result = perfect_blast_match
            workflow_state.set_classification(classification_data)
            pipeline_state.update_workflow_state(tab_sequence_id, workflow_state)

        elif should_classify:
            logger.info(
                f"Running classification workflow for tab {tab_idx} (length: {sequence_length})"
            )
            try:
                # tab_sequence_state = pipeline_state.add_sequence(tab_sequence_id)

                # Use the FASTA file from the SequenceAnalysis
                tmp_fasta = analysis.blast_result.fasta_file

                if not tmp_fasta:
                    logger.error(
                        f"No FASTA file available from analysis for tab {tab_idx}"
                    )
                    raise Exception("No FASTA file available from sequence analysis")

                blast_data = BlastData(
                    seq_type=analysis.sequence_type.value,
                    fasta_file=tmp_fasta,
                    blast_df=analysis.blast_result.blast_hits,
                    processed_sequences=[0],
                    sequence_results={"0": sequence_result},
                    total_sequences=1,
                )
                pipeline_state.update_blast_data(tab_sequence_id, blast_data)

                workflow_state = WorkflowState(
                    stages={
                        stage["id"]: {"progress": 0, "status": "pending"}
                        for stage in WORKFLOW_STAGES
                    },
                    task_id=tab_sequence_id,
                    workflow_started=False,
                )
                workflow_state.fetch_ship_params = FetchShipParams(
                    curated=False, with_sequence=True, dereplicate=True
                )
                workflow_state.fetch_captain_params = FetchCaptainParams(
                    curated=True, with_sequence=True
                )
                pipeline_state.update_workflow_state(tab_sequence_id, workflow_state)
                classification_data = ClassificationData(
                    seq_type=analysis.sequence_type.value, fasta_file=tmp_fasta
                )
                pipeline_state.update_classification_data(
                    tab_sequence_id, classification_data
                )

                meta_df = fetch_meta_data()
                meta_dict = meta_df.to_dict("records") if meta_df is not None else None

                workflow_result = run_classification_workflow_sync(
                    workflow_state=workflow_state.to_dict(),
                    blast_data=blast_data.to_dict(),
                    classification_data=classification_data.to_dict(),
                    meta_dict=meta_dict,
                )
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

                    enriched_classification = _create_enriched_classification(
                        match_stage, match_accession, meta_df
                    )

                    pipeline_state.update_classification_data(
                        tab_sequence_id, enriched_classification
                    )

                    classification_dict = enriched_classification.to_dict()
                    sequence_result["classification"] = classification_dict

                    logger.info(
                        f"Updated tab {tab_idx} with classification: {classification_dict}"
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
    """
    Render the content for the active tab.
    - If we have blast_data, use it to render the content
    - If we don't have blast_data, show an error
    - If we have blast_data, render the content
    - If we don't have blast_data, show an error
    """
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

        # Get centralized state adapter
        adapter = get_dash_adapter()
        pipeline_state = adapter.pipeline_state

        # Get the main sequence ID and create tab-specific ID
        main_sequence_id = pipeline_state._active_sequence_id
        if not main_sequence_id:
            # Try to find the main sequence ID from existing sequences
            # Look for sequences that don't have "_tab_" in their ID
            main_sequences = [seq_id for seq_id in pipeline_state._sequences.keys() 
                            if "_tab_" not in seq_id]
            if main_sequences:
                main_sequence_id = main_sequences[0]
                logger.info(f"Recovered main sequence ID in render_tab_content: {main_sequence_id}")
                pipeline_state.set_active_sequence(main_sequence_id)
            else:
                logger.warning("No main sequence ID found in render_tab_content")
                main_sequence_id = None
                
        if main_sequence_id:
            tab_sequence_id = f"{main_sequence_id}_tab_{tab_idx}"
        else:
            tab_sequence_id = None

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

        # Get classification data from centralized state if available
        classification_data = None
        if tab_sequence_id:
            # Try to get enriched classification data from centralized state first
            centralized_classification = adapter.get_sequence_classification_for_ui(
                tab_sequence_id
            )
            if centralized_classification:
                classification_data = centralized_classification
                logger.debug(
                    f"Using centralized classification data for tab {tab_idx}: {classification_data}"
                )
            else:
                # Fallback to sequence results classification data
                classification_data = sequence_results.get("classification")
                logger.debug(
                    f"Using sequence results classification data for tab {tab_idx}"
                )
        else:
            # Fallback to sequence results classification data
            classification_data = sequence_results.get("classification")
            logger.debug(f"Using fallback classification data for tab {tab_idx}")

        # Render classification results using centralized data
        classification_output = create_classification_output(
            workflow_state=None,  # Workflow state not needed for tab rendering
            classification_data=classification_data,
        )

        # Render BLAST results
        blast_container = create_blast_container(sequence_results, tab_id=tab_idx)

        logger.debug(
            f"Successfully rendered content for tab {tab_idx} using centralized state"
        )
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
        // Check if we have valid data
        if (!data || !data.blast_text) {
            return window.dash_clientside.no_update;
        }
        
        try {            
            // Find the correct container based on active tab
            let containerId = active_tab_idx !== null ? 
                `blast-container-${active_tab_idx}` : 'blast-container';
            
            // Create a function that will attempt to initialize BlasterJS
            const initializeBlasterJS = (attempts = 0, maxAttempts = 5) => {
                // Use standard ID selector
                let container = document.getElementById(containerId);
                
                // If no container is found and we haven't exceeded max attempts, retry
                if (!container && attempts < maxAttempts) {
                    setTimeout(() => initializeBlasterJS(attempts + 1, maxAttempts), 100);
                    return;
                }
                
                // If still no container after max attempts, try fallback options
                if (!container) {
                    // Try to find the default container first since this is likely a single sequence view
                    container = document.getElementById('blast-container');
                    
                    if (!container) {
                        // If still not found, try to find a container with class blast-container 
                        let containers = document.getElementsByClassName('blast-container');
                        if (containers.length > 0) {
                            container = containers[0];
                        } else {
                            console.error("No blast containers found in the DOM");
                            return;
                        }
                    }
                }
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
                    var instance = new blasterjs({
                        string: data.blast_text,
                        multipleAlignments: alignmentsDivId,
                        alignmentsTable: tablesDivId
                    });
                    
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
                    }, 100);
                } catch (blasterError) {
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
        if (!active_tab) {
            return window.dash_clientside.no_update;
        }
        
        if (!blast_data) {
            return window.dash_clientside.no_update;
        }
        
        if (!blast_data.blast_text) {
            return window.dash_clientside.no_update;
        }

        try {
            // Extract the tab index from the active tab ID
            const tabIdx = parseInt(active_tab.split("-")[1]);
            if (isNaN(tabIdx)) {
                return window.dash_clientside.no_update;
            }

            // Find the correct container based on active tab
            const containerId = `blast-container-${tabIdx}`;
            
            // Create a function that will attempt to initialize BlasterJS
            const initializeBlasterJS = (attempts = 0, maxAttempts = 5) => {
                let container = document.getElementById(containerId);
                
                // If no container is found and we haven't exceeded max attempts, retry
                if (!container && attempts < maxAttempts) {
                    setTimeout(() => initializeBlasterJS(attempts + 1, maxAttempts), 100);
                    return;
                }
                
                // If still no container after max attempts, try fallback options
                if (!container) {
                    // Try fallback options
                    const containers = document.getElementsByClassName('blast-container');
                    if (containers.length > 0) {
                        container = containers[0];
                    } else {
                        const mainContainer = document.getElementById('blast-container');
                        if (mainContainer) {
                            container = mainContainer;
                        } else {
                            return;
                        }
                    }
                }

                // Only initialize if empty or not initialized yet
                if (container.children.length === 0 || !container.dataset.initialized) {
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
                          }, 100);                        
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
    """
    REPLACEMENT for update_single_sequence_classification using centralized state.

    This eliminates the complex data merging logic by using a single source of truth.
    """
    from src.utils.classification_utils import create_classification_output

    # Get the centralized state adapter
    adapter = get_dash_adapter()

    # Convert inputs to objects for compatibility
    blast_results = (
        BlastData.from_dict(blast_results_dict) if blast_results_dict else None
    )
    workflow_state = (
        WorkflowState.from_dict(workflow_state_dict) if workflow_state_dict else None
    )

    if workflow_state:
        logger.info(
            f"Workflow state: complete={workflow_state.complete}, found_match={workflow_state.found_match}, match_stage={workflow_state.match_stage}"
        )

    if not blast_results:
        logger.info("No blast_results, returning None")
        return None

    # Get the active sequence ID (assume sequence 0 for single sequence)
    # In a full implementation, this would come from the adapter's active sequence
    sequence_id = adapter.pipeline_state._active_sequence_id
    if not sequence_id:
        logger.info("No active sequence ID, returning None")
        return None

    # Get classification data from centralized state - SINGLE SOURCE OF TRUTH
    classification_data = adapter.get_sequence_classification_for_ui(sequence_id)

    logger.debug(
        f"update_single_sequence_classification: sequence_id={sequence_id}, classification_data={classification_data is not None}"
    )

    if classification_data:
        logger.info(
            f"Found classification data in centralized state: {classification_data}"
        )
    else:
        logger.info("No classification data found in centralized state")

    # Create the classification output using the unified data
    result = create_classification_output(
        workflow_state=workflow_state, classification_data=classification_data
    )

    logger.info(f"create_classification_output returned: {type(result)}")
    return result


@callback(
    [
        Output("workflow-state-store", "data", allow_duplicate=True),
        Output("results-loading-overlay", "visible", allow_duplicate=True),
        Output("blast-data-store", "data", allow_duplicate=True),
    ],
    [
        Input("workflow-state-store", "data"),
        Input("classification-data-store", "data"),
    ],
    [
        State("blast-data-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_classification_workflow_state(
    workflow_state_dict, classification_data_dict, blast_results_dict
):
    """
    Update the classification workflow state using the centralized pipeline state.
    """

    # Early exit conditions
    if workflow_state_dict is None or classification_data_dict is None:
        raise PreventUpdate

    from src.utils.blast_data import ClassificationData

    adapter = get_dash_adapter()
    pipeline_state = adapter.pipeline_state

    workflow_state = WorkflowState.from_dict(workflow_state_dict)
    classification_data = ClassificationData.from_dict(classification_data_dict)
    # blast_results = (
    #     BlastData.from_dict(blast_results_dict) if blast_results_dict else None
    # )

    # Get sequence ID from workflow state
    sequence_id = workflow_state.task_id
    if not sequence_id:
        logger.warning("No sequence ID found in workflow state")
        raise PreventUpdate

    # Check if workflow should run (one-time operation)
    if (
        workflow_state.complete
        or workflow_state.error is not None
        or workflow_state.workflow_started
    ):
        if workflow_state.complete or workflow_state.error is not None:
            # Update centralized state and return final data
            pipeline_state.update_workflow_state(sequence_id, workflow_state)
            store_data = adapter.sync_all_stores(sequence_id)
            return store_data["workflow_state"], False, store_data["blast_data"]
        raise PreventUpdate

    try:
        logger.debug("Running classification workflow via centralized state")

        # Get current state from centralized pipeline
        resolved_sequence_id = pipeline_state.resolve_sequence_id(sequence_id)
        if resolved_sequence_id:
            sequence_id = resolved_sequence_id
            logger.info(f"Resolved sequence ID: {sequence_id}")

        sequence_state = pipeline_state.get_sequence(sequence_id)
        # active_sequence_id = pipeline_state._active_sequence_id

        if not sequence_state:
            logger.error(f"No sequence state found for {sequence_id}")
            # Create a new sequence state to prevent further errors
            pipeline_state.add_sequence_without_activation(sequence_id)
            sequence_state = pipeline_state.get_sequence(sequence_id)
            logger.info(f"Created new sequence state for {sequence_id}")

        # Prepare data for workflow
        blast_data = sequence_state.blast_data
        if not blast_data:
            logger.error(f"No BLAST data found for {sequence_id}")
            raise PreventUpdate

        # Get metadata
        meta_df = fetch_meta_data()
        meta_dict = meta_df.to_dict("records") if meta_df is not None else None

        # Run workflow directly (no Celery)
        result = run_classification_workflow_sync(
            workflow_state=workflow_state.to_dict(),
            blast_data=blast_data.to_dict(),
            classification_data=classification_data.to_dict(),
            meta_dict=meta_dict,
        )

        logger.debug(f"Workflow result type: {type(result)}")

        # Convert result back to workflow state
        if isinstance(result, dict):
            updated_workflow_state = WorkflowState.from_dict(result)
        else:
            updated_workflow_state = result

        # Mark as started and complete
        updated_workflow_state.workflow_started = True
        updated_workflow_state.status = "complete"
        updated_workflow_state.complete = True

        # Update centralized state with results
        pipeline_state.update_workflow_state(sequence_id, updated_workflow_state)

        # Handle workflow completion
        if updated_workflow_state.complete:
            if (
                updated_workflow_state.found_match
                and updated_workflow_state.match_result
            ):
                logger.info(
                    f"Workflow found match: {updated_workflow_state.match_stage} -> {updated_workflow_state.match_result}"
                )

                # Get existing classification data from workflow state
                existing_classification = None
                if (
                    hasattr(updated_workflow_state, "classification_data")
                    and updated_workflow_state.classification_data
                ):
                    from src.utils.blast_data import ClassificationData

                    existing_classification = ClassificationData.from_dict(
                        updated_workflow_state.classification_data
                    )

                # Create enriched classification data with metadata lookup
                enriched_classification = _create_enriched_classification(
                    updated_workflow_state.match_stage,
                    updated_workflow_state.match_result,
                    meta_df,
                    existing_classification,
                )

                # Update centralized state with enriched classification - SINGLE UPDATE
                pipeline_state.update_classification_data(
                    sequence_id, enriched_classification
                )

                logger.info(
                    f"Updated centralized state with enriched classification: {enriched_classification}"
                )
            else:
                logger.info(
                    "Workflow completed but no matches found - clearing classification data"
                )
                # Clear any empty classification data when no matches are found
                sequence_state = pipeline_state.get_sequence(sequence_id)
                if sequence_state:
                    sequence_state.classification_data = None
                    logger.debug(
                        "Cleared empty classification data from centralized state"
                    )

        # Return synchronized data from centralized state
        store_data = adapter.sync_all_stores(sequence_id)
        return store_data["workflow_state"], False, store_data["blast_data"]

    except Exception as e:
        logger.error(f"Error in unified workflow state update: {e}")

        # Update error state in centralized system
        workflow_state.error = str(e)
        workflow_state.status = "failed"
        workflow_state.complete = True
        pipeline_state.update_workflow_state(sequence_id, workflow_state)
        pipeline_state.set_sequence_error(sequence_id, str(e))

        # Return error state
        store_data = adapter.sync_all_stores(sequence_id)
        return store_data["workflow_state"], False, store_data["blast_data"]


def _create_enriched_classification(
    match_stage, match_accession, meta_df, existing_classification=None
):
    """Create ClassificationData with metadata lookup"""
    from src.utils.blast_data import ClassificationData

    # Start with existing classification data if available, otherwise create new
    if existing_classification:
        classification_data = ClassificationData(
            source=existing_classification.source or match_stage,
            closest_match=existing_classification.closest_match or match_accession,
            confidence=existing_classification.confidence
            or ("High" if match_stage in ["exact", "contained"] else "Medium"),
            match_details=existing_classification.match_details,  # Preserve existing match_details
            family=existing_classification.family,
            navis=existing_classification.navis,
            haplotype=existing_classification.haplotype,
        )
    else:
        # Create base classification
        classification_data = ClassificationData(
            source=match_stage,
            closest_match=match_accession,
            confidence="High" if match_stage in ["exact", "contained"] else "Medium",
        )

    # Look up metadata
    try:
        if meta_df is not None and not meta_df.empty:
            meta_match = meta_df[meta_df["accession_display"] == match_accession]
            if not meta_match.empty:
                logger.debug(f"Found metadata for {match_accession}")

                # Add family, navis, haplotype info
                for col, attr in [
                    ("familyName", "family"),
                    ("navis_name", "navis"),
                    ("haplotype_name", "haplotype"),
                ]:
                    if col in meta_match.columns:
                        value = meta_match[col].iloc[0]
                        if value and value != "None":
                            setattr(classification_data, attr, value)

                # Add match details (only if not already set by workflow)
                if (
                    not hasattr(classification_data, "match_details")
                    or not classification_data.match_details
                ):
                    if match_stage == "exact":
                        classification_data.match_details = (
                            f"Exact sequence match to {match_accession}"
                        )
                    elif match_stage == "contained":
                        # Check if it's a perfect match (High confidence indicates perfect match)
                        if (
                            hasattr(classification_data, "confidence")
                            and classification_data.confidence == "High"
                        ):
                            classification_data.match_details = (
                                f"Perfect sequence match to {match_accession}"
                            )
                        else:
                            classification_data.match_details = (
                                f"Query sequence contained within {match_accession}"
                            )
                    elif match_stage == "similar":
                        classification_data.match_details = (
                            f"High similarity to {match_accession}"
                        )
                    else:
                        classification_data.match_details = f"{match_stage.replace('_', ' ').title()} match to {match_accession}"
            else:
                logger.warning(f"No metadata found for {match_accession}")
    except Exception as e:
        logger.error(f"Error looking up metadata for {match_accession}: {e}")

    return classification_data


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

import dash
from dash import dcc, html, callback, clientside_callback
from dash_iconify import DashIconify
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import tempfile
import base64
import pandas as pd
import logging
import time

from src.config.cache import cache
from src.utils.seq_utils import (
    check_input,
    write_temp_fasta,
)
from src.utils.blast_utils import (
    run_blast,
    process_captain_results,
    create_no_matches_alert,
    parse_blast_xml,
    blast_download_button,
)

from src.components.callbacks import (
    curated_switch,
    create_file_upload,
    create_accession_modal,
)
from src.database.sql_manager import fetch_meta_data
from src.utils.telemetry import blast_limit_decorator
from src.components.error_boundary import handle_callback_error, create_error_alert
from src.utils.classification_utils import WORKFLOW_STAGES
from src.tasks import (
    run_blast_search_task,
    run_hmmer_search_task,
    run_classification_workflow_task,
)

dash.register_page(__name__)

logger = logging.getLogger(__name__)


def blast_family_button(family):
    return dbc.Button(
        family,
        color="primary",
        href=f"/wiki?page={family}",
        external_link=False,
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
        dcc.Store(id="captain-results-store", data=None),
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
        # Timeout stores
        dcc.Store(id="blast-timeout-store", data=False),
        dcc.Interval(
            id="blast-timeout-interval", interval=30000, n_intervals=0
        ),  # 30 seconds
        # Interval for polling workflow state
        dcc.Interval(
            id="classification-interval",
            interval=1000,  # 1 second
            disabled=True,
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
                                            dmc.Title("Input Sequence", order=3),
                                            dmc.Textarea(
                                                id="query-text",
                                                placeholder="Paste FASTA sequence here...",
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
                                                            placeholder_text="Drag and drop or click to select a FASTA file",
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
                        html.Div(
                            id="right-column-content",
                            children=[
                                # This div will be replaced with tabs when more than one sequence is in query
                                dmc.Stack(
                                    children=[
                                        # Progress section - initially hidden
                                        html.Div(
                                            id="classification-output", className="mt-4"
                                        ),
                                        dmc.Stack(
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
                                        ),
                                        # BLAST results section
                                        dmc.Stack(
                                            [
                                                html.Div(
                                                    id="blast-container",
                                                    style={
                                                        "width": "100%",
                                                        "display": "flex",
                                                        "flexDirection": "column",
                                                    },
                                                ),
                                                html.Div(id="blast-download"),
                                            ]
                                        ),
                                    ],
                                ),
                            ],
                        ),
                        # Debug output (will be hidden in production)
                        html.Div(
                            id="right-column-content-debug",
                            style={"margin-top": "10px", "color": "red"},
                        ),
                    ],
                    style={"justify-content": "flex-start"},
                ),
            ],
            gutter="xl",
        ),
    ],
)

# Define the clientside JavaScript function directly in the callback
clientside_callback(
    """
    function(data) {
        if (data && data.blast_text) {
            try {
                // Get the blast container
                const container = document.getElementById('blast-container');
                
                // Clear existing content first
                container.innerHTML = '';
                
                // Create the title element
                const titleElement = document.createElement('h2');
                titleElement.innerHTML = 'BLAST Results';
                titleElement.style.marginTop = '15px';
                titleElement.style.marginBottom = '20px';
                titleElement.style.fontWeight = '500';
                titleElement.style.textAlign = 'left';
                
                // Create divs for BlasterJS to use - in vertical stack
                const alignmentsDiv = document.createElement('div');
                alignmentsDiv.id = 'blast-multiple-alignments';
                alignmentsDiv.style.marginBottom = '30px';
                
                const tableDiv = document.createElement('div');
                tableDiv.id = 'blast-alignments-table';
                
                // Add elements to the container in order (vertically stacked)
                container.appendChild(titleElement);
                container.appendChild(alignmentsDiv);
                container.appendChild(tableDiv);
                
                // Initialize BlasterJS
                var blasterjs = require("biojs-vis-blasterjs");
                var instance = new blasterjs({
                    string: data.blast_text,
                    multipleAlignments: "blast-multiple-alignments",
                    alignmentsTable: "blast-alignments-table"
                });
                
                // Add CSS to control width and alignment
                var style = document.createElement('style');
                style.textContent = `
                    #blast-container {
                        width: 100%;
                        display: flex;
                        flex-direction: column;
                        align-items: center;
                    }
                    
                    /* For the title - ensure it's left-aligned */
                    #blast-container > h2 {
                        align-self: flex-start;
                        width: 100%;
                        text-align: left;
                    }
                    
                    /* For alignment visualization section */
                    #blast-multiple-alignments {
                        width: auto !important;
                        display: inline-block !important;
                        margin: 0 auto !important;
                        text-align: center;
                        margin-bottom: 30px !important;
                    }
                    
                    /* Center the title and buttons */
                    #blast-multiple-alignments > div:first-child {
                        text-align: center !important;
                    }
                    
                    /* Make the color key centered */
                    #blast-multiple-alignments div[style*="color: #EEE"] {
                        margin: 0 auto !important;
                        text-align: center !important;
                        display: table !important;
                    }
                    
                    /* Make QUERY and numerical scale centered */
                    #blast-multiple-alignments div[style*="color: #5C6D7E"] {
                        margin: 0 auto !important;
                        display: table !important;
                    }
                    
                    /* Center the buttons area */
                    #blast-multiple-alignments-buttons {
                        text-align: center !important;
                    }
                    
                    /* But left-align the actual alignments */
                    #alignments-container > div:last-child {
                        text-align: left !important;
                    }
                    
                    /* For table section - center the container */
                    #blast-alignments-table {
                        display: inline-block !important;
                        width: auto !important;
                        text-align: center;
                        margin: 0 auto !important;
                    }
                    
                    /* But left-align the actual table content */
                    #blast-alignments-table-img {
                        text-align: left !important;
                        margin: 0 auto !important;
                    }
                    
                    #blast-alignments-table table {
                        width: 100%;
                        table-layout: fixed;
                        margin: 0;
                        text-align: left;
                    }
                    
                    #blast-alignments-table td {
                        max-width: 200px;
                        overflow: hidden;
                        text-overflow: ellipsis;
                        white-space: nowrap;
                        text-align: left;
                    }
                    
                    /* Make the table buttons centered */
                    #blast-alignments-table > div:last-child {
                        text-align: center !important;
                        margin: 15px auto !important;
                    }
                    
                    /* Handle overflow for the table */
                    #blast-alignments-table-img {
                        max-width: 100% !important;
                        overflow-x: auto;
                    }
                `;
                document.head.appendChild(style);
                
                // Add a setTimeout to fix any alignment issues after rendering
                setTimeout(() => {
                    // Handle the alignments container
                    const alignmentsContainer = document.getElementById('alignments-container');
                    if (alignmentsContainer) {
                        // Ensure buttons are centered
                        const buttonsDiv = document.getElementById('blast-multiple-alignments-buttons');
                        if (buttonsDiv) {
                            buttonsDiv.style.textAlign = 'center';
                        }
                    }
                    
                    // Handle the table container
                    const tableImgDiv = document.getElementById('blast-alignments-table-img');
                    if (tableImgDiv) {
                        // Make sure table is left-aligned inside its container
                        tableImgDiv.style.textAlign = 'left';
                        
                        // Handle the download buttons div
                        const buttonsDiv = tableImgDiv.nextElementSibling;
                        if (buttonsDiv) {
                            buttonsDiv.style.textAlign = 'center';
                            buttonsDiv.style.width = '100%';
                            buttonsDiv.style.margin = '15px auto';
                        }
                    }
                }, 500);
                
                return window.dash_clientside.no_update;
            } catch (error) {
                console.error('Error initializing BlasterJS:', error);
                return error.toString();
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("blast-container", "children"),
    Input("processed-blast-store", "data"),
    prevent_initial_call=True,
)


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
@handle_callback_error
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
        logger.info(f"Input type: {input_type}, Number of sequences: {n_seqs}")

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
    [
        Output("submit-button", "disabled", allow_duplicate=True),
        Output("submit-button", "children", allow_duplicate=True),
        Output("upload-error-message", "children", allow_duplicate=True),
        Output("upload-error-store", "data", allow_duplicate=True),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("blast-timeout-interval", "n_intervals"),
    ],
    [State("blast-timeout-store", "data"), State("blast-sequences-store", "data")],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_blast_timeout(n_clicks, n_intervals, timeout_triggered, seq_list):
    triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

    # Default return values
    button_disabled = False
    button_text = "Submit BLAST"  # Set default button text
    error_message = ""
    error_store = None

    # Handle timeout case
    if triggered_id == "blast-timeout-interval":
        if n_clicks and timeout_triggered:
            button_disabled = True
            button_text = "Server timeout"
            error_message = dmc.Alert(
                title="Request Timeout",
                children="The server is taking longer than expected to respond. Please try again later.",
                color="yellow",
                variant="filled",
            )
            error_store = "The server is taking longer than expected to respond. Please try again later."

    return [button_disabled, button_text, error_message, error_store]


@callback(
    [
        Output("blast-sequences-store", "data", allow_duplicate=True),
        Output("upload-error-store", "data", allow_duplicate=True),
        Output("upload-error-message", "children", allow_duplicate=True),
    ],
    [
        Input("submit-button", "n_clicks"),
    ],
    [
        State("query-text", "value"),
        State("blast-fasta-upload", "contents"),
        State("blast-sequences-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def preprocess(n_clicks, query_text_input, query_file_contents, seq_list):
    """
    Process input when the submit button is pressed.
    For file uploads: Use the already parsed sequences from blast-sequences-store
    For text input: Parse the input text now
    """
    if not n_clicks:
        raise PreventUpdate

    error_alert = None

    try:
        # If we have parsed sequences from a file upload, use those
        if seq_list is not None:
            logger.info(
                f"Using pre-parsed sequences from file upload: {len(seq_list)} sequences"
            )
            return seq_list, None, None

        # Otherwise, parse the text input
        if query_text_input:
            logger.info("Processing text input")
            input_type, seq_list, n_seqs, error = check_input(
                query_text_input=query_text_input, query_file_contents=None
            )

            if error:
                logger.error(f"Error parsing text input: {error}")
                error_alert = (
                    error
                    if not isinstance(error, str)
                    else dmc.Alert(
                        title="Error Processing Sequence",
                        children=error,
                        color="red",
                        variant="filled",
                    )
                )
                return (
                    None,
                    str(error) if isinstance(error, str) else "Parsing error",
                    error_alert,
                )

            logger.info(f"Text input processed: {input_type}, {n_seqs} sequences")
            return seq_list, None, None

        # If we have neither text input nor uploaded sequences, show an error
        error_alert = dmc.Alert(
            title="No Input Provided",
            children="Please enter a sequence or upload a FASTA file.",
            color="yellow",
            variant="filled",
        )
        return None, "No input provided", error_alert

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        error_alert = dmc.Alert(
            title="Error Processing Input",
            children=f"An unexpected error occurred: {str(e)}",
            color="red",
            variant="filled",
        )
        return None, str(e), error_alert


@callback(
    Output("submission-id-store", "data"),
    Input("submit-button", "n_clicks"),
    prevent_initial_call=True,
)
def update_submission_id(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return n_clicks  # Use n_clicks as a unique submission ID


# Metadata Processing Callback
@callback(
    Output("processed-metadata-store", "data"),
    Input("curated-input", "value"),
)
@handle_callback_error
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


# BLAST Results Processing Callback
@callback(
    Output("processed-blast-store", "data"),
    Input("blast-results-store", "data"),
)
@handle_callback_error
def process_blast_results(blast_results_file):
    if not blast_results_file:
        return None

    try:
        # Read the BLAST output as text
        with open(blast_results_file, "r") as f:
            blast_results = f.read()

        # Add size limit check
        results_size = len(blast_results)
        if results_size > 5 * 1024 * 1024:  # 5MB limit
            logger.warning(f"BLAST results too large: {results_size} bytes")
            return None

        # Format data for BlasterJS
        data = {
            "blast_text": blast_results  # Pass the raw BLAST text directly
        }

        return data
    except Exception as e:
        logger.error(f"Error processing BLAST results: {str(e)}")
        return None


@blast_limit_decorator
@callback(
    [
        Output("blast-results-store", "data"),
        Output("captain-results-store", "data"),
        Output("classification-interval", "disabled", allow_duplicate=True),
        Output("workflow-state-store", "data"),
        Output("upload-error-message", "children", allow_duplicate=True),
    ],
    [
        Input("submission-id-store", "data"),
        Input("blast-sequences-store", "data"),
        Input("blast-active-tab", "data"),
    ],
    [
        State("evalue-threshold", "value"),
        State("blast-sequences-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def blast_callback(
    submission_id,
    seq_list,
    active_tab,
    evalue_threshold,
    search_type="hmmsearch",
):
    if not all([seq_list, submission_id]):
        return None, None, True, None, None

    try:
        if isinstance(seq_list, list) and len(seq_list) > 0:
            if len(seq_list) > active_tab:
                if active_tab is None:
                    active_tab = 0
                seq_data = seq_list[active_tab]
                query_header = seq_data.get("header", "query")
                query_seq = seq_data.get("sequence", "")
                query_type = seq_data.get("type", "nucl")
                seq_list[active_tab]["processed"] = True
        else:
            # Handle case where seq_list isn't a list (shouldn't happen but just in case)
            error_alert = dmc.Alert(
                title="Invalid Sequence Data",
                children="There was a problem with the sequence data. Please try again.",
                color="red",
                variant="filled",
            )
            return None, None, True, None, error_alert

        # Ensure query_type is a string
        if not isinstance(query_type, str):
            query_type = str(query_type)

        # Check sequence length - skip detailed classification for short sequences
        skip_classification = False
        if query_seq and len(query_seq) < 5000:  # 5kb threshold
            logger.info(
                f"Sequence length {len(query_seq)} is below 5kb threshold, skipping detailed classification"
            )
            skip_classification = True

        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)

        # Start task for BLAST search - Now using celery task
        blast_task = run_blast_search_task.delay(
            query_header=query_header,
            query_seq=query_seq,
            query_type=query_type,
            eval_threshold=evalue_threshold,
        )

        blast_results_file = blast_task.get(timeout=300)  # Wait for BLAST to complete

        # Run Diamond and HMMER using celery tasks
        if query_type == "nucl":
            # For nucleotide sequences, use the HMMER task
            hmmer_task = run_hmmer_search_task.delay(
                query_header=query_header,
                query_seq=query_seq,
                query_type=query_type,
                eval_threshold=evalue_threshold,
            )

            hmmer_results = hmmer_task.get(timeout=300)

            if hmmer_results:
                captain_results_dict = hmmer_results.get("results")
        else:
            # For protein sequences, directly run HMMER task
            hmmer_task = run_hmmer_search_task.delay(
                query_header=query_header,
                query_seq=query_seq,
                query_type=query_type,
                eval_threshold=evalue_threshold,
            )

            hmmer_results = hmmer_task.get(timeout=300)

            if hmmer_results:
                captain_results_dict = hmmer_results.get("results")

        if blast_results_file is None:
            return (
                None,
                None,
                True,
                None,
                None,  # No error
            )

        # Rest of your existing code...
        # Create a new workflow tracker - even if we're skipping classification
        # so we have the proper state in the UI
        globals()["workflow_tracker"] = WorkflowStatus()
        globals()["workflow_tracker"].start_time = time.time()

        if skip_classification:
            # Skip classification for short sequences
            # Mark workflow as complete immediately using a special flag
            globals()["workflow_tracker"].complete = True
            globals()["workflow_tracker"].match_stage = "blast_only"
            globals()[
                "workflow_tracker"
            ].match_result = (
                f"Sequence length ({len(query_seq)}bp) is below 5kb threshold"
            )

            # Return without starting the classification thread
            return (
                blast_results_file,
                captain_results_dict,
                None,
                True,  # Disable interval since we're not doing classification
                None,  # No error
            )

        # For sequences >= 5kb, continue with the regular classification workflow
        upload_data = {
            "seq_type": query_type,
            "fasta": tmp_query_fasta,
            "fetch_ship_params": {
                "curated": False,
                "with_sequence": True,
                "dereplicate": True,
            },
            "fetch_captain_params": {"curated": True, "with_sequence": True},
        }

        # This launches the entire sequential workflow as a single background task
        classification_task = run_classification_workflow_task.delay(upload_data)

        # Return initial state with task ID for polling
        return (
            blast_results_file,
            captain_results_dict,
            False,  # Enable the interval for polling
            {
                "task_id": classification_task.id,
                "status": "started",
                "complete": False,
                "current_stage": None,
                "error": None,
                "start_time": time.time(),
            },
            None,  # No error
        )

    except Exception as e:
        logger.error(f"Error in blast_callback: {str(e)}")
        error_alert = create_error_alert(str(e))
        return None, None, True, None, error_alert


@callback(
    [
        Output("classification-output", "children"),
        Output("blast-download", "children"),
    ],
    [
        Input("blast-results-store", "data"),
        Input("classification-workflow-state", "data"),
        Input("classification-result-store", "data"),
    ],
    [
        State("blast-sequences-store", "data"),
        State("captain-results-store", "data"),
        State("submit-button", "n_clicks"),
        State("evalue-threshold", "value"),
        State("blast-active-tab", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def update_ui_elements(
    blast_results_file,
    classification_state,
    classification_result,
    captain_results_dict,
    seq_list,
    n_clicks,
    evalue,
    active_tab,
):
    if not n_clicks or blast_results_file is None:
        return None, None

    # Check if we're working with more than one sequence
    is_multifasta = seq_list and len(seq_list) > 1

    try:
        if not blast_results_file:
            return html.Div(
                [
                    html.H2("Classification Results", style={"marginBottom": "20px"}),
                    create_no_matches_alert(),
                ]
            ), None

        # Initialize classification data
        classification_data = None

        # Check for classification result data first
        if classification_result is not None:
            result_type = classification_result.get("type")

            if result_type == "match" and classification_result.get("match_stage") in [
                "exact",
                "contained",
                "similar",
            ]:
                # Handle exact/contained/similar matches
                match_stage = classification_result.get("match_stage")
                match_result = classification_result.get("match_result")

                if match_result:
                    # Get metadata for the matched ship
                    meta_df = fetch_meta_data(accession_tag=match_result)

                    if not meta_df.empty:
                        # We found metadata for this match
                        family = (
                            meta_df["familyName"].iloc[0]
                            if "familyName" in meta_df.columns
                            else None
                        )
                        navis = (
                            meta_df["starship_navis"].iloc[0]
                            if "starship_navis" in meta_df.columns
                            else None
                        )
                        haplotype = (
                            meta_df["starship_haplotype"].iloc[0]
                            if "starship_haplotype" in meta_df.columns
                            else None
                        )

                        # Set confidence based on match type
                        confidence = (
                            "High"
                            if match_stage == "exact"
                            else ("Medium" if match_stage == "contained" else "Low")
                        )

                        classification_data = {
                            "title": "Classification based on",
                            "source": f"{match_stage}_match",
                            "family": family,
                            "navis": navis,
                            "haplotype": haplotype,
                            "match_type": f"Matched to {match_result}",
                            "confidence": confidence,
                        }

            # Alternatively if classification state is available and shows a completed process
            elif classification_state and classification_state.get("complete", False):
                # Process based on classification state (existing logic)

                # Check for workflow-based classification (highest priority)
                if classification_state.get("match_stage") in [
                    "exact",
                    "contained",
                    "similar",
                ]:
                    match_stage = classification_state.get("match_stage")
                    match_result = classification_state.get("match_result")

                    if match_result:
                        # Get metadata for the matched ship
                        meta_df = fetch_meta_data(accession_tag=match_result)

                        if not meta_df.empty:
                            # We found metadata for this match
                            family = (
                                meta_df["familyName"].iloc[0]
                                if "familyName" in meta_df.columns
                                else None
                            )
                            navis = (
                                meta_df["starship_navis"].iloc[0]
                                if "starship_navis" in meta_df.columns
                                else None
                            )
                            haplotype = (
                                meta_df["starship_haplotype"].iloc[0]
                                if "starship_haplotype" in meta_df.columns
                                else None
                            )

                            # Set confidence based on match type
                            confidence = (
                                "High"
                                if match_stage == "exact"
                                else ("Medium" if match_stage == "contained" else "Low")
                            )

                            classification_data = {
                                "title": "Classification based on",
                                "source": f"{match_stage}_match",
                                "family": family,
                                "navis": navis,
                                "haplotype": haplotype,
                                "match_type": f"Matched to {match_result}",
                                "confidence": confidence,
                            }

                # Check for family/navis/haplotype classification from workflow
                elif classification_state.get("found_match", False):
                    match_stage = classification_state.get("match_stage")
                    match_result = classification_state.get("match_result")

                    if match_stage == "family" and match_result:
                        # Family classification found
                        if isinstance(match_result, dict) and "family" in match_result:
                            family_name = match_result["family"]
                        else:
                            family_name = match_result

                        classification_data = {
                            "title": "Family Classification",
                            "source": "classification",
                            "family": family_name,
                            "match_type": "Direct classification",
                            "confidence": "Medium",
                        }
                    elif match_stage == "navis" and match_result:
                        # Navis classification found
                        classification_data = {
                            "title": "Navis Classification",
                            "source": "classification",
                            "navis": match_result,
                            "match_type": "Direct classification",
                            "confidence": "Medium",
                        }
                    elif match_stage == "haplotype" and match_result:
                        # Haplotype classification found
                        classification_data = {
                            "title": "Haplotype Classification",
                            "source": "classification",
                            "haplotype": match_result,
                            "match_type": "Direct classification",
                            "confidence": "Medium",
                        }

        # If no classification from workflow, fall back to BLAST results
        if classification_data is None:
            # Process the BLAST results
            blast_tsv = parse_blast_xml(blast_results_file)
            blast_df = pd.read_csv(blast_tsv, sep="\t")

            # Check if dataframe is empty
            if len(blast_df) == 0:
                classification_output = html.Div([create_no_matches_alert()])
                return classification_output, None

            # Sort by evalue (ascending) and pident (descending) to get best hits
            blast_df = blast_df.sort_values(
                ["evalue", "pident"], ascending=[True, False]
            )
            top_hit = blast_df.iloc[0]
            logger.info(f"Top hit: {top_hit}")

            top_evalue = float(top_hit["evalue"])
            top_aln_length = int(top_hit["aln_length"])
            top_pident = float(top_hit["pident"])

            if top_pident >= 90:
                # look up family name from accession tag
                hit_IDs = top_hit["hit_IDs"]
                meta_df = fetch_meta_data(accession_tag=hit_IDs)

                if not meta_df.empty:
                    top_family = meta_df["familyName"].iloc[0]
                    navis = (
                        meta_df["starship_navis"].iloc[0]
                        if "starship_navis" in meta_df.columns
                        else None
                    )
                    haplotype = (
                        meta_df["starship_haplotype"].iloc[0]
                        if "starship_haplotype" in meta_df.columns
                        else None
                    )

                    classification_data = {
                        "title": "Classification from BLAST",
                        "source": "blast_hit",
                        "family": top_family,
                        "navis": navis,
                        "haplotype": haplotype,
                        "match_type": f"BLAST hit length {top_aln_length}bp with {top_pident:.1f}% identity to {hit_IDs}. Evalue: {top_evalue}",
                        "confidence": "Medium" if top_pident >= 95 else "Low",
                    }
                else:
                    # No metadata for this hit, fall back to captain results
                    captain_family = process_captain_results(
                        captain_results_dict, evalue, return_raw=True
                    )
                    if captain_family:
                        classification_data = {
                            "title": "Captain Gene Classification",
                            "source": "classification",
                            "family": captain_family,
                            "match_type": "Captain gene classification",
                            "confidence": "Low",
                        }
            else:
                # BLAST hit below threshold, fall back to captain results
                captain_family = process_captain_results(
                    captain_results_dict, evalue, return_raw=True
                )
                if captain_family:
                    classification_data = {
                        "title": "Captain Gene Classification",
                        "source": "classification",
                        "family": captain_family,
                        "match_type": "Captain gene classification",
                        "confidence": "Low",
                    }

        # Create the classification card with title
        if classification_data:
            classification_output = html.Div(
                [
                    dmc.Title(
                        "Classification Results",
                        order=2,
                        style={"marginBottom": "20px"},
                    ),
                    create_classification_card(classification_data),
                ]
            )
        else:
            # No classification available
            classification_output = html.Div(
                [
                    dmc.Title(
                        "Classification Results",
                        order=2,
                        style={"marginBottom": "20px"},
                    ),
                    dmc.Alert(
                        title="No Classification Available",
                        children="Could not classify this sequence with any available method.",
                        color="yellow",
                        variant="light",
                    ),
                ]
            )

        # Create download button
        download_button = blast_download_button()

        # If we're working with more than one sequence, we may need to handle tabs
        if is_multifasta:
            # Check if we've already created tabs
            ctx = dash.callback_context
            triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

            # Return the outputs to be rendered in the appropriate tab
            if triggered_id == "blast-sequences-store":
                # Initial load - create the tab structure first
                return None, download_button

        # Return regular output for single sequence or active tab
        return classification_output, download_button

    except Exception as e:
        logger.error(f"Error in update_ui_elements: {e}")
        logger.exception("Full traceback:")
        error_div = html.Div(
            [
                dmc.Title(
                    "Classification Results", order=2, style={"marginBottom": "20px"}
                ),
                create_error_alert(str(e)),
            ]
        )
        return error_div, None


@callback(
    Output("blast-timeout-store", "data"),
    [Input("blast-timeout-interval", "n_intervals")],
    [
        State("submit-button", "n_clicks"),
        State("blast-container", "children"),
        State("upload-error-message", "children"),
        State("classification-output", "children"),
    ],
)
def check_timeout(
    n_intervals, n_clicks, blast_container, error_content, classification_output
):
    try:
        if n_clicks and n_intervals > 0:
            has_output = any(
                [
                    blast_container is not None and blast_container != [],
                    error_content is not None,
                    classification_output is not None,
                ]
            )
            return not has_output
        return False
    except Exception as e:
        logger.error(f"Error in check_timeout: {str(e)}")
        return False


clientside_callback(
    """
    function(n_intervals) {
        window.addEventListener('blastAccessionClick', function(event) {
            if (event.detail && event.detail.accession) {
                // Find the store and update it
                var store = document.getElementById('blast-modal-accession');
                if (store) {
                    store.setAttribute('data-dash-store', JSON.stringify(event.detail.accession));
                }
            }
        });
        return window.dash_clientside.no_update;
    }
    """,
    Output("modal-event-handler", "children"),
    Input("interval-component", "n_intervals"),
)


@callback(
    Output("blast-modal-content", "children"),
    Input("blast-modal-accession", "data"),
    prevent_initial_call=True,
)
def update_modal_content(accession):
    if not accession:
        raise PreventUpdate

    try:
        modal_content, _ = create_accession_modal(accession)
        return modal_content
    except Exception as e:
        logger.error(f"Error in update_modal_content: {str(e)}")
        return html.Div(
            f"Error loading details: {str(e)}",
            style={"color": "red", "padding": "20px"},
        )


# Define a class for tracking workflow status
class WorkflowStatus:
    def __init__(self):
        self.start_time = None
        self.current_stage = None
        self.current_stage_idx = None
        self.stage_progress = None
        self.found_match = False
        self.match_stage = None
        self.match_result = None
        self.complete = False
        self.error = None
        self.stage_values = [0] * len(WORKFLOW_STAGES)
        self.stage_animated = [False] * len(WORKFLOW_STAGES)
        self.stage_striped = [True] * len(WORKFLOW_STAGES)

    def to_dict(self):
        """Convert status to dictionary for storing in dcc.Store"""
        return {
            "start_time": self.start_time,
            "current_stage": self.current_stage,
            "current_stage_idx": self.current_stage_idx,
            "stage_progress": self.stage_progress,
            "found_match": self.found_match,
            "match_stage": self.match_stage,
            "match_result": self.match_result,
            "complete": self.complete,
            "error": self.error,
            "stage_values": self.stage_values,
            "stage_animated": self.stage_animated,
            "stage_striped": self.stage_striped,
        }


# Global variable to store workflow status
# This is the key to sharing progress between the long-running process and the UI
workflow_tracker = None


@callback(
    [
        Output("classification-stage-display", "children"),
        Output("classification-progress", "value"),
        Output("classification-progress", "animated"),
        Output("classification-progress", "striped"),
        Output("classification-interval", "disabled", allow_duplicate=True),
    ],
    Input("workflow-state-store", "data"),
    prevent_initial_call=True,
)
@handle_callback_error
def update_ui_from_workflow_state(workflow_state):
    """Update UI elements based on the workflow state."""
    if not workflow_state:
        raise PreventUpdate

    # Set default values
    stage_message = "Starting classification..."
    overall_progress = 0
    animated = True
    striped = True
    disable_interval = False

    # Update based on current state
    if workflow_state.get("complete", False):
        disable_interval = True
        animated = False
        striped = False
        overall_progress = 100

        if workflow_state.get("error"):
            stage_message = f"Error: {workflow_state['error']}"
        elif workflow_state.get("found_match", False):
            match_stage = workflow_state.get("match_stage", "")
            match_result = workflow_state.get("match_result", "")
            stage_info = next(
                (s for s in WORKFLOW_STAGES if s["id"] == match_stage), {}
            )
            stage_message = f"{stage_info.get('label', 'Match')} found: {match_result}"
        else:
            stage_message = "Classification complete - No matches found"
    elif workflow_state.get("current_stage"):
        # Still processing
        current_stage = workflow_state["current_stage"]
        stage_info = next((s for s in WORKFLOW_STAGES if s["id"] == current_stage), {})
        stage_message = f"Processing: {stage_info.get('label', current_stage)}"

        # Calculate progress
        if "stages" in workflow_state:
            completed_stages = sum(
                1 for s in workflow_state["stages"].values() if s.get("complete", False)
            )
            current_progress = (
                workflow_state["stages"].get(current_stage, {}).get("progress", 0)
            )
            overall_progress = ((completed_stages / len(WORKFLOW_STAGES)) * 100) + (
                current_progress / len(WORKFLOW_STAGES)
            )

    return stage_message, overall_progress, animated, striped, disable_interval


# Fix the callback to target the correct interval component
@callback(
    Output("classification-interval", "disabled", allow_duplicate=True),
    Input("classification-submit-button", "n_clicks"),
    prevent_initial_call=True,
)
def enable_interval(n_clicks):
    """Enable the interval when the workflow starts"""
    if n_clicks is None:
        raise PreventUpdate
    return False


# Add a function to create classification card
def create_classification_card(classification_data):
    """
    Create a card displaying classification results from multiple sources

    Args:
        classification_data: Dictionary containing classification information

    Returns:
        dmc.Paper component with classification information
    """
    if not classification_data:
        return None

    title = classification_data.get("title", "Classification Results")
    source = classification_data.get("source", "Unknown")
    family = classification_data.get("family")
    navis = classification_data.get("navis")
    haplotype = classification_data.get("haplotype")
    match_type = classification_data.get("match_type")
    confidence = classification_data.get("confidence", "Low")

    # Define badge colors based on source
    source_colors = {
        "exact_match": "green",
        "contained_match": "teal",
        "similar_match": "cyan",
        "blast_hit": "blue",
        "classification": "violet",
        "unknown": "gray",
    }

    # Define icon based on confidence level
    confidence_icons = {
        "High": "mdi:shield-check",
        "Medium": "mdi:shield-half-full",
        "Low": "mdi:shield-outline",
    }

    source_color = source_colors.get(source, "gray")
    confidence_icon = confidence_icons.get(confidence, "mdi:shield-outline")

    # Create list of classification details
    details = []

    if family:
        details.append(
            dmc.Group(
                [
                    dmc.Text("Family:", fw=700, size="lg"),
                    dmc.Text(family, size="lg", c="indigo"),
                ],
                pos="apart",
            )
        )

    if navis:
        details.append(
            dmc.Group(
                [dmc.Text("Navis:", fw=700), dmc.Text(navis, c="blue")], pos="apart"
            )
        )

    if haplotype:
        details.append(
            dmc.Group(
                [dmc.Text("Haplotype:", fw=700), dmc.Text(haplotype, c="violet")],
                pos="apart",
            )
        )

    if match_type:
        details.append(
            dmc.Group(
                [
                    dmc.Text("Match Type:", fw=700, size="sm"),
                    dmc.Text(match_type, c="dimmed", size="sm"),
                ],
                pos="apart",
            )
        )

    # Create card
    return dmc.Paper(
        children=[
            dmc.Group(
                [
                    dmc.Title(title, order=3),
                    dmc.Badge(source.replace("_", " ").title(), color=source_color),
                ],
                pos="apart",
            ),
            dmc.Space(h=10),
            *details,
            dmc.Space(h=10),
            dmc.Group(
                [
                    dmc.Text(f"Confidence: {confidence}", size="sm", c="dimmed"),
                    dmc.ThemeIcon(
                        DashIconify(icon=confidence_icon, width=16),
                        size="sm",
                        variant="light",
                        color=source_color,
                    ),
                ],
                pos="right",
            ),
        ],
        p="md",
        withBorder=True,
        shadow="sm",
        radius="md",
        style={"marginBottom": "1rem"},
    )


# Add a new callback to generate tabs based on multifasta data
@callback(
    Output("right-column-content", "children"),
    [
        Input("blast-sequences-store", "data"),
        Input("blast-results-store", "data"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def update_tabs_for_multifasta(seq_list, blast_results):
    # Add debug logging
    logger.info(
        f"update_tabs_for_multifasta called with multifasta_data: {type(seq_list)}, blast_results: {type(blast_results)}"
    )

    # If there's no multifasta data or blast results, don't update
    if not seq_list:
        logger.info("No multifasta_data, preventing update")
        raise PreventUpdate

    if not blast_results:
        logger.info("No blast_results, preventing update")
        raise PreventUpdate

    # If we have multifasta data (more than one sequence)
    if isinstance(seq_list, list) and len(seq_list) > 1:
        logger.info(f"Creating tabs for {len(seq_list)} sequences")

        # Create tabs for each sequence
        tabs = []
        for idx, seq in enumerate(seq_list[:10]):  # Limit to 10 sequences
            logger.info(f"Creating tab {idx}: {seq['header'][:20]}")
            tabs.append(
                dbc.Tab(
                    label=f"Sequence {idx + 1}: {seq['header'][:20]}{'...' if len(seq['header']) > 20 else ''}",
                    tab_id=f"tab-{idx}",
                )
            )

        logger.info(f"Created {len(tabs)} tab objects")

        # Create a simple Card with Tabs in the header
        output = dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Tabs(id="blast-tabs", active_tab="tab-0", children=tabs)
                ),
                dbc.CardBody(
                    [
                        # This div will contain the content for the active tab
                        html.Div(
                            id="tab-content",
                            children=[
                                # Initial content is for the first tab
                                html.Div(id="classification-output", className="mt-4"),
                                dmc.Stack(
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
                                ),
                                # BLAST results section
                                dmc.Stack(
                                    [
                                        html.Div(
                                            id="blast-container",
                                            style={
                                                "width": "100%",
                                                "display": "flex",
                                                "flexDirection": "column",
                                            },
                                        ),
                                        html.Div(id="blast-download"),
                                        html.Div(id="subject-seq-button"),
                                    ]
                                ),
                            ],
                        )
                    ]
                ),
            ]
        )

        logger.info("Successfully created card with tabs")
    else:
        # If there's only one sequence, return the default layout
        logger.info("No multifasta data or only one sequence, returning default layout")
        output = dmc.Stack(
            children=[
                # Progress section - initially hidden
                html.Div(id="classification-output", className="mt-4"),
                dmc.Stack(
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
                ),
                # BLAST results section
                dmc.Stack(
                    [
                        html.Div(
                            id="blast-container",
                            style={
                                "width": "100%",
                                "display": "flex",
                                "flexDirection": "column",
                            },
                        ),
                        html.Div(id="blast-download"),
                        html.Div(id="subject-seq-button"),
                    ]
                ),
            ],
        )
    return output


@callback(
    [
        Output("blast-active-tab", "data"),
        Output("blast-processed-tabs", "data"),
        Output("tab-content", "children"),
        Output(
            "blast-sequences-store", "data", allow_duplicate=True
        ),  # Add this to update the processed flag
    ],
    Input("blast-tabs", "active_tab"),
    [
        State("blast-processed-tabs", "data"),
        State("blast-sequences-store", "data"),
        State("blast-results-store", "data"),
        State("classification-result-store", "data"),
        State("evalue-threshold", "value"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_tab_switch(
    active_tab,
    processed_tabs,
    seq_list,
    blast_results,
    classification_result,
    evalue_threshold,
):
    """Track which tab is active and update processed tabs list"""
    logger.info(f"handle_tab_switch called with active_tab: {active_tab}")

    if active_tab is None or not seq_list:
        logger.info("No active_tab or multifasta_data, preventing update")
        raise PreventUpdate

    # Extract the numeric part from the tab ID (e.g., "tab-2" -> 2)
    tab_idx = int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0
    logger.info(f"Extracted tab index: {tab_idx}")

    # Add this tab to processed tabs list if not already there
    if processed_tabs is None:
        processed_tabs = []

    if tab_idx not in processed_tabs:
        logger.info(f"Adding tab {tab_idx} to processed tabs list")
        processed_tabs.append(tab_idx)

    # Get the sequence data for this tab
    if tab_idx >= len(seq_list):
        logger.error(
            f"Tab index {tab_idx} is out of range for multifasta_data length {len(seq_list)}"
        )
        raise PreventUpdate

    seq_data = seq_list[tab_idx]

    # Update the processed flag for this sequence
    updated_seq_list = seq_list.copy()
    updated_seq_list[tab_idx]["processed"] = True

    header = seq_data["header"]
    sequence = seq_data["sequence"]
    seq_type = seq_data.get("type", "nucl")
    logger.info(f"Processing sequence for tab {tab_idx}: {header[:30]}...")

    # Process this sequence for the tab
    # First write sequence to temporary FASTA file
    tmp_query_fasta = write_temp_fasta(header, sequence)
    logger.info(f"Wrote sequence to temporary file: {tmp_query_fasta}")

    # Run BLAST search for this sequence
    logger.info(f"Running BLAST search for tab {tab_idx}")
    blast_results_file = run_blast(
        query_type=seq_type,
        query_fasta=tmp_query_fasta,
        tmp_blast=tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name,
        input_eval=evalue_threshold,
        threads=2,
    )

    if blast_results_file is None:
        # Handle error case
        logger.error(f"No BLAST results were returned for tab {tab_idx}")
        classification_output = html.Div(
            create_error_alert("No BLAST results were returned")
        )
        blast_container = html.Div(
            [
                html.H2(
                    "BLAST Results",
                    style={
                        "marginTop": "15px",
                        "marginBottom": "20px",
                        "fontWeight": "500",
                        "textAlign": "left",
                    },
                ),
                create_error_alert("No BLAST results were returned"),
            ]
        )
        download_button = None
    else:
        logger.info(f"BLAST results returned: {blast_results_file}")
        # Process BLAST results for classification
        blast_tsv = parse_blast_xml(blast_results_file)
        blast_df = pd.read_csv(blast_tsv, sep="\t")
        logger.info(f"Parsed BLAST results with {len(blast_df)} hits")

        # Create classification output
        classification_data = None

        # Check if dataframe is empty
        if len(blast_df) == 0:
            logger.info("No BLAST hits for classification")
            classification_output = html.Div([create_no_matches_alert()])
        else:
            # Sort by evalue (ascending) and pident (descending) to get best hits
            blast_df = blast_df.sort_values(
                ["evalue", "pident"], ascending=[True, False]
            )
            top_hit = blast_df.iloc[0]
            logger.info(
                f"Top hit: {top_hit['hit_IDs']} with e-value {top_hit['evalue']}"
            )

            top_evalue = float(top_hit["evalue"])
            top_aln_length = int(top_hit["aln_length"])
            top_pident = float(top_hit["pident"])

            if top_pident >= 90:
                # look up family name from accession tag
                hit_IDs = top_hit["hit_IDs"]
                logger.info(f"Looking up metadata for hit ID: {hit_IDs}")
                meta_df = fetch_meta_data(accession_tag=hit_IDs)

                if not meta_df.empty:
                    logger.info(f"Metadata found for {hit_IDs}")
                    top_family = meta_df["familyName"].iloc[0]
                    navis = (
                        meta_df["starship_navis"].iloc[0]
                        if "starship_navis" in meta_df.columns
                        else None
                    )
                    haplotype = (
                        meta_df["starship_haplotype"].iloc[0]
                        if "starship_haplotype" in meta_df.columns
                        else None
                    )

                    classification_data = {
                        "title": "Classification from BLAST",
                        "source": "blast_hit",
                        "family": top_family,
                        "navis": navis,
                        "haplotype": haplotype,
                        "match_type": f"BLAST hit length {top_aln_length}bp with {top_pident:.1f}% identity to {hit_IDs}. Evalue: {top_evalue}",
                        "confidence": "Medium" if top_pident >= 95 else "Low",
                    }
                    logger.info(f"Created classification data: {classification_data}")
                else:
                    logger.info(f"No metadata found for hit ID: {hit_IDs}")

            # Create the classification card with title
            if classification_data:
                logger.info("Creating classification card")
                classification_output = html.Div(
                    [
                        dmc.Title(
                            "Classification Results",
                            order=2,
                            style={"marginBottom": "20px"},
                        ),
                        create_classification_card(classification_data),
                    ]
                )
            else:
                logger.info("No classification data available")
                # No classification available
                classification_output = html.Div(
                    [
                        dmc.Title(
                            "Classification Results",
                            order=2,
                            style={"marginBottom": "20px"},
                        ),
                        dmc.Alert(
                            title="No Classification Available",
                            children="Could not classify this sequence with any available method.",
                            color="yellow",
                            variant="light",
                        ),
                    ]
                )

        # Read BLAST output for visualization
        with open(blast_results_file, "r") as f:
            blast_results_text = f.read()
        logger.info(f"Read BLAST results file with {len(blast_results_text)} bytes")

        # Create the BlasterJS container (will be filled by clientside callback)
        blast_container = html.Div(
            id="blast-container",
            style={
                "width": "100%",
                "display": "flex",
                "flexDirection": "column",
            },
        )

        # Store the BLAST results for BlasterJS
        blast_data_store = dcc.Store(
            id="processed-blast-store", data={"blast_text": blast_results_text}
        )
        logger.info("Created BLAST data store")

        # Create download button
        download_button = blast_download_button()

    # Put together the tab content
    logger.info("Assembling complete tab content")
    tab_content = [
        # Classification output
        classification_output,
        # Progress section (hidden initially)
        dmc.Stack(
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
        ),
        # BLAST results section
        dmc.Stack(
            [
                blast_container,
                html.Div(id="blast-download", children=download_button),
                # Add the blast data store to trigger BlasterJS
                blast_data_store if blast_results_file is not None else None,
            ]
        ),
    ]

    logger.info(
        f"Returning data for tab {tab_idx} with content length {len(tab_content)}"
    )
    return tab_idx, processed_tabs, tab_content, updated_seq_list


@callback(
    Output("workflow-state-store", "data", allow_duplicate=True),
    Input("submit-button", "n_clicks"),
    State("blast-sequences-store", "data"),
    State("evalue-threshold", "value"),
    prevent_initial_call=True,
)
@handle_callback_error
def start_classification_workflow(n_clicks, seq_list, evalue_threshold):
    """Start the classification workflow and initialize the state store."""
    if not n_clicks or not seq_list:
        raise PreventUpdate

    # Prepare data for the workflow
    if isinstance(seq_list, list) and len(seq_list) > 0:
        seq_data = seq_list[0]  # Process the first sequence
        query_header = seq_data.get("header", "query")
        query_seq = seq_data.get("sequence", "")
        query_type = seq_data.get("type", "nucl")

        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)

        # Prepare upload data for classification
        upload_data = {
            "seq_type": query_type,
            "fasta": tmp_query_fasta,
            "fetch_ship_params": {
                "curated": False,
                "with_sequence": True,
                "dereplicate": True,
            },
            "fetch_captain_params": {"curated": True, "with_sequence": True},
        }

        # Start background task with celery
        classification_task = run_classification_workflow_task.delay(upload_data)

        # Return initial state with task ID
        return {
            "task_id": classification_task.id,
            "status": "started",
            "complete": False,
            "current_stage": None,
            "error": None,
            "start_time": time.time(),
        }

    # No valid sequences
    return None


@callback(
    Output("workflow-state-store", "data", allow_duplicate=True),
    Input("classification-interval", "n_intervals"),
    State("workflow-state-store", "data"),
    prevent_initial_call=True,
)
@handle_callback_error
def poll_classification_task(n_intervals, current_state):
    """Poll the classification task and update the state store."""
    if not current_state or not current_state.get("task_id"):
        raise PreventUpdate

    # Skip polling if task is already complete
    if current_state.get("complete", False):
        return dash.no_update

    # Get task status from celery
    task_id = current_state["task_id"]
    task = run_classification_workflow_task.AsyncResult(task_id)

    if task.state == "PENDING":
        # Task hasn't started yet
        return current_state
    elif task.state == "SUCCESS":
        # Task completed - return the full result
        # This contains the final state from run_classification_workflow
        return task.result
    elif task.state == "FAILURE":
        # Task failed
        return {
            **current_state,
            "complete": True,
            "error": str(task.result) if task.result else "Task failed",
        }

    # Still pending or unknown state
    return current_state

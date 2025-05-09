import dash
from dash import dcc, html, callback, clientside_callback
from dash_iconify import DashIconify
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import os
import base64
import pandas as pd
import logging
import time
import tempfile

from src.config.cache import cache
from src.utils.seq_utils import (
    check_input,
    write_temp_fasta,
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
            max_intervals=300,  # Maximum 5 minutes of polling
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
        console.log("BlasterJS callback triggered", data);
        
        // Check if we have valid data
        if (!data || !data.blast_text) {
            console.log("No valid BLAST data available");
            // Return empty div to avoid undefined errors
            const container = document.getElementById('blast-container');
            return window.dash_clientside.no_update;
        }
        
        try {
            console.log("Initializing BlasterJS with text length:", data.blast_text.length);
            
            // Use standard ID selector
            const container = document.getElementById('blast-container');
            const loader = document.getElementById('blast-loader');
            
            if (!container) {
                console.error("No blast container found in the DOM");
                console.log("Available IDs:", 
                    Array.from(document.querySelectorAll('[id]'))
                        .map(el => el.id)
                        .join(', ')
                );
                return window.dash_clientside.no_update;
            }
            
            console.log("Found container:", container);
            
            // Clear existing content first
            container.innerHTML = '';
            
            // If we have empty blast text, show a message
            if (!data.blast_text.trim()) {
                container.innerHTML = '<div style="padding: 20px; text-align: center; color: #666;">No BLAST results to display</div>';
                return window.dash_clientside.no_update;
            }
            
            // Create the title element
            const titleElement = document.createElement('h2');
            titleElement.innerHTML = 'BLAST Results';
            titleElement.style.marginTop = '15px';
            titleElement.style.marginBottom = '20px';
            titleElement.style.textAlign = 'left';
            titleElement.style.width = '100%';
            
            // Create simple divs for BlasterJS with explicit left alignment
            const alignmentsDiv = document.createElement('div');
            alignmentsDiv.id = 'blast-multiple-alignments';
            alignmentsDiv.style.textAlign = 'left';
            alignmentsDiv.style.width = '100%';
            
            const tableDiv = document.createElement('div');
            tableDiv.id = 'blast-alignments-table';
            tableDiv.style.textAlign = 'left';
            tableDiv.style.width = '100%';
            
            // Add elements to the container
            container.appendChild(titleElement);
            container.appendChild(alignmentsDiv);
            container.appendChild(tableDiv);
            
            // Add a style element to ensure BlasterJS output is properly aligned
            const styleEl = document.createElement('style');
            styleEl.textContent = `
                #blast-multiple-alignments, #blast-alignments-table {
                    text-align: left !important;
                    margin-left: 0 !important;
                    padding-left: 0 !important;
                }
                #blast-multiple-alignments div, #blast-alignments-table div,
                #blast-multiple-alignments table, #blast-alignments-table table {
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
                    multipleAlignments: "blast-multiple-alignments",
                    alignmentsTable: "blast-alignments-table"
                });
                console.log("BlasterJS initialized successfully");
                
                // Additional styling fix after BlasterJS renders
                setTimeout(function() {
                    const tables = container.querySelectorAll('table');
                    tables.forEach(function(table) {
                        table.style.marginLeft = '0';
                        table.style.textAlign = 'left';
                    });
                }, 100);
                
                // Hide the loader after successful initialization
                if (loader) {
                    // Add a small delay to ensure content is rendered
                    setTimeout(function() {
                        // Access the data-dashloaderisloading attribute and set it to false
                        loader.setAttribute('data-dashloaderisloading', 'false');
                    }, 200);
                }
            } catch (blasterError) {
                console.error("Error initializing BlasterJS library:", blasterError);
                container.innerHTML += "<div style='color:red;'>Error initializing BLAST viewer: " + blasterError + "</div>";
                
                // Hide the loader even if there's an error
                if (loader) {
                    loader.setAttribute('data-dashloaderisloading', 'false');
                }
            }
            
            return window.dash_clientside.no_update;
        } catch (error) {
            console.error('Overall error in callback:', error);
            // Error handling - display error in container
            const container = document.getElementById('blast-container');
            if (container) {
                container.innerHTML = "<div style='color:red; padding: 20px;'>Error displaying BLAST results: " + error.toString() + "</div>";
            }
            return window.dash_clientside.no_update;
        }
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
        Output("right-column-content", "children", allow_duplicate=True),
    ],
    [
        Input("submit-button", "n_clicks"),
    ],
    [
        State("query-text", "value"),
        State("blast-sequences-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
    id="preprocess-callback",
)
def preprocess(n_clicks, query_text_input, seq_list):
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
        # If we have parsed sequences from a file upload, use those
        if seq_list is not None:
            logger.info(
                f"Using pre-parsed sequences from file upload: {len(seq_list)} sequences"
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

            return seq_list, None, None, ui_content

        # Otherwise, parse the text input with max_sequences=10
        if query_text_input:
            logger.info("Processing text input")
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
                    None,
                )

            logger.info(f"Text input processed: {input_type}, {n_seqs} sequences")

            # Create UI structure based on number of sequences from text input
            ui_content = (
                create_tab_layout(seq_list)
                if len(seq_list) > 1
                else create_single_layout()
            )

            return seq_list, None, None, ui_content

        # If we have neither text input nor uploaded sequences, show an error
        error_alert = dmc.Alert(
            title="No Input Provided",
            children="Please enter a sequence or upload a FASTA file.",
            color="yellow",
            variant="filled",
        )
        return None, "No input provided", error_alert, None

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        error_alert = dmc.Alert(
            title="Error Processing Input",
            children=f"An unexpected error occurred: {str(e)}",
            color="red",
            variant="filled",
        )
        return None, str(e), error_alert, None


# Helper functions to create layouts
def create_single_layout():
    return dmc.Stack(
        children=[
            html.Div(id="classification-output", className="mt-4"),
            # Progress section
            dmc.Stack(
                [
                    dmc.Group([dbc.Progress(id="classification-progress", value=0)]),
                    dmc.Group(
                        [
                            dmc.Text("Classification Status:", size="lg", fw=500),
                            dmc.Text(
                                id="classification-stage-display", size="lg", c="blue"
                            ),
                        ]
                    ),
                ],
                id="classification-progress-section",
                style={"display": "none"},
            ),
            # BLAST results section
            dmc.Stack([html.Div(id="blast-container")]),
        ],
    )


def create_tab_layout(seq_list):
    # Create tabs for multiple sequences
    tabs = []
    for idx, seq in enumerate(seq_list[:10]):  # Limit to 10 sequences
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
    [Input("blast-results-store", "data"), Input("blast-active-tab", "data")],
)
def process_blast_results(blast_results_store, active_tab_idx):
    if not blast_results_store:
        logger.warning("No blast_results_store data")
        return None

    # Convert active_tab_idx to string (since keys in sequence_results are strings)
    tab_idx = str(active_tab_idx or 0)
    logger.info(f"Processing BLAST results for tab index: {tab_idx}")

    # Get the correct sequence results based on active tab
    sequence_results = blast_results_store.get("sequence_results", {}).get(tab_idx)
    if not sequence_results:
        logger.warning(f"No sequence results for tab index {tab_idx}")
        return None

    blast_results_file = sequence_results.get("blast_file")
    if not blast_results_file:
        logger.warning(f"No blast file in sequence results for tab {tab_idx}")
        return None

    logger.info(f"Reading BLAST file: {blast_results_file}")
    try:
        # Check if file exists
        if not os.path.exists(blast_results_file):
            logger.error(f"BLAST file does not exist: {blast_results_file}")
            return None

        # Read the BLAST output as text
        with open(blast_results_file, "r") as f:
            blast_results = f.read()

        # Add size limit check
        results_size = len(blast_results)
        logger.info(f"Read BLAST results, size: {results_size} bytes")
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
        # Always return consistent data format, even on error
        return {"blast_text": ""}  # Return empty string instead of None


@callback(
    [
        Output("blast-results-store", "data"),
        Output("classification-interval", "disabled", allow_duplicate=True),
        Output("workflow-state-store", "data"),
    ],
    [
        Input("submission-id-store", "data"),
    ],
    [
        State("blast-sequences-store", "data"),
        State("evalue-threshold", "value"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
def process_sequences(submission_id, seq_list, evalue_threshold):
    """Process only the first sequence initially, for both text input and file uploads"""
    if not all([seq_list, submission_id]):
        return None, True, None

    try:
        blast_results = None
        classification_interval_disabled = True
        workflow_state = None

        # Process only the first sequence to start
        if seq_list and len(seq_list) > 0:
            # Process sequence 0
            sequence_result = process_single_sequence(seq_list[0], evalue_threshold)

            if sequence_result:
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

                # Default to disabled
                classification_interval_disabled = True

                # Set up workflow state if needed
                if not skip_classification and sequence_result.get("blast_file"):
                    # Create upload_data with file contents instead of file paths
                    # Read the file content to include in upload_data
                    fasta_content = None
                    tmp_query_fasta = None

                    try:
                        tmp_query_fasta = write_temp_fasta(
                            seq_list[0].get("header", "query"),
                            seq_list[0].get("sequence", ""),
                        )

                        if tmp_query_fasta and os.path.exists(tmp_query_fasta):
                            with open(tmp_query_fasta, "r") as f:
                                fasta_content = f.read()
                    except Exception as e:
                        logger.error(f"Error reading temp FASTA file: {e}")
                        # If we can't read the file, skip classification
                        skip_classification = True

                    if not skip_classification:
                        # Initialize a complete workflow state structure
                        workflow_state = {
                            "task_id": str(submission_id),  # Ensure it's a string
                            "status": "initialized",
                            "complete": False,
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
                            "upload_data": {
                                "seq_type": seq_list[0].get("type", "nucl"),
                                "fasta": {"content": fasta_content}
                                if fasta_content
                                else (tmp_query_fasta or ""),
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

                        # Only enable the interval if we have a proper workflow state
                        if workflow_state is not None:
                            classification_interval_disabled = False

        return blast_results, classification_interval_disabled, workflow_state
    except Exception as e:
        logger.error(f"Error in process_sequences: {e}")
        # Return basic data on error
        error_state = {
            "complete": True,
            "error": str(e),
            "status": "failed",
            "task_id": str(submission_id) if submission_id else None,
        }
        return None, True, error_state


def process_single_sequence(seq_data, evalue_threshold):
    """Process a single sequence and return structured results"""
    query_header = seq_data.get("header", "query")
    query_seq = seq_data.get("sequence", "")
    query_type = seq_data.get("type", "nucl")

    # Write sequence to temporary FASTA file
    tmp_query_fasta = write_temp_fasta(query_header, query_seq)
    if not tmp_query_fasta:
        logger.error("Failed to write temporary FASTA file")
        return {
            "blast_file": None,
            "classification": None,
            "processed": False,
            "sequence": query_seq,
            "error": "Failed to create temporary file",
        }

    try:
        # Run BLAST search
        blast_task = run_blast_search_task.delay(
            query_header=query_header,
            query_seq=query_seq,
            query_type=query_type,
            eval_threshold=evalue_threshold,
        )
        blast_results = blast_task.get(timeout=300)

        # Handle case where blast_results is None or invalid
        if not blast_results:
            logger.warning("BLAST search returned no results")
            blast_results_file = None
        elif isinstance(blast_results, dict) and "content" in blast_results:
            # Handle the new format (content-based instead of file path)
            blast_results_file = tempfile.NamedTemporaryFile(
                suffix=".blast", delete=False
            ).name
            with open(blast_results_file, "w") as f:
                f.write(blast_results["content"])
            logger.debug(f"Created temporary BLAST results file: {blast_results_file}")
        elif isinstance(blast_results, str) and os.path.exists(blast_results):
            # Handle the old format (file path)
            blast_results_file = blast_results
            logger.debug(f"Using existing BLAST results file: {blast_results_file}")
        else:
            # Handle unexpected format
            logger.error(f"Unexpected BLAST results format: {type(blast_results)}")
            blast_results_file = None
    except Exception as e:
        logger.error(f"Error running BLAST search: {e}")
        blast_results_file = None

    # If BLAST search failed, return basic structure with error
    if not blast_results_file:
        return {
            "blast_file": None,
            "classification": None,
            "processed": False,
            "sequence": query_seq,
            "error": "BLAST search failed",
        }

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

    # Initialize classification data
    classification_data = None

    try:
        # Run classification workflow
        classification_task = run_classification_workflow_task.delay(
            upload_data=upload_data,
        )
        workflow_result = classification_task.get(timeout=300)

        # Check for classification result data first
        if workflow_result is not None:
            result_type = workflow_result.get("type")

            if result_type == "match" and workflow_result.get("match_stage") in [
                "exact",
                "contained",
                "similar",
            ]:
                # Handle exact/contained/similar matches
                match_stage = workflow_result.get("match_stage")
                match_result = workflow_result.get("match_result")

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
            elif workflow_result and workflow_result.get("complete", False):
                # Process based on classification state (existing logic)

                # Check for workflow-based classification (highest priority)
                if workflow_result.get("match_stage") in [
                    "exact",
                    "contained",
                    "similar",
                ]:
                    match_stage = workflow_result.get("match_stage")
                    match_result = workflow_result.get("match_result")

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
                elif workflow_result.get("found_match", False):
                    match_stage = workflow_result.get("match_stage")
                    match_result = workflow_result.get("match_result")

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
    except Exception as e:
        logger.error(f"Error in classification workflow: {e}")

    # If no classification from workflow, fall back to BLAST results
    if classification_data is None:
        # Process the BLAST results
        try:
            # Verify the file exists before processing
            if blast_results_file and os.path.exists(blast_results_file):
                blast_tsv = parse_blast_xml(blast_results_file)
                if blast_tsv and os.path.exists(blast_tsv):
                    blast_df = pd.read_csv(blast_tsv, sep="\t")

                    # Check if dataframe is empty
                    if len(blast_df) == 0:
                        logger.info("No BLAST hits found")
                    else:
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
                                    "confidence": "Medium"
                                    if top_pident >= 95
                                    else "Low",
                                }
        except Exception as e:
            logger.error(f"Error processing BLAST results: {e}")

    # If still no classification, try HMMER
    if classification_data is None:
        try:
            hmmer_task = run_hmmer_search_task.delay(
                query_header=query_header,
                query_seq=query_seq,
                query_type=query_type,
                eval_threshold=evalue_threshold,
            )
            hmmer_results = hmmer_task.get(timeout=300)

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
                captain_family, aln_length, evalue = select_ship_family(first_result)

                if captain_family:
                    confidence = "High" if aln_length > 80 and evalue < 0.001 else "Low"

                    classification_data = {
                        "title": "Captain Gene Classification",
                        "source": "hmmsearch",
                        "family": captain_family,
                        "match_type": f"Captain gene match with length {aln_length}, Evalue: {evalue}",
                        "confidence": confidence,
                    }
        except Exception as e:
            logger.error(f"Error processing HMMER results: {e}")

    # Return structured results
    return {
        "blast_file": blast_results_file,
        "classification": classification_data,
        "processed": True,
        "sequence": query_seq,  # Include sequence for length checks
    }


# 2. Add a callback to process additional sequences when tabs are clicked
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
    if not active_tab or not results_store or not seq_list:
        raise PreventUpdate

    # Extract tab index from tab id (e.g., "tab-2" -> 2)
    tab_idx = int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0

    # Check if this sequence has already been processed
    processed_sequences = results_store.get("processed_sequences", [])

    if tab_idx in processed_sequences:
        # Already processed this sequence
        raise PreventUpdate

    # Ensure tab_idx is valid
    if tab_idx >= len(seq_list):
        logger.error(
            f"Tab index {tab_idx} out of range (seq_list length: {len(seq_list)})"
        )
        raise PreventUpdate

    # Process this sequence
    logger.info(f"Processing sequence for tab {tab_idx}")
    sequence_result = process_single_sequence(seq_list[tab_idx], evalue_threshold)

    if not sequence_result:
        logger.error(f"Failed to process sequence for tab {tab_idx}")
        raise PreventUpdate

    # Update the results store
    updated_results = dict(results_store)
    updated_results["processed_sequences"].append(tab_idx)
    updated_results["sequence_results"][str(tab_idx)] = sequence_result

    logger.info(f"Updated results for tab {tab_idx}")
    return updated_results


# 4. Render the content for the current tab
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
    if not active_tab or not results_store:
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
    if tab_idx not in processed_sequences:
        # Not processed yet - show loading state
        return dmc.Stack(
            [
                dmc.Center(dmc.Loader(size="xl")),
                dmc.Text("Processing sequence...", size="lg"),
            ]
        )

    # Get results for this sequence
    sequence_results = results_store.get("sequence_results", {}).get(str(tab_idx))
    if not sequence_results:
        return dmc.Alert(
            title="Error",
            children="No results found for this sequence",
            color="red",
            variant="filled",
        )

    # Check for errors in the results
    if "error" in sequence_results:
        return dmc.Alert(
            title="Error Processing Sequence",
            children=sequence_results["error"],
            color="red",
            variant="filled",
        )

    # Render classification results
    classification_output = create_classification_output(sequence_results)

    # Render BLAST results
    blast_container = create_blast_container(sequence_results)

    return [
        classification_output,
        # Progress section (hidden initially)
        dmc.Stack(
            [
                dmc.Group([dbc.Progress(id="classification-progress", value=0)]),
                dmc.Group(
                    [
                        dmc.Text("Classification Status:", size="lg", fw=500),
                        dmc.Text(
                            id="classification-stage-display", size="lg", c="blue"
                        ),
                    ]
                ),
            ],
            id="classification-progress-section",
            style={"display": "none"},
        ),
        # BLAST results section
        dmc.Stack(
            [
                blast_container,
            ]
        ),
    ]


def create_classification_output(sequence_results):
    """Create the classification output component"""
    # Extract classification data
    classification_data = sequence_results.get("classification")

    if classification_data:
        return html.Div(
            [
                dmc.Title(
                    "Classification Results", order=2, style={"marginBottom": "20px"}
                ),
                create_classification_card(classification_data),
            ]
        )
    else:
        # No classification available
        return html.Div(
            [
                dmc.Title(
                    "Classification Results", order=2, style={"marginBottom": "20px"}
                ),
                dmc.Alert(
                    title="No Classification Available",
                    children="Could not classify this sequence with any available method.",
                    color="yellow",
                    variant="light",
                ),
            ]
        )


def create_blast_container(sequence_results):
    """Create the BLAST container with a spinner wrapper"""
    blast_file = sequence_results.get("blast_file")
    if not blast_file:
        return html.Div(
            [
                html.H2(
                    "BLAST Results", style={"marginTop": "15px", "marginBottom": "20px"}
                ),
                create_no_matches_alert(),
            ]
        )

    # Create with a regular string ID and wrap in a spinner
    return html.Div(
        [
            dbc.Spinner(
                children=html.Div(
                    id="blast-container",  # Use regular string ID
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


def load_blast_text(blast_file):
    """Load BLAST output as text"""
    if not blast_file or not os.path.exists(blast_file):
        return None

    with open(blast_file, "r") as f:
        return f.read()


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


@callback(
    Output("classification-output", "children"),
    Input("blast-results-store", "data"),
    prevent_initial_call=True,
)
def update_single_sequence_classification(blast_results_store):
    """Update classification output for single sequence submissions"""
    if not blast_results_store:
        return None

    # Get the results for sequence 0
    sequence_results = blast_results_store.get("sequence_results", {}).get("0")
    if not sequence_results:
        return None

    # Create the classification output using the existing function
    return create_classification_output(sequence_results)


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
    Output("workflow-state-store", "data", allow_duplicate=True),
    [
        Input("classification-interval", "n_intervals"),
        Input("workflow-state-store", "data"),
    ],
    prevent_initial_call=True,
)
def update_classification_workflow_state(n_intervals, workflow_state):
    """Poll for updated classification workflow state"""
    # First check if we should skip this update
    if workflow_state is None:
        logger.warning("Workflow state is None")
        return {"complete": True, "error": "Missing workflow state", "status": "failed"}

    if workflow_state.get("complete", False):
        raise PreventUpdate

    # Check for timing out after too many intervals (e.g., 300 = 5 minutes with 1s interval)
    max_intervals = 300  # 5 minutes at 1 second intervals
    if n_intervals and n_intervals > max_intervals:
        logger.warning(
            f"Classification polling timed out after {n_intervals} intervals"
        )
        workflow_state = (
            workflow_state.copy()
        )  # Create a copy to avoid modifying the input
        workflow_state["error"] = "Classification timed out"
        workflow_state["complete"] = True
        workflow_state["status"] = "timeout"
        return workflow_state

    # Check if task is running
    task_id = workflow_state.get("task_id")
    if not task_id:
        logger.warning("No task ID found in workflow state")
        workflow_state = workflow_state.copy()
        workflow_state["error"] = "Missing task ID"
        workflow_state["complete"] = True
        workflow_state["status"] = "failed"
        return workflow_state

    try:
        # Start the classification workflow task if needed
        if workflow_state.get("task_id") and not workflow_state.get("celery_task_id"):
            # Use the upload_data from the workflow_state
            upload_data = workflow_state.get("upload_data", {})
            if not upload_data:
                logger.error("No upload_data found in workflow state")
                workflow_state = workflow_state.copy()
                workflow_state["error"] = "Missing upload data"
                workflow_state["complete"] = True
                workflow_state["status"] = "failed"
                return workflow_state

            # Start the task
            task = run_classification_workflow_task.delay(upload_data=upload_data)
            workflow_state = workflow_state.copy()
            workflow_state["celery_task_id"] = task.id
            # Initialize stages if not present
            if "stages" not in workflow_state:
                workflow_state["stages"] = {
                    stage["id"]: {"progress": 0, "complete": False}
                    for stage in WORKFLOW_STAGES
                }
            return workflow_state

        # Check existing task status
        celery_task_id = workflow_state.get("celery_task_id")
        if not celery_task_id:
            logger.warning("No Celery task ID found in workflow state")
            workflow_state = workflow_state.copy()
            workflow_state["error"] = "Missing Celery task ID"
            workflow_state["complete"] = True
            workflow_state["status"] = "failed"
            return workflow_state

        # Get AsyncResult for task
        task = run_classification_workflow_task.AsyncResult(celery_task_id)
        task_status = task.status

        # Create a copy of workflow state before modifying
        workflow_state = workflow_state.copy()

        # Update workflow status based on task state
        if task_status == "PENDING":
            workflow_state["status"] = "pending"
        elif task_status == "STARTED":
            workflow_state["status"] = "running"
        elif task_status == "SUCCESS":
            # Merge the results from the task
            try:
                task_result = task.get(timeout=5)  # Add timeout to prevent hanging
                if isinstance(task_result, dict):
                    # Make sure we don't overwrite crucial fields with None values
                    for key, value in task_result.items():
                        if value is not None or key not in workflow_state:
                            workflow_state[key] = value
            except Exception as e:
                logger.error(f"Error getting task result: {e}")

            workflow_state["status"] = "complete"
            workflow_state["complete"] = True
        elif task_status in ("FAILURE", "REVOKED"):
            workflow_state["status"] = "failed"
            try:
                workflow_state["error"] = (
                    str(task.result) if task.result else "Task failed"
                )
            except Exception as e:
                workflow_state["error"] = f"Error accessing task result: {str(e)}"
            workflow_state["complete"] = True
        else:
            # For unknown status, log it but keep current state
            logger.warning(f"Unknown task status: {task_status}")
            workflow_state["status"] = f"unknown ({task_status})"

        return workflow_state

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
        return workflow_state


@callback(
    [
        Output("classification-progress", "value"),
        Output("classification-stage-display", "children"),
        Output("classification-progress-section", "style"),
    ],
    Input("workflow-state-store", "data"),
    prevent_initial_call=True,
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

    # Show section if classification is in progress
    style = (
        {"display": "block"} if status and status != "pending" else {"display": "none"}
    )

    return progress, stage_text, style


@callback(
    Output("classification-interval", "disabled", allow_duplicate=True),
    Input("workflow-state-store", "data"),
    prevent_initial_call=True,
)
def disable_interval_when_complete(workflow_state):
    """Disable the interval when classification is complete"""
    if workflow_state is None:
        # If no workflow state, keep the interval disabled
        return True

    # Always disable if complete flag is set
    if workflow_state.get("complete", False):
        return True

    # Disable if there's an error
    if workflow_state.get("error") is not None:
        return True

    # Disable if status indicates completion
    if workflow_state.get("status") in ["complete", "failed", "timeout"]:
        return True

    # Otherwise, keep the interval running
    return False

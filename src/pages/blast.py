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

from src.config.cache import cache
from src.utils.seq_utils import (
    check_input,
    write_temp_fasta,
)
from src.utils.blast_utils import (
    create_no_matches_alert,
    parse_blast_xml,
    blast_download_button,
    select_ship_family,
)

from src.components.callbacks import (
    curated_switch,
    create_file_upload,
)
from src.database.sql_manager import fetch_meta_data
from src.utils.telemetry import blast_limit_decorator
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
        console.log("BlasterJS callback triggered", data);
        if (data && data.blast_text) {
            try {
                console.log("Initializing BlasterJS with text length:", data.blast_text.length);
                
                // Use standard ID selector
                const container = document.getElementById('blast-container');
                
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
                
                // Create the title element
                const titleElement = document.createElement('h2');
                titleElement.innerHTML = 'BLAST Results';
                titleElement.style.marginTop = '15px';
                titleElement.style.marginBottom = '20px';
                
                // Create simple divs for BlasterJS
                const alignmentsDiv = document.createElement('div');
                alignmentsDiv.id = 'blast-multiple-alignments';
                
                const tableDiv = document.createElement('div');
                tableDiv.id = 'blast-alignments-table';
                
                // Add elements to the container
                container.appendChild(titleElement);
                container.appendChild(alignmentsDiv);
                container.appendChild(tableDiv);
                
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
                } catch (blasterError) {
                    console.error("Error initializing BlasterJS library:", blasterError);
                    container.innerHTML += "<div style='color:red;'>Error initializing BLAST viewer: " + blasterError + "</div>";
                }
                
                return window.dash_clientside.no_update;
            } catch (error) {
                console.error('Overall error in callback:', error);
                return error.toString();
            }
        }
        console.log("No BLAST data available");
        return window.dash_clientside.no_update;
    }
    """,
    Output("blast-container", "children"),  # Use standard ID here
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
        return None


@blast_limit_decorator
@callback(
    [
        Output("blast-results-store", "data"),
        Output("classification-interval", "disabled"),
        Output("workflow-state-store", "data"),
    ],
    [
        Input("submission-id-store", "data"),
        Input("blast-sequences-store", "data"),
    ],
    [
        State("evalue-threshold", "value"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
def process_sequences(submission_id, seq_list, evalue_threshold):
    """Process only the first sequence initially"""
    if not all([seq_list, submission_id]):
        return None, None, None

    blast_results = None
    classification_interval_disabled = None
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
            }

            # Determine if we should enable classification interval
            skip_classification = len(sequence_result.get("sequence", "")) < 5000
            classification_interval_disabled = skip_classification

            # Set up workflow state if needed
            if not skip_classification and sequence_result.get("blast_file"):
                # This should be set up for classification workflow
                workflow_state = {
                    "task_id": submission_id,  # Use submission ID for now
                    "status": "started",
                    "complete": False,
                    "current_stage": None,
                    "error": None,
                    "start_time": time.time(),
                }

    return blast_results, classification_interval_disabled, workflow_state


def process_single_sequence(seq_data, evalue_threshold):
    """Process a single sequence and return structured results"""
    query_header = seq_data.get("header", "query")
    query_seq = seq_data.get("sequence", "")
    query_type = seq_data.get("type", "nucl")

    # Write sequence to temporary FASTA file
    tmp_query_fasta = write_temp_fasta(query_header, query_seq)

    # Run BLAST search
    blast_task = run_blast_search_task.delay(
        query_header=query_header,
        query_seq=query_seq,
        query_type=query_type,
        eval_threshold=evalue_threshold,
    )
    blast_results_file = blast_task.get(timeout=300)

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

    # Run classification workflow
    classification_task = run_classification_workflow_task.delay(
        upload_data=upload_data,
    )

    # Initialize classification data
    classification_data = None
    try:
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
            blast_tsv = parse_blast_xml(blast_results_file)
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
                            "confidence": "Medium" if top_pident >= 95 else "Low",
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
            captain_results = hmmer_results.get("results") if hmmer_results else None

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

    # Process this sequence
    if tab_idx < len(seq_list):
        sequence_result = process_single_sequence(seq_list[tab_idx], evalue_threshold)

        # Update the results store
        updated_results = dict(results_store)
        updated_results["processed_sequences"].append(tab_idx)
        updated_results["sequence_results"][str(tab_idx)] = sequence_result

        return updated_results

    raise PreventUpdate


# 3. Create tabs based on number of sequences
@callback(
    Output("right-column-content", "children"),
    Input("blast-sequences-store", "data"),
    prevent_initial_call=True,
)
def create_tabs_structure(seq_list):
    """Create tabs structure based on number of sequences"""
    if not seq_list:
        raise PreventUpdate

    # Single sequence case - no tabs needed
    if len(seq_list) == 1:
        return dmc.Stack(
            children=[
                html.Div(id="classification-output", className="mt-4"),
                # Progress section
                dmc.Stack(
                    [
                        dmc.Group(
                            [dbc.Progress(id="classification-progress", value=0)]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Classification Status:", size="lg", fw=500),
                                dmc.Text(
                                    id="classification-stage-display",
                                    size="lg",
                                    c="blue",
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
                        html.Div(id="blast-container"),
                        html.Div(id="blast-download"),
                    ]
                ),
            ],
        )

    # Multiple sequences - create tabs
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
    tab_idx = int(active_tab.split("-")[1]) if active_tab and "-" in active_tab else 0

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

    # Render classification results
    classification_output = create_classification_output(sequence_results)

    # Render BLAST results
    blast_container = create_blast_container(sequence_results)

    # Render download button
    download_button = (
        blast_download_button() if sequence_results.get("blast_file") else None
    )

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
                html.Div(id="blast-download", children=download_button),
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
    """Create the BLAST container"""
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

    # Create with a regular string ID
    return html.Div(
        id="blast-container",  # Use regular string ID
        style={
            "width": "100%",
            "display": "flex",
            "flexDirection": "column",
        },
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
    Output("right-column-content-debug", "children"),
    Input("processed-blast-store", "data"),
    prevent_initial_call=True,
)
def debug_blast_container(processed_blast_data):
    if processed_blast_data:
        blast_text_size = (
            len(processed_blast_data.get("blast_text", ""))
            if processed_blast_data.get("blast_text")
            else 0
        )
        return f"BLAST data size: {blast_text_size} bytes"
    return "No BLAST data"

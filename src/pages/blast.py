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
from threading import Thread

from src.config.cache import cache
from src.utils.seq_utils import (
    guess_seq_type,
    check_input,
    write_temp_fasta,
    parse_fasta,
    parse_fasta_from_file,
    parse_multifasta_from_file,
    write_fasta,
)
from src.utils.blast_utils import (
    run_blast,
    run_hmmer,
    run_diamond,
    process_captain_results,
    create_no_matches_alert,
    parse_blast_xml,
    blast_download_button,
    get_blast_db,
)

from src.components.callbacks import (
    curated_switch,
    create_file_upload,
    create_accession_modal,
)
from src.database.sql_manager import fetch_meta_data, fetch_captains, fetch_ships
from src.config.settings import BLAST_DB_PATHS
from src.utils.telemetry import blast_limit_decorator
from src.components.error_boundary import handle_callback_error, create_error_alert
from src.utils.classification_utils import WORKFLOW_STAGES

from src.tasks import (
    check_exact_matches_task,
    check_contained_matches_task,
    check_similar_matches_task,
    run_family_classification_task,
    run_navis_classification_task,
    run_haplotype_classification_task,
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
        dcc.Store(id="query-header-store"),
        dcc.Store(id="query-seq-store"),
        dcc.Store(id="query-type-store"),
        dcc.Store(id="upload-error-store"),
        # Results stores
        dcc.Store(id="blast-results-store"),
        dcc.Store(id="captain-results-store"),
        dcc.Store(id="submission-id-store"),
        # Processed data stores
        dcc.Store(id="processed-metadata-store"),
        dcc.Store(id="processed-blast-store"),
        # Multifasta stores
        dcc.Store(id="multifasta-sequences-store"),  # Store list of sequences
        dcc.Store(id="multifasta-active-tab", data=0),  # Store active tab index
        dcc.Store(
            id="multifasta-processed-tabs", data=[0]
        ),  # Store which tabs have been processed
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
        dcc.Store(id="classification-result-store", data=None),
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
                                                            output_id="blast-fasta-sequence-upload",
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
                                # This div will be replaced with tabs when multifasta is detected
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
                                                html.Div(id="subject-seq-button"),
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
    Output("blast-fasta-sequence-upload", "children"),
    [
        Input("blast-fasta-upload", "contents"),
        Input("blast-fasta-upload", "filename"),
    ],
)
@handle_callback_error
def update_fasta_details(seq_content, seq_filename):
    if seq_content is None:
        return [
            html.Div(
                html.P(
                    ["Select a FASTA file to upload"],
                )
            )
        ]
    else:
        try:
            # "," is the delimeter for splitting content_type from content_string
            content_type, content_string = seq_content.split(",")
            query_string = base64.b64decode(content_string).decode("utf-8")

            # Try to parse as multifasta
            sequences, seq_count, error = parse_multifasta_from_file(seq_content)

            if error and not sequences:
                # If parse as multifasta failed, show the error
                return html.Div([error])

            if seq_count > 1:
                # Multiple sequences found - display count and first few
                return html.Div(
                    [
                        html.P(f"File name: {seq_filename}"),
                        html.P(
                            [
                                f"Found {seq_count} sequences in the file",
                                html.Span(
                                    " (maximum of 10 will be processed)",
                                    style={"color": "orange", "fontStyle": "italic"}
                                    if seq_count > 10
                                    else {"display": "none"},
                                ),
                            ]
                        ),
                        html.Ul(
                            [
                                html.Li(f"{seq[0]} ({len(seq[1])} bp)")
                                for seq in sequences[:3]
                            ]
                            + ([html.Li("...")] if seq_count > 3 else [])
                        ),
                        html.P(
                            "Only the first 10 sequences will be processed."
                            if seq_count > 10
                            else "All sequences will be processed when you submit."
                        ),
                    ]
                )
            else:
                # Single sequence (backward compatible)
                children = parse_fasta(query_string, seq_filename)
                return children

        except Exception as e:
            logger.error(e)
            return html.Div(["There was an error processing this file."])


@callback(
    [
        Output("submit-button", "disabled"),
        Output("submit-button", "children"),
        Output("upload-error-message", "children"),
        Output("upload-error-store", "data"),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("blast-fasta-upload", "contents"),
        Input("blast-timeout-interval", "n_intervals"),
    ],
    [State("blast-fasta-upload", "filename"), State("blast-timeout-store", "data")],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_submission_and_upload(
    n_clicks, contents, n_intervals, filename, timeout_triggered
):
    triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

    # Default return values
    button_disabled = False
    button_text = "Submit BLAST"
    error_message = ""
    error_store = None

    # Handle timeout case
    if triggered_id == "blast-timeout-interval":
        if n_clicks and timeout_triggered:
            button_disabled = True
            button_text = "Server timeout"
            error_message = "The server is taking longer than expected to respond. Please try again later."
            error_store = error_message

    # Handle file upload
    if triggered_id == "blast-fasta-upload" and contents is not None:
        max_size = 10 * 1024 * 1024  # 10 MB
        content_type, content_string = contents.split(",")

        # Use our updated parse_fasta_from_file function
        header, seq, fasta_error = parse_fasta_from_file(contents)

        decoded = base64.b64decode(content_string)
        file_size = len(decoded)

        if fasta_error:
            error_message = fasta_error
            error_store = error_message
            button_disabled = True
        elif file_size > max_size:
            error_message = f"Error: The file '{filename}' exceeds the 10 MB limit."
            error_store = error_message
            button_disabled = True

    return [button_disabled, button_text, error_message, error_store]


@blast_limit_decorator
@callback(
    [
        Output("query-header-store", "data"),
        Output("query-seq-store", "data"),
        Output("query-type-store", "data"),
        Output("multifasta-sequences-store", "data"),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("query-text", "value"),
        Input("blast-fasta-upload", "contents"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
)
@handle_callback_error
def preprocess(n_clicks, query_text_input, query_file_contents):
    if not n_clicks:
        raise PreventUpdate

    try:
        # Check if we have a multifasta file
        multifasta_data = None

        if query_file_contents:
            logger.info("Processing file contents in preprocess")
            sequences, seq_count, error = parse_multifasta_from_file(
                query_file_contents
            )
            logger.info(
                f"Multifasta parsing result: {seq_count} sequences, error: {error}"
            )

            if seq_count > 1:
                # We have a valid multifasta file
                logger.info(f"Processing multifasta with {seq_count} sequences")
                multifasta_data = []

                # Store all sequences
                for idx, (header, seq) in enumerate(sequences):
                    seq_type = guess_seq_type(seq)
                    multifasta_data.append(
                        {
                            "index": idx,
                            "header": header,
                            "sequence": seq,
                            "type": seq_type,
                            "processed": False,
                        }
                    )
                logger.info(
                    f"Created multifasta_data list with {len(multifasta_data)} entries"
                )

                # For initial processing, return the first sequence
                first_seq = multifasta_data[0]
                logger.info(
                    f"Using first sequence for initial processing: {first_seq['header'][:30]}..."
                )
                return (
                    first_seq["header"],
                    first_seq["sequence"],
                    first_seq["type"],
                    multifasta_data,
                )
            else:
                # Only one sequence, process as single sequence
                logger.info("Only one sequence detected, processing as single sequence")
                input_type, query_header, query_seq = check_input(
                    query_text_input, query_file_contents
                )

                if input_type in ("none", "both"):
                    logger.info("Invalid input: none or both")
                    return None, None, None, None

                query_type = guess_seq_type(query_seq)
                logger.info(
                    f"Single sequence processed: {query_header[:30]}... of type {query_type}"
                )
                return query_header, query_seq, query_type, None
        else:
            # No file contents, process text input
            logger.info("No file contents, checking text input")
            input_type, query_header, query_seq = check_input(
                query_text_input, query_file_contents
            )

            if input_type in ("none", "both"):
                logger.info("Invalid input: none or both")
                return None, None, None, None

            query_type = guess_seq_type(query_seq)
            logger.info(
                f"Single sequence processed: {query_header[:30]}... of type {query_type}"
            )
            return query_header, query_seq, query_type, None

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        return None, None, None, None


@callback(
    Output("submission-id-store", "data"),
    Input("submit-button", "n_clicks"),
    prevent_initial_call=True,
)
def update_submission_id(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return n_clicks  # Use n_clicks as a unique submission ID


@callback(
    Output("subject-seq-dl-package", "data"),
    [Input("subject-seq-button", "n_clicks"), Input("subject-seq", "data")],
)
@handle_callback_error
def subject_seq_download(n_clicks, filename):
    try:
        if n_clicks:
            logger.info(f"Download initiated for file: {filename}")
            return dcc.send_file(filename)
        else:
            return dash.no_update
    except Exception as e:
        logger.error(f"Error in subject_seq_download: {str(e)}")
        return dash.no_update


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


@callback(
    [
        Output("blast-results-store", "data"),
        Output("captain-results-store", "data"),
        Output("subject-seq-button", "children"),
        Output("classification-interval", "disabled"),
        Output("classification-progress-section", "style", allow_duplicate=True),
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
        Input("submission-id-store", "data"),
        Input("multifasta-active-tab", "data"),
    ],
    [
        State("evalue-threshold", "value"),
        State("multifasta-sequences-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def fetch_captain(
    query_header,
    query_seq,
    query_type,
    submission_id,
    active_tab,
    evalue_threshold,
    multifasta_data,
    search_type="hmmsearch",
):
    # First, check if we have a valid submission
    if not all([query_header, query_seq, query_type, submission_id]):
        return None, None, None, True, {"display": "none"}

    # Check if we're dealing with multifasta data and need to switch to a different sequence
    if (
        multifasta_data
        and len(multifasta_data) > 1
        and active_tab is not None
        and active_tab > 0
    ):
        # Get the sequence for the active tab
        if active_tab < len(multifasta_data):
            tab_data = multifasta_data[active_tab]

            # If this tab hasn't been processed yet, use its sequence data
            if not tab_data.get("processed", False):
                logger.info(
                    f"Processing sequence for tab {active_tab}: {tab_data['header']}"
                )
                query_header = tab_data["header"]
                query_seq = tab_data["sequence"]
                query_type = tab_data["type"]

                # Mark this tab as processed in our local copy
                multifasta_data[active_tab]["processed"] = True

    try:
        captain_results_dict = None

        # Check sequence length - skip detailed classification for short sequences
        skip_classification = False
        if query_seq and len(query_seq) < 5000:  # 5kb threshold
            logger.info(
                f"Sequence length {len(query_seq)} is below 5kb threshold, skipping detailed classification"
            )
            skip_classification = True

        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)

        # Get the appropriate database configuration
        db = get_blast_db(query_type)

        if query_type == "nucl":
            # For nucleotide sequences, run diamond blastx first
            diamond_results = run_diamond(
                db_list=BLAST_DB_PATHS,  # Pass full paths dict
                query_type=query_type,
                input_gene="tyr",
                input_eval=evalue_threshold,
                query_fasta=tmp_query_fasta,
                threads=2,
            )

            if diamond_results:
                # Extract protein sequence from diamond results
                protein_seq = diamond_results[0].get("qseq_translated")
                if protein_seq:
                    # Write protein sequence to temp file
                    tmp_protein_fasta = write_temp_fasta(query_header, protein_seq)

                    # Run hmmsearch with protein sequence
                    captain_results_dict = run_hmmer(
                        db_list=BLAST_DB_PATHS,  # Pass full paths dict
                        query_type="prot",
                        input_gene="tyr",
                        input_eval=evalue_threshold,
                        query_fasta=tmp_protein_fasta,
                        threads=2,
                    )
        else:
            # For protein sequences, run diamond blastp
            diamond_results = run_diamond(
                db_list=BLAST_DB_PATHS,
                query_type="prot",
                input_gene="tyr",
                input_eval=evalue_threshold,
                query_fasta=tmp_query_fasta,
                threads=2,
            )

            if diamond_results:
                captain_results_dict = run_hmmer(
                    db_list=BLAST_DB_PATHS,
                    query_type="prot",
                    input_gene="tyr",
                    input_eval=evalue_threshold,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )

        # Run BLAST search for visualization
        blast_results_file = run_blast(
            db_list=db,  # Pass string path
            query_type=query_type,
            query_fasta=tmp_query_fasta,
            tmp_blast=tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name,
            input_eval=evalue_threshold,
            threads=2,
        )

        if blast_results_file is None:
            return (
                None,
                None,
                create_error_alert("No BLAST results were returned"),
                True,
                {"display": "none"},
            )

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
                {
                    "display": "none"
                },  # Hide progress since we're not doing classification
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

        def run_classification_workflow():
            try:
                # Process each stage in sequence
                for i, stage in enumerate(WORKFLOW_STAGES):
                    stage_id = stage["id"]

                    # Update status to show progress
                    globals()["workflow_tracker"].current_stage = stage_id
                    globals()["workflow_tracker"].current_stage_idx = i
                    globals()["workflow_tracker"].stage_progress = 0
                    globals()["workflow_tracker"].stage_values[i] = (
                        10  # Starting progress
                    )
                    globals()["workflow_tracker"].stage_animated[i] = True

                    logger.info(
                        f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}"
                    )

                    # Define task functions based on stage
                    result = None

                    if stage_id == "exact":
                        result = check_exact_matches_task.delay(
                            fasta=upload_data["fasta"],
                            ships_dict=fetch_ships(
                                **upload_data["fetch_ship_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Exact match result: {result}")

                        # If we found an exact match, exit early with a completed state
                        if result:
                            logger.info(
                                f"Exact match found: {result}, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = (
                                100  # Mark this stage as completed
                            )
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].found_match = True
                            globals()["workflow_tracker"].match_stage = "exact"
                            globals()["workflow_tracker"].match_result = result
                            globals()["workflow_tracker"].complete = True
                            break

                    elif stage_id == "contained":
                        globals()["workflow_tracker"].stage_progress = 30
                        globals()["workflow_tracker"].stage_values[i] = 30
                        result = check_contained_matches_task.delay(
                            fasta=upload_data["fasta"],
                            ships_dict=fetch_ships(
                                **upload_data["fetch_ship_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Contained match result: {result}")

                        # If we found a contained match, exit early
                        if result:
                            logger.info(
                                f"Contained match found: {result}, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = 100
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].found_match = True
                            globals()["workflow_tracker"].match_stage = "contained"
                            globals()["workflow_tracker"].match_result = result
                            globals()["workflow_tracker"].complete = True
                            break

                    elif stage_id == "similar":
                        globals()["workflow_tracker"].stage_progress = 30
                        globals()["workflow_tracker"].stage_values[i] = 30
                        result = check_similar_matches_task.delay(
                            fasta=upload_data["fasta"],
                            ships_dict=fetch_ships(
                                **upload_data["fetch_ship_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Similar match result: {result}")

                        # If we found a similar match, exit early
                        if result:
                            logger.info(
                                f"Similar match found: {result}, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = 100
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].found_match = True
                            globals()["workflow_tracker"].match_stage = "similar"
                            globals()["workflow_tracker"].match_result = result
                            globals()["workflow_tracker"].complete = True
                            break

                    elif stage_id == "family":
                        globals()["workflow_tracker"].stage_progress = 30
                        globals()["workflow_tracker"].stage_values[i] = 30
                        result = run_family_classification_task.delay(
                            fasta=upload_data["fasta"],
                            seq_type=upload_data["seq_type"],
                            db_list=fetch_ships(
                                **upload_data["fetch_ship_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Family classification result: {result}")

                        # For family stage, if result is None, stop the pipeline with appropriate message
                        if result is None:
                            logger.info(
                                "Family classification did not find a match, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = (
                                100  # Mark this stage as completed
                            )
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].complete = True
                            globals()[
                                "workflow_tracker"
                            ].match_stage = (
                                "family"  # Not really a match but for display purposes
                            )
                            globals()[
                                "workflow_tracker"
                            ].match_result = "No family match found"
                            break

                        # Store family result for use in subsequent stages
                        if isinstance(result, dict) and "protein" in result:
                            upload_data["protein_file"] = result["protein"]

                    elif stage_id == "navis":
                        globals()["workflow_tracker"].stage_progress = 30
                        globals()["workflow_tracker"].stage_values[i] = 30

                        # Check if we have a protein file from the family stage
                        if (
                            "protein_file" not in upload_data
                            or not upload_data["protein_file"]
                        ):
                            logger.warning(
                                "No protein file available for navis classification"
                            )
                            globals()["workflow_tracker"].stage_values[i] = 100
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].complete = True
                            globals()["workflow_tracker"].match_stage = "navis"
                            globals()[
                                "workflow_tracker"
                            ].match_result = "No protein data available for navis"
                            break

                        result = run_navis_classification_task.delay(
                            fasta=upload_data["protein_file"],
                            existing_ships=fetch_captains(
                                **upload_data["fetch_captain_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Navis classification result: {result}")

                        # For navis stage, if result is None, stop the pipeline with appropriate message
                        if result is None:
                            logger.info(
                                "Navis classification did not find a match, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = 100
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].complete = True
                            globals()["workflow_tracker"].match_stage = "navis"
                            globals()[
                                "workflow_tracker"
                            ].match_result = "No navis match found"
                            break

                    elif stage_id == "haplotype":
                        globals()["workflow_tracker"].stage_progress = 30
                        globals()["workflow_tracker"].stage_values[i] = 30
                        result = run_haplotype_classification_task.delay(
                            fasta=upload_data["fasta"],
                            existing_ships=fetch_ships(
                                **upload_data["fetch_ship_params"]
                            ).to_dict("records"),
                            navis=fetch_captains(
                                **upload_data["fetch_captain_params"]
                            ).to_dict("records"),
                        ).get(timeout=180)
                        logger.debug(f"Haplotype classification result: {result}")

                        # For haplotype stage, if result is None, stop the pipeline with appropriate message
                        if result is None:
                            logger.info(
                                "Haplotype classification did not find a match, stopping pipeline"
                            )
                            globals()["workflow_tracker"].stage_values[i] = 100
                            globals()["workflow_tracker"].stage_animated[i] = False
                            globals()["workflow_tracker"].complete = True
                            globals()["workflow_tracker"].match_stage = "haplotype"
                            globals()[
                                "workflow_tracker"
                            ].match_result = "No haplotype match found"
                            break

                    else:
                        logger.warning(
                            f"No task defined for stage {stage_id}, skipping"
                        )
                        globals()["workflow_tracker"].stage_values[i] = 100
                        globals()["workflow_tracker"].stage_animated[i] = False
                        continue

                    # Update the status to show completion for this stage
                    globals()["workflow_tracker"].stage_progress = 100
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False

                # Mark as complete if we haven't already done so
                globals()["workflow_tracker"].complete = True

            except Exception as e:
                # Handle errors - only unexpected errors
                error_message = f"Error during classification: {str(e)}"
                logger.error(error_message)
                logger.exception("Full traceback:")  # Log full traceback for debugging
                globals()["workflow_tracker"].error = error_message
                globals()["workflow_tracker"].complete = True

                # Mark current stage as failed
                if globals()["workflow_tracker"].current_stage_idx is not None:
                    i = globals()["workflow_tracker"].current_stage_idx
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].stage_striped[i] = False

        # Start background thread to run classification
        classification_thread = Thread(target=run_classification_workflow)
        classification_thread.daemon = True
        classification_thread.start()

        return (
            blast_results_file,
            captain_results_dict,
            None,
            False,
            {"display": "block"},
        )

    except Exception as e:
        logger.error(f"Error in fetch_captain: {str(e)}")
        return None, None, create_error_alert(str(e)), True, {"display": "none"}


# 3. UI Update Callback - Modified to use classification results
@callback(
    [
        Output("classification-output", "children", allow_duplicate=True),
        Output("blast-download", "children"),
    ],
    [
        Input("blast-results-store", "data"),
        Input("captain-results-store", "data"),
        Input("classification-workflow-state", "data"),
        Input("classification-result-store", "data"),
        Input("multifasta-sequences-store", "data"),
    ],
    [
        State("submit-button", "n_clicks"),
        State("evalue-threshold", "value"),
        State("multifasta-active-tab", "data"),
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
    captain_results_dict,
    classification_state,
    classification_result,
    multifasta_data,
    n_clicks,
    evalue,
    active_tab,
):
    if not n_clicks or blast_results_file is None:
        return None, None

    # Check if we're working with multifasta
    is_multifasta = multifasta_data and len(multifasta_data) > 1

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

        # If we're working with multifasta, we may need to handle tabs
        if is_multifasta:
            # Check if we've already created tabs
            ctx = dash.callback_context
            triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

            # Return the outputs to be rendered in the appropriate tab
            if triggered_id == "multifasta-sequences-store":
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
        State("subject-seq-button", "children"),
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


# First callback: Handle file upload
@callback(
    [
        Output("classification-upload", "data"),
        Output("classification-file-info", "children"),
        Output("classification-file-info", "style"),
        Output("classification-submit-button", "style"),
    ],
    [
        Input("classification-fasta-upload", "contents"),
        Input("classification-fasta-sequence-upload", "contents"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_file_upload(seq_content, upload_contents):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    logger.info(f"File upload triggered by: {trigger_id}")

    # Use whichever content is provided
    content_to_use = seq_content if seq_content is not None else upload_contents

    if content_to_use is None:
        raise PreventUpdate

    # Parse and validate FASTA
    header, seq, fasta_error = parse_fasta_from_file(content_to_use)
    if fasta_error or not seq:
        raise ValueError(f"FASTA parsing error: {fasta_error or 'No sequence found'}")

    seq_type = guess_seq_type(seq)
    if not seq_type:
        raise ValueError("Could not determine sequence type")

    # Save sequence to temporary file
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_fasta({"query_sequence": seq}, tmp_fasta)

    # Prepare upload data
    upload_data = {
        "seq_type": seq_type,
        "fasta": tmp_fasta,
        "fetch_ship_params": {
            "curated": False,
            "with_sequence": True,
            "dereplicate": True,
        },
        "fetch_captain_params": {"curated": True, "with_sequence": True},
    }

    # Create file info display
    file_info = dmc.Alert(
        title="Sequence uploaded successfully",
        children=[
            html.P(f"Sequence header: {header}"),
            html.P(f"Sequence type: {seq_type}"),
            html.P(f"Sequence length: {len(seq)} bp"),
        ],
        color="green",
        variant="light",
    )

    return [
        upload_data,
        file_info,
        {"display": "block"},  # Show file info
        {"display": "block"},  # Show submit button
    ]


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


# Second callback: Initialize workflow and start the background process
@callback(
    [
        Output("classification-workflow-state", "data", allow_duplicate=True),
        Output("classification-progress-section", "style", allow_duplicate=True),
    ],
    Input("classification-submit-button", "n_clicks"),
    [
        State("classification-upload", "data"),
    ],
    prevent_initial_call=True,
)
def run_workflow_background(n_clicks, upload_data):
    if n_clicks is None or upload_data is None:
        raise PreventUpdate

    global workflow_tracker
    workflow_tracker = WorkflowStatus()
    workflow_tracker.start_time = time.time()

    logger.info("Starting classification workflow in background")

    try:
        # Process each stage in sequence
        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]

            # Update status to show progress
            workflow_tracker.current_stage = stage_id
            workflow_tracker.current_stage_idx = i
            workflow_tracker.stage_progress = 0
            workflow_tracker.stage_values[i] = 10  # Starting progress
            workflow_tracker.stage_animated[i] = True

            logger.info(f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}")

            # Define task functions based on stage
            result = None

            if stage_id == "exact":
                result = check_exact_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict(
                        "records"
                    ),
                ).get(timeout=180)
                logger.debug(f"Exact match result: {result}")

                # If we found an exact match, exit early with a completed state
                if result:
                    logger.info(f"Exact match found: {result}, stopping pipeline")
                    globals()["workflow_tracker"].stage_values[i] = (
                        100  # Mark this stage as completed
                    )
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].found_match = True
                    globals()["workflow_tracker"].match_stage = "exact"
                    globals()["workflow_tracker"].match_result = result
                    globals()["workflow_tracker"].complete = True
                    break

            elif stage_id == "contained":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = check_contained_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict(
                        "records"
                    ),
                ).get(timeout=180)
                logger.debug(f"Contained match result: {result}")

                # If we found a contained match, exit early
                if result:
                    logger.info(f"Contained match found: {result}, stopping pipeline")
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].found_match = True
                    globals()["workflow_tracker"].match_stage = "contained"
                    globals()["workflow_tracker"].match_result = result
                    globals()["workflow_tracker"].complete = True
                    break

            elif stage_id == "similar":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = check_similar_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict(
                        "records"
                    ),
                ).get(timeout=180)
                logger.debug(f"Similar match result: {result}")

                # If we found a similar match, exit early
                if result:
                    logger.info(f"Similar match found: {result}, stopping pipeline")
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].found_match = True
                    globals()["workflow_tracker"].match_stage = "similar"
                    globals()["workflow_tracker"].match_result = result
                    globals()["workflow_tracker"].complete = True
                    break

            elif stage_id == "family":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = run_family_classification_task.delay(
                    fasta=upload_data["fasta"],
                    seq_type=upload_data["seq_type"],
                    db_list=fetch_ships(**upload_data["fetch_ship_params"]).to_dict(
                        "records"
                    ),
                ).get(timeout=180)
                logger.debug(f"Family classification result: {result}")

                # For family stage, if result is None, stop the pipeline with appropriate message
                if result is None:
                    logger.info(
                        "Family classification did not find a match, stopping pipeline"
                    )
                    globals()["workflow_tracker"].stage_values[i] = (
                        100  # Mark this stage as completed
                    )
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].complete = True
                    globals()[
                        "workflow_tracker"
                    ].match_stage = (
                        "family"  # Not really a match but for display purposes
                    )
                    globals()["workflow_tracker"].match_result = "No family match found"
                    break

                # Store family result for use in subsequent stages
                if isinstance(result, dict) and "protein" in result:
                    upload_data["protein_file"] = result["protein"]

            elif stage_id == "navis":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30

                # Check if we have a protein file from the family stage
                if "protein_file" not in upload_data or not upload_data["protein_file"]:
                    logger.warning("No protein file available for navis classification")
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].complete = True
                    globals()["workflow_tracker"].match_stage = "navis"
                    globals()[
                        "workflow_tracker"
                    ].match_result = "No protein data available for navis"
                    break

            elif stage_id == "haplotype":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = run_haplotype_classification_task.delay(
                    fasta=upload_data["fasta"],
                    existing_ships=fetch_ships(
                        **upload_data["fetch_ship_params"]
                    ).to_dict("records"),
                    navis=fetch_captains(**upload_data["fetch_captain_params"]).to_dict(
                        "records"
                    ),
                ).get(timeout=180)
                logger.debug(f"Haplotype classification result: {result}")

                # For haplotype stage, if result is None, stop the pipeline with appropriate message
                if result is None:
                    logger.info(
                        "Haplotype classification did not find a match, stopping pipeline"
                    )
                    globals()["workflow_tracker"].stage_values[i] = 100
                    globals()["workflow_tracker"].stage_animated[i] = False
                    globals()["workflow_tracker"].complete = True
                    globals()["workflow_tracker"].match_stage = "haplotype"
                    globals()[
                        "workflow_tracker"
                    ].match_result = "No haplotype match found"
                    break

            else:
                logger.warning(f"No task defined for stage {stage_id}, skipping")
                globals()["workflow_tracker"].stage_values[i] = 100
                globals()["workflow_tracker"].stage_animated[i] = False
                continue

            # Update the status to show completion for this stage
            globals()["workflow_tracker"].stage_progress = 100
            globals()["workflow_tracker"].stage_values[i] = 100
            globals()["workflow_tracker"].stage_animated[i] = False

        # Mark as complete if we haven't already done so
        globals()["workflow_tracker"].complete = True

    except Exception as e:
        # Handle errors - only unexpected errors
        error_message = f"Error during classification: {str(e)}"
        logger.error(error_message)
        logger.exception("Full traceback:")  # Log full traceback for debugging
        globals()["workflow_tracker"].error = error_message
        globals()["workflow_tracker"].complete = True

        # Mark current stage as failed
        if globals()["workflow_tracker"].current_stage_idx is not None:
            i = globals()["workflow_tracker"].current_stage_idx
            globals()["workflow_tracker"].stage_values[i] = 100
            globals()["workflow_tracker"].stage_animated[i] = False
            globals()["workflow_tracker"].stage_striped[i] = False

    # Return the final state - no longer controlling the interval here
    return workflow_tracker.to_dict(), {"display": "block"}


# Third callback to update the UI periodically based on classification status
@callback(
    [
        Output("classification-stage", "data", allow_duplicate=True),
        Output("classification-stage-display", "children", allow_duplicate=True),
        Output("classification-progress", "value", allow_duplicate=True),
        Output("classification-progress", "animated", allow_duplicate=True),
        Output("classification-progress", "striped", allow_duplicate=True),
        Output("classification-result-store", "data", allow_duplicate=True),
        Output("classification-interval", "disabled", allow_duplicate=True),
    ],
    Input("classification-interval", "n_intervals"),
    prevent_initial_call=True,
)
@handle_callback_error
def update_ui_from_status(n_intervals):
    """Update UI based on the current workflow status."""
    global workflow_tracker

    if workflow_tracker is None:
        raise PreventUpdate

    # Get current status
    status = workflow_tracker

    # Default values
    stage_message = "Starting classification..."
    result_data = None

    # Determine if we should disable the interval after this update
    disable_interval = False

    if status.complete:
        # Disable interval after this update since workflow is complete
        disable_interval = True

        if status.error:
            # True error condition
            stage_message = "Error occurred during classification"
            result_data = {"type": "error", "message": status.error}
        elif status.found_match:
            # Successful match
            stage_info = next(
                (s for s in WORKFLOW_STAGES if s["id"] == status.match_stage), {}
            )
            stage_message = (
                f"{stage_info.get('label', 'Match')} found: {status.match_result}"
            )

            # Create result data for the store
            result_data = {
                "type": "match",
                "match_stage": status.match_stage,
                "match_result": status.match_result,
                "stage_color": stage_info.get("color", "blue"),
                "stage_label": stage_info.get("label", "Match"),
            }
        elif status.match_stage:
            # No match found, but we know which stage stopped the pipeline
            stage_info = next(
                (s for s in WORKFLOW_STAGES if s["id"] == status.match_stage), {}
            )
            no_match_message = f"No {status.match_stage} match found"

            if status.match_result and "No" in status.match_result:
                no_match_message = status.match_result

            stage_message = f"Classification stopped: {no_match_message}"

            result_data = {
                "type": "no_match",
                "match_stage": status.match_stage,
                "message": no_match_message,
                "stage_color": stage_info.get("color", "yellow"),
            }
        else:
            # Generic no match case (all stages completed)
            stage_message = "Classification complete - No matches found"
            result_data = {
                "type": "no_match",
                "message": "No matches were found for this sequence",
            }
    elif status.current_stage:
        # Still processing a stage
        stage_info = next(
            (s for s in WORKFLOW_STAGES if s["id"] == status.current_stage), {}
        )
        stage_message = f"Processing: {stage_info.get('label', status.current_stage)}"

    # Calculate overall progress
    if status.complete:
        overall_progress = 100
        overall_animated = False
        overall_striped = False
    elif status.current_stage_idx is not None:
        # Calculate progress based on current stage
        overall_progress = ((status.current_stage_idx) / len(WORKFLOW_STAGES)) * 100
        if status.stage_progress:
            # Add partial progress from current stage
            stage_contribution = (1 / len(WORKFLOW_STAGES)) * (
                status.stage_progress / 100
            )
            overall_progress += stage_contribution * 100
        overall_animated = True
        overall_striped = True
    else:
        overall_progress = 0
        overall_animated = True
        overall_striped = True

    # Log the final result for debugging
    if status.complete and result_data is not None:
        logger.info(f"Classification complete: {stage_message}")

    return (
        stage_message,  # classification-stage data
        stage_message,  # classification-stage-display
        overall_progress,  # overall progress value
        overall_animated,  # overall animated state
        overall_striped,  # overall striped state
        result_data,  # Data for the result store
        disable_interval,  # Whether to disable the interval
    )


# Add a separate callback just to enable the interval when workflow starts
@callback(
    Output("classification-workflow-interval", "disabled"),
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
        Input("multifasta-sequences-store", "data"),
        Input("blast-results-store", "data"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def update_tabs_for_multifasta(multifasta_data, blast_results):
    # Add debug logging
    logger.info(
        f"update_tabs_for_multifasta called with multifasta_data: {type(multifasta_data)}, blast_results: {type(blast_results)}"
    )

    # If there's no multifasta data or blast results, don't update
    if not multifasta_data:
        logger.info("No multifasta_data, preventing update")
        raise PreventUpdate

    if not blast_results:
        logger.info("No blast_results, preventing update")
        raise PreventUpdate

    # If we have multifasta data (more than one sequence)
    if isinstance(multifasta_data, list) and len(multifasta_data) > 1:
        logger.info(f"Creating tabs for {len(multifasta_data)} sequences")

        # Create tabs for each sequence
        tabs = []
        for idx, seq in enumerate(multifasta_data[:10]):  # Limit to 10 sequences
            logger.info(f"Creating tab {idx}: {seq['header'][:20]}")
            tabs.append(
                dbc.Tab(
                    label=f"Sequence {idx + 1}: {seq['header'][:20]}{'...' if len(seq['header']) > 20 else ''}",
                    tab_id=f"tab-{idx}",
                )
            )

        logger.info(f"Created {len(tabs)} tab objects")

        # Create a simple Card with Tabs in the header
        card = dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Tabs(id="multifasta-tabs", active_tab="tab-0", children=tabs)
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
        return card

    # If there's only one sequence, return the default layout
    logger.info("No multifasta data or only one sequence, returning default layout")
    return dmc.Stack(
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


@callback(
    [
        Output("multifasta-active-tab", "data"),
        Output("multifasta-processed-tabs", "data"),
        Output("tab-content", "children"),
    ],
    [Input("multifasta-tabs", "active_tab")],
    [
        State("multifasta-processed-tabs", "data"),
        State("multifasta-sequences-store", "data"),
        State("blast-results-store", "data"),
        State("captain-results-store", "data"),
        State("evalue-threshold", "value"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_tab_switch(
    active_tab,
    processed_tabs,
    multifasta_data,
    blast_results,
    captain_results,
    evalue_threshold,
):
    """Track which tab is active and update processed tabs list"""
    logger.info(f"handle_tab_switch called with active_tab: {active_tab}")

    if active_tab is None or not multifasta_data:
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
    if tab_idx >= len(multifasta_data):
        logger.error(
            f"Tab index {tab_idx} is out of range for multifasta_data length {len(multifasta_data)}"
        )
        raise PreventUpdate

    seq_data = multifasta_data[tab_idx]
    header = seq_data["header"]
    sequence = seq_data["sequence"]
    seq_type = seq_data["type"]
    logger.info(f"Processing sequence for tab {tab_idx}: {header[:30]}...")

    # Process this sequence for the tab
    # First write sequence to temporary FASTA file
    tmp_query_fasta = write_temp_fasta(header, sequence)
    logger.info(f"Wrote sequence to temporary file: {tmp_query_fasta}")

    # Get the appropriate database configuration
    db = get_blast_db(seq_type)
    logger.info(f"Using database for {seq_type} sequence type")

    # Run BLAST search for this sequence
    logger.info(f"Running BLAST search for tab {tab_idx}")
    blast_results_file = run_blast(
        db_list=db,
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
                html.Div(id="subject-seq-button"),
                # Add the blast data store to trigger BlasterJS
                blast_data_store if blast_results_file is not None else None,
            ]
        ),
    ]

    logger.info(
        f"Returning data for tab {tab_idx} with content length {len(tab_content)}"
    )
    return tab_idx, processed_tabs, tab_content


@callback(
    Output("right-column-content-debug", "children"),
    Input("multifasta-sequences-store", "data"),
    prevent_initial_call=True,
)
def debug_multifasta_data(multifasta_data):
    """Debug callback to check if multifasta data is being correctly stored"""
    if not multifasta_data:
        logger.info("No multifasta data")

    if isinstance(multifasta_data, list):
        logger.info(f"Multifasta data contains {len(multifasta_data)} sequences")
    else:
        logger.info(f"Multifasta data is not a list: {type(multifasta_data)}")

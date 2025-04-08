import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import dcc, html, callback, clientside_callback, ClientsideFunction

from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import tempfile
import base64
import pandas as pd
import logging
import subprocess
import os


from src.config.cache import cache
from src.utils.seq_utils import ( guess_seq_type,
    check_input,
    write_temp_fasta,
    parse_fasta,
    parse_fasta_from_file
    )
from src.utils.blast_utils import (
    run_blast,
    run_hmmer,
    run_diamond,
    make_captain_alert,
    process_captain_results,
    create_no_matches_alert,
    parse_blast_xml,
    blast_download_button,
)

from src.components.callbacks import curated_switch, create_file_upload
from src.database.sql_manager import fetch_meta_data
from src.config.settings import BLAST_DB_PATHS
from src.utils.telemetry import get_client_ip, get_blast_limit_info, blast_limit_decorator
from src.components.error_boundary import handle_callback_error, create_error_alert
from src.utils.tree import plot_tree

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
            dcc.Store(id="query-header-store"),
            dcc.Store(id="query-seq-store"),
            dcc.Store(id="query-type-store"),
            dcc.Store(id="blast-results-store"),
            dcc.Store(id="captain-results-store"),
            dcc.Store(id="upload-error-store"),
            dcc.Store(id="processed-metadata-store"),
            dcc.Store(id="processed-blast-store"),
            dcc.Store(id="submission-id-store"),            
            dcc.Store(id="timeout-store", data=False),
            dcc.Interval(id="timeout-interval", interval=30000, n_intervals=0),  # 30 seconds
        dmc.Space(h=20),       
        dmc.Grid(
            children=[
                dmc.GridCol(
                    span={"sm": 12, "lg": 4},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                dmc.Title("BLAST Search", order=1),
                                dmc.Text(
                                    "Search protein/nucleotide sequences for Starships and Starship-associated genes",
                                    c="dimmed",
                                    size="lg",
                                ),
                                # Input Section
                                dmc.Stack([
                                    dmc.Title("Input Sequence", order=3),
                                    dmc.Textarea(
                                        id="query-text",
                                        placeholder="Paste FASTA sequence here...",
                                        minRows=5,
                                        style={"width": "100%"},
                                    ),
                                ], gap="xs"),
                                
                                # Upload Section
                                dmc.Stack([
                                    dmc.Center(
                                        dmc.Text("Or", size="lg"),
                                    ),
                                    dmc.Paper(
                                        children=create_file_upload(
                                            upload_id="blast-fasta-upload",
                                            output_id="blast-fasta-sequence-upload",
                                            accept_types=[".fa", ".fas", ".fasta", ".fna"],
                                            placeholder_text="Drag and drop or click to select a FASTA file"
                                        ),
                                        withBorder=False,
                                        shadow="sm",
                                        radius="md",
                                        style={"cursor": "pointer"}
                                    ),
                                    html.Div(
                                        id="upload-error-message",
                                        style={"color": "red"}
                                    ),
                                ], gap="md"),
                                
                                # Options Section
                                dmc.Stack([
                                    dmc.Title("Search Options", order=3),
                                    curated_switch(
                                        text="Only search curated Starships",
                                        size="sm"
                                    ),
                                ], gap="xs"),
                                
                                # Submit Section
                                dmc.Stack([
                                    dmc.Center(
                                        dmc.Button(
                                            "Submit BLAST",
                                            id="submit-button",
                                            variant="gradient",
                                            gradient={"from": "indigo", "to": "cyan"},
                                            size="lg",
                                            loading=False,
                                            loaderProps={"variant": "dots", "color": "white"},
                                        ),
                                    ),
                                ], gap="xs"),                                
                            ], gap="md"),
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
                            children=[
                                html.Div(id="ship-family"),
                                html.Div(id="phylogeny-plot"),
                                html.Div([
                                    html.Div(id='blast-multiple-alignments'),
                                    html.Div(id='blast-alignments-table'),
                                ], id='blast-container'),
                                html.Div(id="blast-download"),
                                html.Div(id="subject-seq-button"),
                            ],
                        ),
                    ],
                    style={"justify-content": "flex-start"}
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
                // Clear existing content
                document.getElementById('blast-multiple-alignments').innerHTML = '';
                document.getElementById('blast-alignments-table').innerHTML = '';
                
                // Initialize BlasterJS
                var blasterjs = require("biojs-vis-blasterjs");
                var instance = new blasterjs({
                    string: data.blast_text,
                    multipleAlignments: "blast-multiple-alignments",
                    alignmentsTable: "blast-alignments-table"
                });
                
                // Add CSS to control table width
                var style = document.createElement('style');
                style.textContent = `
                    #blast-alignments-table {
                        max-width: 100%;
                        overflow-x: auto;
                    }
                    #blast-alignments-table table {
                        width: 100%;
                        table-layout: fixed;
                    }
                    #blast-alignments-table td {
                        max-width: 200px;
                        overflow: hidden;
                        text-overflow: ellipsis;
                        white-space: nowrap;
                    }
                `;
                document.head.appendChild(style);
                
                return window.dash_clientside.no_update;
            } catch (error) {
                console.error('Error initializing BlasterJS:', error);
                return error.toString();
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output('blast-container', 'children'),
    Input('processed-blast-store', 'data'),
    prevent_initial_call=True
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
        Input("timeout-interval", "n_intervals")
    ],
    [
        State("blast-fasta-upload", "filename"),
        State("timeout-store", "data")
    ],
    prevent_initial_call=True
)
@handle_callback_error
def handle_submission_and_upload(n_clicks, contents, n_intervals, filename, timeout_triggered):
    triggered_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    
    # Default return values
    button_disabled = False
    button_text = "Submit BLAST"
    error_message = ""
    error_store = None
    
    # Handle timeout case
    if triggered_id == "timeout-interval":
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
    
    return [
        button_disabled,
        button_text,
        error_message,
        error_store
    ]

@blast_limit_decorator
@callback(
    [
        Output("query-header-store", "data"),
        Output("query-seq-store", "data"),
        Output("query-type-store", "data"),
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
        input_type, query_header, query_seq = check_input(
            query_text_input, query_file_contents
        )
        
        if input_type in ("none", "both"):
            return None, None, None

        query_type = guess_seq_type(query_seq)
        return query_header, query_seq, query_type

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        return None, None, None


@callback(
    Output("submission-id-store", "data"),
    Input("submit-button", "n_clicks"),
    prevent_initial_call=True
)
def update_submission_id(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return n_clicks  # Use n_clicks as a unique submission ID


@callback(
    [
        Output("blast-results-store", "data"),
        Output("captain-results-store", "data"),
        Output("subject-seq-button", "children"),
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
        Input("submission-id-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def fetch_captain(query_header, query_seq, query_type, submission_id, search_type="hmmsearch"):
    if not all([query_header, query_seq, query_type, submission_id]):
        return None, None, None
        
    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        
        # Instead of using timeout context manager, we'll rely on subprocess timeout
        try:
            blast_results_file = run_blast(
                db_list=BLAST_DB_PATHS,
                query_type=query_type,
                query_fasta=tmp_query_fasta,
                tmp_blast=tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name,
                input_eval=0.01,
                threads=2,
            )
            

            if blast_results_file is None:
                error_div = html.Div([
                    dmc.Alert(
                        title="BLAST Error",
                        color="red",
                        children="No BLAST results were returned. Please try again with a different sequence.",
                    )
                ])
                return None, None, error_div
           
            if search_type == "diamond":
                captain_results_dict = run_diamond(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
            else:
                captain_results_dict = run_hmmer(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
                
            return blast_results_file, captain_results_dict, None

        except subprocess.TimeoutExpired:
            error_div = html.Div([
                dmc.Alert(
                    title="Operation Timeout",
                    color="red",
                    children="The BLAST search took too long to complete. Please try with a shorter sequence or try again later.",
                )
            ])
            return None, None, error_div
            
    except Exception as e:
        logger.error(f"Error in fetch_captain: {str(e)}")
        error_div = html.Div([
            dmc.Alert(
                title="Error",
                color="red",
                children=f"An error occurred: {str(e)}"
            )
        ])
        return None, None, error_div


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


# 1. Metadata Processing Callback
@cache.memoize()
@callback(
    Output("processed-metadata-store", "data"),
    Input("curated-input", "value"),
)
@handle_callback_error
def process_metadata(curated):
    try:
        initial_df = cache.get("meta_data")
        if initial_df is None:
            initial_df = fetch_meta_data(curated)
        
        return initial_df[["accession_tag", "familyName"]].drop_duplicates().to_dict('records')
    except Exception as e:
        logger.error(f"Error processing metadata: {str(e)}")
        return None

# 2. BLAST Results Processing Callback
@callback(
    Output("processed-blast-store", "data"),
    [Input("blast-results-store", "data")],
)
@handle_callback_error
def process_blast_results(blast_results_file):
    if not blast_results_file:
        return None
        
    try:
        # Read the BLAST output as text
        with open(blast_results_file, 'r') as f:
            blast_results = f.read()

        # Add size limit check
        results_size = len(blast_results)
        if results_size > 5 * 1024 * 1024:  # 5MB limit
            logger.warning(f"BLAST results too large: {results_size} bytes")
            return None

        # Format data for BlasterJS
        data = {
            'blast_text': blast_results  # Pass the raw BLAST text directly
        }
        
        return data
        
    except Exception as e:
        logger.error(f"Error processing BLAST results: {str(e)}")
        return None

# 3. UI Update Callback
@callback(
    [
        Output("ship-family", "children"),
        Output("blast-download", "children"),
        Output("phylogeny-plot", "children")
    ],
    [Input("blast-results-store", "data"), 
     Input("captain-results-store", "data")],
    State("submit-button", "n_clicks"),
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def update_ui_elements(blast_results_file, captain_results_dict, n_clicks):
    if not n_clicks or blast_results_file is None:
        return None, None, None
        
    try:
        if not blast_results_file:
            return html.Div(create_no_matches_alert()), None, None

        # Process blast results first for more basic family identification
        blast_tsv = parse_blast_xml(blast_results_file)
        blast_df = pd.read_csv(blast_tsv, sep="\t")
        
        # Check if dataframe is empty
        if blast_df.empty:
            return html.Div(create_no_matches_alert()), None, None
                    
        # Simple selection of top hit
        if len(blast_df) > 0:
            # Sort by evalue (ascending) and pident (descending) to get best hits
            blast_df = blast_df.sort_values(['evalue', 'pident'], ascending=[True, False])
            top_hit = blast_df.iloc[0]
            logger.info(f"Top hit: {top_hit}")

            top_evalue = float(top_hit['evalue'])
            top_aln_length = int(top_hit['aln_length'])
            top_pident = float(top_hit['pident'])

            if top_pident >= 90:
                # look up family name from accession tag
                query_accession = top_hit['query_id']
                top_family = fetch_meta_data(accession_tag=query_accession)['familyName']
                ship_family = make_captain_alert(top_family, top_aln_length, top_evalue, "blast")
            else:
                # Process captain results
                ship_family = process_captain_results(captain_results_dict)
        else:
            return html.Div(create_no_matches_alert()), None, None

        # Create download button
        download_button = blast_download_button()
        
        # Create phylogeny plot if captain results are available
        phylogeny_plot = None
        if captain_results_dict:
            try:
                # Handle case where captain_results_dict is a list of results
                if isinstance(captain_results_dict, list):
                    # Get the first result if available
                    family_name = captain_results_dict[0].get('family_name') if captain_results_dict else None
                else:
                    # Original behavior for dict
                    family_name = captain_results_dict.get('family_name')
                
                if family_name:
                    fig = plot_tree(
                        highlight_families=family_name,
                        tips=None  # We don't have specific tips to highlight
                    )
                    
                    phylogeny_plot = html.Div([
                        dbc.Accordion(
                            children=[
                                dbc.AccordionItem(
                                    children=[
                                        dcc.Graph(
                                            id="phylogeny-graph",
                                            figure=fig,
                                            style={'height': '800px'}
                                        )
                                    ],
                                    title="Captain Phylogeny",
                                    item_id="phylogeny",
                                )
                            ],
                            always_open=False,
                            active_item=None,
                            id="phylogeny-accordion"
                        )
                    ])
            except Exception as e:
                logger.error(f"Error creating phylogeny plot: {e}")
                phylogeny_plot = html.Div(create_error_alert("Could not create phylogeny plot"))
        
        return ship_family, download_button, phylogeny_plot
        
    except Exception as e:
        logger.error(f"Error in update_ui_elements: {e}")
        error_div = html.Div(create_error_alert(str(e)))
        return error_div, None, None

@callback(
    Output("timeout-store", "data"),
    [Input("timeout-interval", "n_intervals")],
    [
        State("submit-button", "n_clicks"),
        State("blast-container", "children"),
        State("ship-family", "children"),
        State("subject-seq-button", "children"),
    ]
)
def check_timeout(n_intervals, n_clicks, blast_container, family_content, error_content):
    if n_clicks and n_intervals > 0:
        has_output = any([
            blast_container is not None,
            family_content is not None,
            error_content is not None
        ])
        return not has_output
    return False

clientside_callback(
    """
    function(n_intervals) {
        window.addEventListener('blastAccessionClick', function(event) {
            if (event.detail && event.detail.accession) {
                // Find the store and update it
                var store = document.getElementById('blast-modal-data');
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
    Input("blast-modal-data", "data"),
    prevent_initial_call=True
)
def update_modal_content(accession):
    from src.components.callbacks import create_accession_modal
    if not accession:
        raise PreventUpdate
    
    try:
        modal_content, _ = create_accession_modal(accession)
        return modal_content
    except Exception as e:
        logger.error(f"Error in update_modal_content: {str(e)}")
        return html.Div(
            f"Error loading details: {str(e)}", 
            style={"color": "red", "padding": "20px"}
        )
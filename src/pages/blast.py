from dash import Dash
import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import dcc, html, callback, MATCH, no_update
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State
import dash_bio as dashbio

import io
import re
import tempfile
import base64
from datetime import date
import pandas as pd
import plotly.graph_objects as go
import logging
from flask import jsonify, request
import subprocess
import os
import json


from src.config.cache import cache
from src.utils.seq_utils import ( guess_seq_type,
    check_input,
    write_temp_fasta,
                                 )
from src.utils.blast_utils import (
    run_blast,
    run_hmmer,
    run_diamond,
    blast_table,
    run_lastz,
    select_ship_family,
    parse_lastz_output,
    blast_chords,
)

from src.utils.blast_utils import (
    run_blast,
    run_hmmer,
    run_diamond,
    blast_table,
    run_lastz,
    select_ship_family,
    parse_lastz_output,
    blast_chords,
)
from src.components.callbacks import curated_switch, create_accession_modal, create_modal_callback
from src.utils.seq_utils import parse_fasta, parse_fasta_from_file

from src.database.sql_manager import fetch_meta_data
from src.database.blastdb import db_list

from src.utils.telemetry import get_client_ip, get_blast_limit_info, blast_limit_decorator


dash.register_page(__name__)

logger = logging.getLogger(__name__)


def blast_family_button(family):
    return dbc.Button(
        family,
        color="primary",
        href=f"/wiki?page={family}",
        external_link=False,
    )


modal = dmc.Modal(
    id="blast-modal",
    opened=False,
    centered=True,
    overlayProps={"blur": 3},
    size="lg",
    children=[
        dmc.Title(id="blast-modal-title", order=3),
        dmc.Space(h="md"),
        html.Div(id="blast-modal-content"),
    ],
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
        
        dmc.Space(h=20),
        dmc.Paper(
            children=[
                dmc.Title("BLAST Search", order=1, mb="md"),
                dmc.Text(
                    "Search protein/nucleotide sequences for Starships and Starship-associated genes",
                    c="dimmed",
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        
        dmc.Grid(
            children=[
                dmc.GridCol(
                    span={"sm": 12, "lg": 4},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
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
                                        children=dcc.Upload(
                                            id="blast-fasta-upload",
                                            children=html.Div(
                                                id="blast-fasta-sequence-upload",
                                                children="Drag and drop or click to select a FASTA file",
                                                style={"textAlign": "center", "padding": "20px"}
                                            ),
                                            multiple=False,
                                            accept=".fa, .fas, .fasta, .fna",
                                            className="upload-box text-center",
                                        ),
                                        withBorder=False,
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
                                    dmc.Text(
                                        id="rate-limit-info",
                                        size="sm",
                                        c="dimmed",
                                    ),
                                    html.Div(
                                        id="rate-limit-alert",
                                        style={"display": "none"}
                                    ),
                                ], gap="xs"),
                                
                                # Submit Section
                                dmc.Stack([
                                    dmc.Button(
                                        "Submit BLAST",
                                        id="submit-button",
                                        variant="gradient",
                                        gradient={"from": "indigo", "to": "cyan"},
                                        size="lg",
                                        fullWidth=True,
                                    ),
                                    dmc.Text(
                                        id="rate-limit-info",
                                        size="sm",
                                        c="dimmed",
                                        # align="center",
                                    ),
                                ], gap="xs"),
                                
                                # Results Preview
                                dcc.Loading(
                                    id="family-loading",
                                    type="circle",
                                    children=html.Div(id="ship-family"),
                                ),
                            ], gap="xl"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            style={"height": "100%"},  # Make paper fill grid column
                        ),
                    ],
                ),
                
                # Right Column - Results Panel
                dmc.GridCol(
                    span={"sm": 12, "lg": 8},
                    children=[
                        dmc.Paper(
                            children=[
                                # Family identification alert
                                html.Div(id="ship-family"),
                                dmc.Space(h=20),
                                # BLAST visualization
                                html.Div(id="blast-visualization"),
                                # Hidden elements for blasterjs
                                html.Div(id="blast-multiple-alignments"),
                                html.Div(id="blast-alignments-table"),
                                html.Div(id="blast-single-alignment"),
                                dcc.Upload(id="blastinput", style={"display": "none"}),
                            ],
                            p="xl",
                            radius="md",
                            withBorder=True,
                        ),
                    ],
                ),
            ],
            gutter="xl",
        ),
        modal,
    ],
)


@callback(
    Output("blast-fasta-sequence-upload", "children"),
    [
        Input("blast-fasta-upload", "contents"),
        Input("blast-fasta-upload", "filename"),
    ],
)
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
        Output("rate-limit-info", "children"),
        Output("rate-limit-alert", "children"),
        Output("rate-limit-alert", "style"),
        Output("upload-error-message", "children"),
        Output("upload-error-store", "data"),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("blast-fasta-upload", "contents")
    ],
    State("blast-fasta-upload", "filename"),
    prevent_initial_call=True
)
def handle_submission_and_upload(n_clicks, contents, filename):
    triggered_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    
    # Default return values
    button_disabled = False
    button_text = "Submit BLAST"
    limit_info = ""
    limit_alert = None
    alert_style = {"display": "none"}
    error_message = ""
    error_store = None
    
    # Handle file upload
    if triggered_id == "blast-fasta-upload" and contents is not None:
        max_size = 10 * 1024 * 1024  # 10 MB
        content_type, content_string = contents.split(",")
        
        header, seq, fasta_length_error_message = parse_fasta_from_file(contents)
        
        decoded = base64.b64decode(content_string)
        file_size = len(decoded)
        
        if fasta_length_error_message:
            error_message = dbc.Alert(f"Error: {fasta_length_error_message}", color="danger")
            error_store = error_message
            button_disabled = True
        elif file_size > max_size:
            error_message = dbc.Alert(f"Error: The file '{filename}' exceeds the 10 MB limit.", color="danger")
            error_store = error_message
            button_disabled = True
    
    # Handle rate limit check
    if triggered_id == "submit-button" and n_clicks:
        try:
            ip_address = get_client_ip()
            limit_info_data = get_blast_limit_info(ip_address)
            
            if limit_info_data["remaining"] <= 0:
                button_disabled = True
                button_text = "Rate limit exceeded"
                limit_info = f"Limit reached: {limit_info_data['submissions']}/{limit_info_data['limit']} submissions this hour"
                limit_alert = dbc.Alert(
                    [
                        html.I(className="bi bi-exclamation-triangle-fill me-2"),
                        f"Rate limit exceeded. You have used {limit_info_data['submissions']}/{limit_info_data['limit']} submissions this hour. Please try again later.",
                    ],
                    color="warning",
                    className="d-flex align-items-center",
                )
                alert_style = {"display": "block"}
            else:
                limit_info = f"Remaining: {limit_info_data['remaining']}/{limit_info_data['limit']} submissions"
        
        except Exception as e:
            logger.error(f"Error checking rate limit: {str(e)}")
    
    return [
        button_disabled,
        button_text,
        limit_info,
        limit_alert,
        alert_style,
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
)
def preprocess(n_clicks, query_text_input, query_file_contents):
    if not n_clicks:
        raise PreventUpdate

    try:
        # logger.info(
        #     f"preprocess called with n_clicks={n_clicks}, query_text_input={query_text_input}, query_file_contents={query_file_contents}"
        # )

        input_type, query_header, query_seq = check_input(
            query_text_input, query_file_contents
        )
        
        if input_type in ("none", "both"):
            logger.info("Invalid input type; returning None.")
            return None, None, None

        query_type = guess_seq_type(query_seq)
        logger.info(f"guess_seq_type returned query_type={query_type}")

        return query_header, query_seq, query_type

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        return None, None, None


@callback(
    [
        Output("blast-results-store", "data"),
        Output("captain-results-store", "data"),
    ],  # Removed Output("subject-seq-button", "children")
    [
        Input("submit-button", "n_clicks"),
        State("query-text", "value"),
        State("blast-fasta-upload", "contents"),
    ],
    prevent_initial_call=True
)
def process_blast_search(n_clicks, query_text_input, query_file_contents):
    if not n_clicks:
        return dash.no_update, dash.no_update
        
    try:
        input_type, query_header, query_seq = check_input(
            query_text_input, query_file_contents
        )
        
        if input_type in ("none", "both"):
            return dash.no_update, dash.no_update
            
        blast_results = run_blast_search(
            sequence=query_seq,
            db_path=db_list['ship']['nucl']
        )
        
        return blast_results, None
        
    except Exception as e:
        logger.error(f"Error in process_blast_search: {str(e)}")
        return dash.no_update, dash.no_update


@callback(
    Output("blast-visualization", "children"),
    Input("blast-results-store", "data"),
    prevent_initial_call=True
)
def update_blast_visualization(blast_results):
    if not blast_results:
        return dash.no_update
        
    try:
        # Convert the XML string into an array of lines and escape special characters
        blast_lines = blast_results.split('\n')
        escaped_lines = [line.replace('"', '\\"').replace("'", "\\'") for line in blast_lines]
        lines_js = ',\n'.join(f'"{line}"' for line in escaped_lines)
        
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
            <script src="/assets/js/blaster.min.js"></script>
            <script src="/assets/js/html2canvas.js"></script>
            <style>
                /* Prevent text selection on clicks */
                * {{
                    -webkit-user-select: none;
                    -moz-user-select: none;
                    -ms-user-select: none;
                    user-select: none;
                }}
            </style>
        </head>
        <body>
            <div id="blast-multiple-alignments"></div>
            <div id="blast-alignments-table"></div>
            <div id="blast-single-alignment"></div>
            
            <script>
                // Prevent default click behavior
                document.addEventListener('click', function(e) {{
                    e.preventDefault();
                    e.stopPropagation();
                }}, true);
                
                var alignments = [
                    {lines_js}
                ].join('\\n');
                
                var blasterjs = require("biojs-vis-blasterjs");
                var instance = new blasterjs({{
                    string: alignments,
                    multipleAlignments: "blast-multiple-alignments",
                    alignmentsTable: "blast-alignments-table",
                    singleAlignment: "blast-single-alignment",
                }});
            </script>
        </body>
        </html>
        """
        
        return html.Div([
            html.Iframe(
                srcDoc=html_content,
                style={
                    "width": "100%",
                    "height": "800px",
                    "border": "none",
                    "overflow": "auto"
                },
            )
        ])
    except Exception as e:
        logger.error(f"Error in update_blast_visualization: {str(e)}")
        return dash.no_update

@callback(
    Output("ship-family", "children"),
    [
        Input("blast-results-store", "data"),
        Input("captain-results-store", "data"),
        Input("curated-input", "value"),
    ],
    State("submit-button", "n_clicks"),
)
def update_ui(blast_results_xml, captain_results_dict, curated, n_clicks):
    if blast_results_xml is None and captain_results_dict is None:
        raise PreventUpdate
    if n_clicks:
        logger.info(f"Updating UI with n_clicks={n_clicks}")
        try:
            from Bio.Blast import NCBIXML
            from io import StringIO
            
            # Parse XML BLAST results
            blast_records = NCBIXML.parse(StringIO(blast_results_xml))
            record = next(blast_records)  # Get first BLAST record
            
            if record.alignments:  # If we have any hits
                best_hit = record.alignments[0]  # Get first (best) alignment
                best_hsp = best_hit.hsps[0]      # Get first (best) HSP
                
                return dmc.Alert(
                    title="BLAST Results",
                    children=[
                        f"Best hit: {best_hit.hit_def}",
                        dmc.Space(h=5),
                        dmc.Text(
                            f"Alignment length = {best_hsp.align_length}, "
                            f"Identity = {(best_hsp.identities / best_hsp.align_length) * 100:.1f}%, "
                            f"E-value = {best_hsp.expect}",
                            size="sm",
                            c="dimmed"
                        ),
                    ],
                    color="yellow",
                    variant="light",
                    withCloseButton=False,
                )

            return dmc.Alert(
                title="No Hits Found",
                children="No significant BLAST hits were found for your sequence.",
                color="yellow",
                variant="light",
                withCloseButton=False,
            )

        except Exception as e:
            logger.error(f"Error in update_ui: {str(e)}")
            return no_update

toggle_modal = create_modal_callback(
    "blast-table",
    "blast-modal",
    "blast-modal-content",
    "blast-modal-title"
)

def run_blast_search(sequence, db_path):
    # Create temporary file for query sequence
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write(sequence)
        query_file = temp_file.name

    try:
        # Modified BLAST command to output XML format
        blast_cmd = f"blastn -query {query_file} -db {db_path} -outfmt 5"
        
        process = subprocess.Popen(
            blast_cmd.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            # Return raw XML output
            return stdout.decode('utf-8')
        else:
            raise Exception(f"BLAST error: {stderr.decode('utf-8')}")
            
    finally:
        # Clean up temporary file
        os.unlink(query_file)

def init_blast_routes(server):
    """Initialize Flask routes for BLAST API"""
    
    @server.route('/blast-api', methods=['POST'])
    @blast_limit_decorator
    def blast_api():
        sequence = request.form.get('sequence', '').strip()
        if not sequence:
            return jsonify({'error': 'No sequence provided'})
            
        try:
            # Create temporary file for the sequence
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
                temp_file.write(sequence)
                query_file = temp_file.name

            try:
                # Use the ships database path from your db_list
                blast_db = db_list['ship']['nucl']
                
                # Modified BLAST command to output XML format
                blast_cmd = f"blastn -query {query_file} -db {blast_db} -outfmt 5"
                
                process = subprocess.Popen(
                    blast_cmd.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                
                stdout, stderr = process.communicate()
                
                if process.returncode == 0:
                    return stdout.decode('utf-8')  # Return raw XML
                else:
                    raise Exception(f"BLAST error: {stderr.decode('utf-8')}")
                    
            finally:
                # Clean up temporary file
                os.unlink(query_file)
                
        except Exception as e:
            logger.error(f"BLAST API error: {str(e)}")
            return jsonify({'error': str(e)})

app = dash.get_app()
init_blast_routes(app.server)
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

from dash import get_app

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
    blast_table,
    make_captain_alert,
    process_captain_results,
    create_no_matches_alert,
    blast_download_button,
)

from src.components.callbacks import curated_switch, create_modal_callback, create_file_upload
from src.database.sql_manager import fetch_meta_data
from src.config.settings import BLAST_DB_PATHS, PHYLOGENY_PATHS
from src.utils.telemetry import get_client_ip, get_blast_limit_info, blast_limit_decorator
from src.components.error_boundary import handle_callback_error, create_error_alert
from src.utils.tree import plot_tree, run_mafft, add_to_tree, gappa

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
        dcc.Store(id="processed-metadata-store"),
        dcc.Store(id="processed-blast-store"),
        
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
                                    dmc.Center(
                                        dmc.Button(
                                            "Submit BLAST",
                                            id="submit-button",
                                            variant="gradient",
                                            gradient={"from": "indigo", "to": "cyan"},
                                            size="lg",
                                        ),
                                    ),
                                    dmc.Text(
                                        id="rate-limit-info",
                                        size="sm",
                                        c="dimmed",
                                        # align="center",
                                    ),
                                ], gap="xs"),                                
                            ], gap="md"),
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
                        dcc.Loading(
                            children=[
                                dmc.Stack(
                                    children=[
                                        html.Div(id="ship-family"),
                                        html.Div(id="blast-table"),
                                        html.Div(id="blast-download"),
                                        html.Div(id="subject-seq-button"),
                                        html.Div(id="phylogeny-plot"),
                                        html.Div(
                                            id="ship-aln",
                                            style={
                                                'minHeight': '200px',
                                                'width': '100%',
                                                'marginTop': '20px'
                                            }
                                        ),
                                    ],
                                ),
                            ],
                            type="default",
                            color="#0066cc",
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
@handle_callback_error
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
        
        # Use our updated parse_fasta_from_file function
        header, seq, fasta_error = parse_fasta_from_file(contents)
        
        decoded = base64.b64decode(content_string)
        file_size = len(decoded)
        
        if fasta_error:
            error_message = fasta_error  # This will now be a dmc.Alert component
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
@handle_callback_error
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
        logger.info(
            f"check_input returned input_type={input_type}, query_header={query_header}, query_seq={query_seq}"
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
        Output("subject-seq-button", "children"),
        Output("phylogeny-plot", "children"),
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def fetch_captain(query_header, query_seq, query_type, search_type="hmmsearch"):
    if not all([query_header, query_seq, query_type]):
        return None, None, None, None
        
    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        logger.info(f"Temp FASTA written: {tmp_query_fasta}")

        # Run BLAST
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name

        try:
            blast_results = run_blast(
                db_list=BLAST_DB_PATHS,
                query_type=query_type,
                query_fasta=tmp_query_fasta,
                tmp_blast=tmp_blast,
                input_eval=0.01,
                threads=2,
            )
            logger.info(f"BLAST results: {blast_results.head()}")
            if blast_results is None:
                raise ValueError("BLAST returned no results!")
        except Exception as e:
            logger.error(f"BLAST error: {str(e)}")
            raise

        blast_results_dict = blast_results.to_dict("records")

        # Run HMMER
        logger.info(f"Running HMMER")

        subject_seq_button = None
        # subject_seq = None

        try:
            if search_type == "diamond":
                results_dict = run_diamond(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
            if search_type == "hmmsearch":
                results_dict = run_hmmer(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )

            if results_dict is None or len(results_dict) == 0:
                logger.error("Diamond/HMMER returned no results!")
                raise
        except Exception as e:
            logger.error(f"Diamond/HMMER error: {str(e)}")
            raise

        # Run phylogenetic placement if protein sequence
        phylogeny_plot = None
        if query_type == "prot" and results_dict:
            try:
                with tempfile.TemporaryDirectory() as tmp_dir:
                    # Write query to temporary file
                    tmp_query_fasta = write_temp_fasta(query_header, query_seq)
                    
                    logger.info(f"Running MAFFT alignment for query: {query_header}")
                    # Run MAFFT alignment
                    tmp_headers, query_msa = run_mafft(
                        query=tmp_query_fasta,
                        ref_msa=PHYLOGENY_PATHS["msa"]
                    )
                    
                    logger.info("Running EPA-NG placement")
                    # Run EPA-NG placement
                    new_tree = add_to_tree(
                        query_msa=query_msa,
                        tree=PHYLOGENY_PATHS["tree"],
                        ref_msa=PHYLOGENY_PATHS["msa"],
                        model="PROTGTR+G+F",
                        tmp_dir=tmp_dir
                    )
                    
                    logger.info("Running Gappa to get final tree")
                    # Run Gappa to get final tree
                    out_tree = gappa(new_tree)
                    
                    # Create phylogeny plot
                    if out_tree is not None:
                        logger.info("Creating phylogeny plot")
                        fig = plot_tree(
                            highlight_families="all",
                            tips=[query_header]  # Add the query sequence header
                        )
                        
                        phylogeny_plot = html.Div([
                            dmc.Paper(
                                children=[
                                    dmc.Title("Phylogenetic Placement", order=3),
                                    dcc.Graph(
                                        id="phylogeny-graph",
                                        figure=fig,
                                        style={'height': '800px'}  # Set explicit height
                                    )
                                ],
                                p="md",
                                withBorder=True,
                                shadow="sm",
                                radius="md",
                                style={'marginTop': '20px'}
                            )
                        ])
                        
                        logger.info("Phylogeny plot created successfully")
                    else:
                        logger.warning("No tree output from Gappa")
                        
            except Exception as e:
                logger.error(f"Error in phylogenetic placement: {str(e)}")
                phylogeny_plot = html.Div([
                    dmc.Alert(
                        title="Phylogenetic Placement Error",
                        color="yellow",
                        children=f"Could not generate phylogenetic tree: {str(e)}"
                    )
                ])

        return blast_results_dict, results_dict, subject_seq_button, phylogeny_plot

    except Exception as e:
        logger.error(f"Error in fetch_captain: {str(e)}")
        error_div = html.Div([
            dmc.Alert(
                title="Error",
                color="red",
                children=str(e)
            )
        ])
        return None, None, error_div, None


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
@cache.memoize(timeout=86400)
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
    [Input("blast-results-store", "data"), Input("processed-metadata-store", "data")],
)
@handle_callback_error
def process_blast_results(blast_results_dict, metadata_dict):
    if not blast_results_dict or not metadata_dict:
        return None
        
    try:
        # Add size limit check
        results_size = len(str(blast_results_dict))
        if results_size > 5 * 1024 * 1024:  # 5MB limit
            logger.warning(f"BLAST results too large: {results_size} bytes")
            return None

        blast_df = pd.DataFrame(blast_results_dict)
        meta_df = pd.DataFrame(metadata_dict)
                
        # Process in chunks if dataset is large
        chunk_size = 1000
        if len(blast_df) > chunk_size:
            chunks = []
            for i in range(0, len(blast_df), chunk_size):
                chunk = blast_df.iloc[i:i + chunk_size]
                processed_chunk = pd.merge(
                    meta_df[["accession_tag", "familyName"]],
                    chunk,
                    left_on="accession_tag",
                    right_on="sseqid",
                    how="right",
                    suffixes=('', '_blast')
                )
                chunks.append(processed_chunk)
            df_for_table = pd.concat(chunks)
        else:
            df_for_table = pd.merge(
                meta_df[["accession_tag", "familyName"]],
                blast_df,
                left_on="accession_tag",
                right_on="sseqid",
                how="right",
                suffixes=('', '_blast')
            )
                
        df_for_table = (df_for_table
            .drop_duplicates(subset=["accession_tag", "pident", "length"])
            .dropna(subset=["accession_tag"])
            .fillna("")
            .sort_values("evalue"))
        
        # Limit number of results
        df_for_table = df_for_table.head(1000)  # Only return top 1000 matches
        
        return df_for_table.to_dict('records') if not df_for_table.empty else None
        
    except Exception as e:
        logger.error(f"Error processing BLAST results: {str(e)}")
        return None

# 3. UI Update Callback
@callback(
    [
        Output("ship-family", "children"),
        Output("blast-table", "children"),
        Output("blast-download", "children")
    ],
    [Input("processed-blast-store", "data"), Input("captain-results-store", "data")],
    State("submit-button", "n_clicks"),
    prevent_initial_call=True
)
@handle_callback_error
def update_ui_elements(processed_blast_results, captain_results_dict, n_clicks):
    if not n_clicks or processed_blast_results is None:
        return None, None, None
        
    try:
        if not processed_blast_results:
            return html.Div(create_no_matches_alert()), None, None
            
        try:
            blast_df = pd.DataFrame(processed_blast_results)
            if blast_df.empty:
                return html.Div(create_no_matches_alert()), None, None
        except Exception as e:
            logger.error(f"Error converting BLAST results to DataFrame: {e}")
            return html.Div(create_error_alert(str(e))), None, None
        
        try:
            best_match = blast_df.nsmallest(1, 'evalue').iloc[0]
        except Exception as e:
            logger.error(f"Error getting best match: {e}")
            return html.Div(create_error_alert("Could not determine best match")), None, None
        
        try:
            ship_family = (
                make_captain_alert(
                    best_match["familyName"], 
                    best_match["length"], 
                    best_match["evalue"], 
                    search_type="blast"
                ) if best_match["pident"] > 50
                else process_captain_results(captain_results_dict)
            )
        except Exception as e:
            logger.error(f"Error creating family alert: {e}")
            return html.Div(create_error_alert("Could not create family alert")), None, None
            
        try:
            table = blast_table(blast_df)
            download_button = blast_download_button()
        except Exception as e:
            logger.error(f"Error creating BLAST table: {e}")
            return ship_family, create_error_alert("Could not create results table"), None
        
        return ship_family, table, download_button
        
    except Exception as e:
        logger.error(f"Error in update_ui_elements: {e}")
        return create_error_alert(str(e)), None, None

@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("blast-table", "derived_virtual_data")],
    prevent_initial_call=True
)
@handle_callback_error
def download_tsv(n_clicks, table_data):
    if not n_clicks:
        return None
        
    if not table_data:
        logger.error("Error: No data available for download.")
        return None

    try:
        df = pd.DataFrame(table_data)
        
        # Format the data for download
        download_columns = ['accession_tag', 'familyName', 'pident', 'length', 'evalue', 'bitscore']
        df = df[download_columns]
        
        tsv_string = df.to_csv(sep="\t", index=False)
        tsv_bytes = io.BytesIO(tsv_string.encode())
        b64 = base64.b64encode(tsv_bytes.getvalue()).decode()

        today = date.today().strftime("%Y-%m-%d")

        return dict(
            content=f"data:text/tab-separated-values;base64,{b64}",
            filename=f"starbase_blast_{today}.tsv",
        )

    except Exception as e:
        logger.error(f"Error preparing download: {str(e)}")
        return None

toggle_modal = create_modal_callback(
    "blast-table",
    "blast-modal",
    "blast-modal-content",
    "blast-modal-title"
)
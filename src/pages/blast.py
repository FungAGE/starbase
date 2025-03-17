import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import dcc, html, callback
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State

import io
import tempfile
import base64
from datetime import date
import pandas as pd
import logging
import subprocess


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
        
        dcc.Store(id="timeout-store", data=False),
        dcc.Interval(id="timeout-interval", interval=300000, n_intervals=0),  # 5 minutes
        
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
                                            loading=False,
                                            loaderProps={"variant": "dots", "color": "white"},
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
                                html.Div(id="blast-table"),
                                html.Div(id="blast-download"),
                                html.Div(id="subject-seq-button"),
                            ],
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
        Input("blast-fasta-upload", "contents"),
        Input("timeout-interval", "n_intervals")  # Add timeout interval input
    ],
    [
        State("blast-fasta-upload", "filename"),
        State("timeout-store", "data")  # Add timeout store state
    ],
    prevent_initial_call=True
)
@handle_callback_error
def handle_submission_and_upload(n_clicks, contents, n_intervals, filename, timeout_triggered):
    triggered_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
    
    # Default return values
    button_disabled = False
    button_text = "Submit BLAST"
    limit_info = ""
    limit_alert = None
    alert_style = {"display": "none"}
    error_message = ""
    error_store = None
    
    # Handle timeout case
    if triggered_id == "timeout-interval":
        # Only show timeout if button was clicked and we're still waiting for results
        if n_clicks and timeout_triggered and not dash.callback_context.inputs.get("ship-family.children"):
            limit_alert = dbc.Alert(
                [
                    html.I(className="bi bi-exclamation-triangle-fill me-2"),
                    "The server is taking longer than expected to respond. Please try again later.",
                ],
                color="warning",
                className="d-flex align-items-center",
            )
            alert_style = {"display": "block"}
            button_disabled = True
            button_text = "Server timeout"
            return [
                button_disabled,
                button_text,
                limit_info,
                limit_alert,
                alert_style,
                error_message,
                error_store
            ]
        # If button wasn't clicked or results are already shown, don't show timeout
        return [
            button_disabled,
            button_text,
            limit_info,
            limit_alert,
            {"display": "none"},
            error_message,
            error_store
        ]
    
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
        # Remove timeout wrapper since it uses signals
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
    [
        Output("blast-results-store", "data"),
        Output("captain-results-store", "data"),
        Output("subject-seq-button", "children"),
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
    ],
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def fetch_captain(query_header, query_seq, query_type, search_type="hmmsearch"):
    if not all([query_header, query_seq, query_type]):
        return None, None, None
        
    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        
        # Instead of using timeout context manager, we'll rely on subprocess timeout
        try:
            blast_results = run_blast(
                db_list=BLAST_DB_PATHS,
                query_type=query_type,
                query_fasta=tmp_query_fasta,
                tmp_blast=tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name,
                input_eval=0.01,
                threads=2,
            )
            
            if blast_results is None:
                error_div = html.Div([
                    dmc.Alert(
                        title="BLAST Error",
                        color="red",
                        children="No BLAST results were returned. Please try again with a different sequence.",
                    )
                ])
                return None, None, error_div

            blast_results_dict = blast_results.to_dict("records")
            
            # Run HMMER/Diamond with timeout
            if search_type == "diamond":
                results_dict = run_diamond(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
            else:
                results_dict = run_hmmer(
                    db_list=BLAST_DB_PATHS,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
                
            return blast_results_dict, results_dict, None

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
        Output("blast-download", "children"),
        Output("phylogeny-plot", "children")
    ],
    [Input("processed-blast-store", "data"), 
     Input("captain-results-store", "data")],
    State("submit-button", "n_clicks"),
    running=[
        (Output("submit-button", "loading"), True, False),
        (Output("submit-button", "disabled"), True, False),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def update_ui_elements(processed_blast_results, captain_results_dict, n_clicks):
    if not n_clicks or processed_blast_results is None:
        return None, None, None, None
        
    try:
        if not processed_blast_results:
            return html.Div(create_no_matches_alert()), None, None, None
            
        try:
            blast_df = pd.DataFrame(processed_blast_results)
            if blast_df.empty:
                return html.Div(create_no_matches_alert()), None, None, None
        except Exception as e:
            logger.error(f"Error converting BLAST results to DataFrame: {e}")
            return html.Div(create_error_alert(str(e))), None, None, None
        
        ship_family = None
        table = None
        download_button = None
        phylogeny_plot = None
        
        try:
            best_match = blast_df.nsmallest(1, 'evalue').iloc[0]
            
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
            ship_family = html.Div(create_error_alert(str(e)))
            
        try:
            # Create phylogeny plot
            logger.info("Creating phylogeny plot")
            fig = plot_tree(
                highlight_families=best_match["familyName"],
                tips=[best_match["accession_tag"]]
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
        
        try:
            blast_table_content = blast_table(blast_df)
            download_button = blast_download_button()
            
            # Wrap table and download button in accordion
            table = html.Div([
                dbc.Accordion(
                    children=[
                        dbc.AccordionItem(
                            children=[
                                blast_table_content,
                                html.Div(
                                    download_button,
                                    style={"margin-top": "1rem"}
                                )
                            ],
                            title="BLAST Results",
                            item_id="blast-results",
                        )
                    ],
                    always_open=False,
                    id="blast-accordion"
                )
            ])
            # Set download_button to None since it's now included in the table accordion
            download_button = None
            
        except Exception as e:
            logger.error(f"Error creating BLAST table: {e}")
            table = html.Div(create_error_alert("Could not create results table"))
            download_button = None
        
        return ship_family, table, download_button, phylogeny_plot
        
    except Exception as e:
        logger.error(f"Error in update_ui_elements: {e}")
        error_div = html.Div(create_error_alert(str(e)))
        return error_div, None, None, None

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

@callback(
    Output("timeout-store", "data"),
    [Input("timeout-interval", "n_intervals")],
    [State("submit-button", "n_clicks")]
)
def check_timeout(n_intervals, n_clicks):
    if n_clicks and n_intervals > 0:
        return True
    return False
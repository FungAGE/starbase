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
    make_captain_alert
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
from src.components.callbacks import curated_switch, create_modal_callback, create_file_upload
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
                            children=dmc.Stack([
                                dmc.Title("BLAST Results", order=3),
                                dcc.Loading(
                                    id="family-loading",
                                    type="circle",
                                    children=html.Div(id="ship-family"),
                                ),
                                dcc.Loading(
                                    id="blast-table-loading",
                                    type="circle",
                                    children=html.Div(id="blast-table"),
                                ),
                                dcc.Loading(
                                    id="subject-seq-button-loading",
                                    type="circle",
                                    children=html.Div(id="subject-seq-button"),
                                ),
                                dcc.Loading(
                                    id="ship-aln-loading",
                                    type="circle",
                                    children=html.Div(
                                        id="ship-aln",
                                        style={
                                            'minHeight': '200px',
                                            'width': '100%',
                                            'marginTop': '20px'
                                        }
                                    ),
                                ),
                            ], gap="xl"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            style={"height": "100%"},  # Make paper fill grid column
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
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
    ],
)
def fetch_captain(query_header, query_seq, query_type, search_type="hmmsearch"):
    try:
        if not query_header or not query_seq:
            logger.error("Missing query header or sequence.")
            return None, None, None

        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        logger.info(f"Temp FASTA written: {tmp_query_fasta}")

        # Run BLAST
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name

        try:
            blast_results = run_blast(
                db_list=db_list,
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
                    db_list=db_list,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
            if search_type == "hmmsearch":
                results_dict = run_hmmer(
                    db_list=db_list,
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

        return blast_results_dict, results_dict, subject_seq_button

    except Exception as e:
        logger.error(f"Error in fetch_captain: {str(e)}")
        return None, None, None


@callback(
    Output("subject-seq-dl-package", "data"),
    [Input("subject-seq-button", "n_clicks"), Input("subject-seq", "data")],
)
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


no_captain_alert = dbc.Alert(
    "No captain sequence found (e-value threshold 0.01).",
    color="warning",
)


@callback(
    [
        Output("ship-family", "children"),
        Output("blast-table", "children"),
    ],
    [
        Input("blast-results-store", "data"),
        Input("captain-results-store", "data"),
        Input("curated-input", "value"),
    ],
    State("submit-button", "n_clicks"),
)
def update_ui(blast_results_dict, captain_results_dict, curated, n_clicks):
    if not n_clicks:
        raise PreventUpdate

    if blast_results_dict is None and captain_results_dict is None:
        return [
            dmc.Alert(
                title="Search Error",
                children=[
                    "No results were returned from the BLAST search.",
                    dmc.Space(h=5),
                    dmc.Text(
                        "This could be due to:",
                        size="sm",
                        mb=5,
                    ),
                    dmc.List(
                        [
                            "Invalid sequence format",
                            "No significant matches found",
                            "Server-side processing error",
                        ],
                        size="sm",
                    ),
                ],
                color="red",
                variant="light",
                withCloseButton=False,
            ),
            None
        ]

    logger.info(f"Updating UI with n_clicks={n_clicks}")
    try:
        ship_family = no_update
        ship_table = no_update

        blast_results_df = pd.DataFrame(blast_results_dict)

        # TODO: caching the curated dataset makes no sense. filter the full dataset based on curated flag after loading from cache.
        initial_df = cache.get("meta_data")

        if initial_df is None:
            initial_df = fetch_meta_data(curated)

        initial_df = initial_df[["accession_tag", "familyName"]].drop_duplicates()

        if blast_results_dict:
            logger.info("Rendering BLAST table")
            df_for_table = pd.merge(
                initial_df,
                blast_results_df,
                left_on="accession_tag",
                right_on="sseqid",
                how="right",
            )
            
            if len(df_for_table) == 0:
                return [
                    dmc.Alert(
                        title="No Matches Found",
                        children=[
                            "Your sequence did not match any Starships in our database.",
                            dmc.Space(h=5),
                            dmc.Text(
                                "Suggestions:",
                                size="sm",
                                mb=5,
                            ),
                            dmc.List(
                                [
                                    "Check if your sequence is in the correct format",
                                    "Try searching with a different region of your sequence",
                                    "Consider using a less stringent E-value threshold",
                                ],
                                size="sm",
                            ),
                        ],
                        color="yellow",
                        variant="light",
                        withCloseButton=False,
                    ),
                    None
                ]

            # we want to remove true duplicates, while keeping other hits which may be at a different location in the same ship
            df_for_table = df_for_table.drop_duplicates(
                subset=["accession_tag", "pident", "length"]
            )
            df_for_table = df_for_table[df_for_table["accession_tag"].notna()]
            df_for_table.fillna("", inplace=True)

            if len(df_for_table) > 0:
                ship_table = blast_table(df_for_table)
            else:
                ship_table = dbc.Alert(
                    "No BLAST results found.",
                    color="danger",
                )

            min_evalue_rows = df_for_table.loc[
                df_for_table.groupby("qseqid")["evalue"].idxmin()
            ]
            if min_evalue_rows.empty:
                logger.warning(
                    "min_evalue_rows is empty after grouping by qseqid and selecting min evalue."
                )
            else:
                logger.info(f"min_evalue_rows contents: {min_evalue_rows}")

            if not min_evalue_rows.empty and "pident" in min_evalue_rows.columns:
                if min_evalue_rows["pident"].iloc[0] > 95:
                    family_name = min_evalue_rows["familyName"].iloc[0]
                    aln_len = min_evalue_rows["length"].iloc[0]
                    ev = min_evalue_rows["evalue"].iloc[0]
                    ship_family = make_captain_alert(family_name, aln_len, ev)
                else:
                    if captain_results_dict:
                        logger.info("Processing Diamond/HMMER results")
                        captain_results_df = pd.DataFrame(captain_results_dict)
                        # captain_results_df["sseqid"] = captain_results_df["sseqid"].apply(
                        #     clean_shipID
                        # )
                        if len(captain_results_df) > 0:
                            try:
                                superfamily, family_aln_length, family_evalue = (
                                    select_ship_family(captain_results_df)
                                )
                                if superfamily:
                                    family = initial_df[
                                        initial_df["familyName"] == superfamily
                                    ]["familyName"].unique()[0]
                                    if family:
                                        ship_family = make_captain_alert(family, family_aln_length, family_evalue)
                                    else:
                                        ship_family = dmc.Alert(
                                            title="Starship Family Not Determined",
                                            children=[
                                                "Starship family could not be determined.",
                                                dmc.Space(h=5),
                                                dmc.Text(
                                                    "Please try a different query or increase the e-value threshold.",
                                                    size="sm",
                                                    c="dimmed"
                                                ),
                                            ],
                                            color="red",
                                            variant="light",
                                            withCloseButton=False,
                                        )

                            except Exception as e:
                                logger.error(
                                    f"Error selecting ship family: {str(e)}"
                                )
                                ship_family = html.Div(f"Error: {str(e)}")
                        else:
                            ship_family = no_captain_alert
                    else:
                        ship_family = no_captain_alert
            else:
                ship_family = dbc.Alert(
                    "No matching rows found for minimum evalue selection.",
                    color="danger",
                )
        return ship_family, ship_table

    except Exception as e:
        logger.error(f"Error in update_ui: {str(e)}")
        return [
            dmc.Alert(
                title="Processing Error",
                children=[
                    "An error occurred while processing your search.",
                    dmc.Space(h=5),
                    dmc.Text(
                        f"Error details: {str(e)}",
                        size="sm",
                        c="dimmed",
                    ),
                    dmc.Space(h=5),
                    dmc.Text(
                        "Please try again or contact support if the problem persists.",
                        size="sm",
                    ),
                ],
                color="red",
                variant="light",
                withCloseButton=False,
            ),
            None
        ]


@cache.memoize()
@callback(
    Output("ship-aln", "children"),
    [
        Input("blast-table", "derived_virtual_data"),
        Input("blast-table", "derived_virtual_selected_rows"),
        Input("blast-table", "selected_rows"),
        Input("curated-input", "value"),
    ],
    State("query-type-store", "data"),
)
def blast_alignments(ship_blast_results, derived_selected_rows, selected_rows, curated, query_type):
    # Use whichever selection is available
    active_selection = derived_selected_rows if derived_selected_rows else selected_rows
    
    try:
        if not active_selection or len(active_selection) == 0:
            return None

        logger.info(
            f"blast_alignments called with selection: {active_selection}, query_type: {query_type}"
        )
        
        if not ship_blast_results or len(ship_blast_results) == 0:
            logger.error(
                "No BLAST results available because ship_blast_results is empty or None."
            )
            raise

        ship_blast_results_df = pd.DataFrame(ship_blast_results)

        row_idx = active_selection[0]

        try:
            row = ship_blast_results_df.iloc[row_idx]
            qseq = str(row["qseq"])
            qseqid = str(row["qseqid"])
            sseq = str(row["sseq"])
            sseqid = str(row["sseqid"])

        except IndexError:
            logger.error(f"Error: Row index {row_idx} out of bounds.")
            raise
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)

        try:
            with open(tmp_fasta.name, "w") as f:
                f.write(f">{qseqid}\n{qseq}\n>{sseqid}\n{sseq}\n")
        except Exception as file_error:
            logger.error(f"Error writing to FASTA file: {file_error}")
            raise

        try:
            with open(tmp_fasta.name, "r") as file:
                data = file.read()
        except Exception as read_error:
            logger.error(f"Error reading FASTA file: {read_error}")
            raise

        color = "nucleotide" if query_type == "nucl" else "clustal2"

        aln = dashbio.AlignmentChart(
            id="alignment-viewer",
            data=data,
            height=200,
            tilewidth=30,
            colorscale=color,
            showconsensus=False,
            showconservation=False,
            showgap=False,
            showid=False,
            ticksteps=5,
            tickstart=0,
        )

        # Add download button below alignment
        download_section = html.Div([
            aln,
            dmc.Space(h="md"),
            dmc.Button(
                "Download Alignment FASTA",
                id="download-alignment-button",
                leftSection=[html.I(className="bi bi-download")],
                variant="gradient",
                gradient={"from": "indigo", "to": "cyan"},
                size="lg",
            ),
            dcc.Download(id="download-alignment-fasta")
        ])

        return download_section

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        raise

@callback(
    Output("download-alignment-fasta", "data"),
    Input("download-alignment-button", "n_clicks"),
    [
        State("blast-table", "derived_virtual_data"),
        State("blast-table", "derived_virtual_selected_rows"),
        State("blast-table", "selected_rows"),
    ],
    prevent_initial_call=True
)
def download_alignment_fasta(n_clicks, ship_blast_results, derived_selected_rows, selected_rows):
    if not n_clicks:
        return None
        
    active_selection = derived_selected_rows if derived_selected_rows else selected_rows
    
    if not active_selection or len(active_selection) == 0:
        return None
        
    try:
        ship_blast_results_df = pd.DataFrame(ship_blast_results)
        row = ship_blast_results_df.iloc[active_selection[0]]
        
        qseq = str(row["qseq"]).replace("\n", "")
        qseqid = str(row["qseqid"])
        sseq = str(row["sseq"]).replace("\n", "")
        sseqid = str(row["sseqid"])
        
        fasta_content = f">{qseqid}\n{qseq}\n>{sseqid}\n{sseq}\n"
        
        return dict(
            content=fasta_content,
            filename=f"alignment_{qseqid}_{sseqid}.fasta",
            type="text/plain"
        )
        
    except Exception as e:
        logger.error(f"Error preparing alignment download: {str(e)}")
        return None

@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("blast-table", "derived_virtual_data")],
    prevent_initial_call=True
)
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

@callback(
    Output("lastz-plot", "figure"),
    [
        Input("blast-table", "derived_virtual_data"),
        Input("blast-table", "derived_virtual_selected_rows"),
    ],
)
def create_alignment_plot(ship_blast_results, selected_row):
    """
    Creates a Plotly scatter plot from the alignment DataFrame and saves it as a PNG.
    """
    tmp_fasta_clean = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)
    lastz_output = tempfile.NamedTemporaryFile(suffix=".tsv", delete=True)

    try:
        ship_blast_results_df = pd.DataFrame(ship_blast_results)
        if ship_blast_results_df.empty:
            logger.error("Error: No blast results available for plotting.")
            return None
    except Exception as e:
        logger.error(f"Error converting blast results to DataFrame: {e}")
        return None

    if selected_row is None or not selected_row:
        logger.error("No row selected for alignment.")
        return None

    try:
        row = ship_blast_results_df.iloc[selected_row]
        qseq = re.sub("-", "", row["qseq"])
        qseqid = row["qseqid"]
        sseq = re.sub("-", "", row["sseq"])
        sseqid = row["sseqid"]

        logger.info(f"Selected query ID: {qseqid}, subject ID: {sseqid}")

        with open(tmp_fasta_clean.name, "w") as f:
            f.write(f">{qseqid}\n{qseq}\n>{sseqid}\n{sseq}\n")

        logger.info("Running LASTZ...")
        run_lastz(tmp_fasta_clean.name, lastz_output.name)

        lastz_df = parse_lastz_output(lastz_output.name)

        if lastz_df.empty:
            logger.error("Error: No alignment data from LASTZ.")
            return None

        x_values = []
        y_values = []

        for _, lastz_row in lastz_df.iterrows():
            x_values.extend([lastz_row["qstart"], lastz_row["qend"]])
            y_values.extend([lastz_row["sstart"], lastz_row["send"]])

        fig = go.Figure(data=go.Scatter(x=x_values, y=y_values, mode="lines"))

        fig.update_layout(
            title="LASTZ Alignment",
            xaxis_title=f"Query: {qseqid}",
            yaxis_title=f"Subject: {sseqid}",
        )

        return fig

    except IndexError as e:
        logger.error(f"Error: Selected row index is out of bounds. Details: {e}")
        return None

    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return None


toggle_modal = create_modal_callback(
    "blast-table",
    "blast-modal",
    "blast-modal-content",
    "blast-modal-title"
)
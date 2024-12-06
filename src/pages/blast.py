import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import dcc, html, callback, MATCH, no_update
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State
import dash_bio as dashbio

import io
import os
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
    blast2html,
    parse_blast_pairwise,
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
                            children=dmc.Stack([
                                dmc.Title("BLAST Results", order=3),
                                dcc.Loading(
                                    id="ship-blast-table-loading",
                                    type="circle",
                                    children=html.Div(id="ship-blast-table"),
                                ),
                                dcc.Loading(
                                    id="subject-seq-button-loading",
                                    type="circle",
                                    children=html.Div(id="subject-seq-button"),
                                ),
                                dcc.Loading(
                                    id="ship-aln-loading",
                                    type="circle",
                                    children=html.Div(id="ship-aln"),
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
        # logger.info(
        #     f"check_input returned input_type={input_type}, query_header={query_header}, query_seq={query_seq}"
        # )

        if input_type in ("none", "both"):
            logger.info("Invalid input type; returning None.")
            return None, None, None

        query_type = guess_seq_type(query_seq)
        logger.info(f"guess_seq_type returned query_type={query_type}")

        return query_header, query_seq, query_type

    except Exception as e:
        logger.error(f"Error in preprocess: {str(e)}")
        return None, None, None


def captain_family_classification(
    query_header, query_seq, query_type, tmp_query_fasta, search_type="hmmsearch"
):
    try:
        if not query_header or not query_seq:
            logger.error("Missing query header or sequence.")
            return None, None, None

        logger.info(f"Running HMMER")

        subject_seq_button = None
        try:
            if search_type == "diamond":
                classification_results_dict = run_diamond(
                    db_list=db_list,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )
            if search_type == "hmmsearch":
                classification_results_dict = run_hmmer(
                    db_list=db_list,
                    query_type=query_type,
                    input_genes="tyr",
                    input_eval=0.01,
                    query_fasta=tmp_query_fasta,
                    threads=2,
                )

            if (
                classification_results_dict is None
                or len(classification_results_dict) == 0
            ):
                logger.error("Diamond/HMMER returned no results!")
                raise
        except Exception as e:
            logger.error(f"Diamond/HMMER error: {str(e)}")
            raise
    except Exception as e:
        logger.error(f"Error in fetch_captain: {str(e)}")
        return None
    return classification_results_dict, subject_seq_button


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
def blast(query_header, query_seq, query_type):
    tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name
    logger.info(f"Temp BLAST written: {tmp_blast}")
    tmp_query_fasta = write_temp_fasta(query_header, query_seq)
    logger.info(f"Temp FASTA written: {tmp_query_fasta}")

    try:
        blast_out_file = run_blast(
            db_list=db_list,
            query_type=query_type,
            query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
            outfmt="pairwise",
        )
        if isinstance(blast_out_file, str) and os.path.exists(blast_out_file):
            blast_htmls = blast2html(input=blast_out_file)
        else:
            raise ValueError("BLAST returned no results!")
            # raise ValueError(
            #     "Invalid input: blast_results must be a DataFrame or a valid file path."
            # )
    except Exception as e:
        logger.error(f"BLAST error: {str(e)}")
        raise

    # blast_results_dict = parse_blast_pairwise(blast_out_file).to_dict("records")

    classification_results_dict, subject_seq_button = captain_family_classification(
        query_header, query_seq, query_type, tmp_query_fasta, search_type="hmmsearch"
    )
    return (
        # blast_results_dict,
        blast_htmls,
        classification_results_dict,
        subject_seq_button,
    )


@callback(
    [
        Output("ship-family", "children"),
        Output("blast-results-container", "children"),

    ],
    [
        Input("blast-results-store", "data"),
        Input("captain-results-store", "data"),
        Input("curated-input", "value"),
    ],
    State("submit-button", "n_clicks"),
)
def update_ui(blast_results_dict, captain_results_dict, curated, n_clicks):
    if blast_results_dict is None and captain_results_dict is None:
        raise PreventUpdate
        
    logger.info(f"Updating UI with n_clicks={n_clicks}")
    
    try:
        ship_family = no_update
        blast_viewer = no_update
        
        no_captain_alert = dbc.Alert(
            "No captain sequence found (e-value threshold 0.01).",
            color="warning",
        )

        # Convert blast results to DataFrame
        if isinstance(blast_results_dict, dict):
            # Convert single result to list format
            blast_results_dict = {k: [v] if not isinstance(v, list) else v 
                                for k, v in blast_results_dict.items()}
        blast_results_df = pd.DataFrame(blast_results_dict)
        
        # Load metadata
        initial_df = load_from_cache("meta_data")
        if initial_df is None:
            initial_df = fetch_meta_data(curated)
        initial_df = initial_df[["accession_tag", "familyName"]].drop_duplicates()
        
        df_for_table = pd.merge(
            initial_df,
            blast_results_df,
            left_on="accession_tag",
            right_on="sseqid",
            how="right",
        )
        
        # Remove duplicates while keeping other hits
        df_for_table = df_for_table.drop_duplicates(
            subset=["accession_tag", "pident", "length"]
        )
        df_for_table = df_for_table[df_for_table["accession_tag"].notna()]
        df_for_table.fillna("", inplace=True)

        if blast_results_dict:
            blast_viewer = BlastViewer(blast_results_dict)
        else:
            blast_viewer = dbc.Alert(
                "No BLAST results found.",
                color="danger",
            )
            
        # Process min evalue rows
        min_evalue_rows = df_for_table.loc[
            df_for_table.groupby("qseqid")["evalue"].idxmin()
        ]
        
        if min_evalue_rows.empty:
            logger.warning(
                "min_evalue_rows is empty after grouping by qseqid and selecting min evalue."
            )
        else:
            if not min_evalue_rows.empty and "pident" in min_evalue_rows.columns:
                if min_evalue_rows["pident"].iloc[0] > 95:
                    family_name = min_evalue_rows["familyName"].iloc[0]
                    aln_len = min_evalue_rows["length"].iloc[0]
                    ev = min_evalue_rows["evalue"].iloc[0]
                    ship_family = dbc.Alert(
                        [
                            f"Your sequence is likely in Starship family: {family_name} (Alignment length = {aln_len}, evalue = {ev})",
                        ],
                        color="warning",
                    )
                else:
                    # Process captain results if available
                    if captain_results_dict:
                        logger.info("Processing Diamond/HMMER results")
                        # Same conversion for captain results if needed
                        if isinstance(captain_results_dict, dict):
                            captain_results_dict = {k: [v] if not isinstance(v, list) else v 
                                                  for k, v in captain_results_dict.items()}
                        captain_results_df = pd.DataFrame(captain_results_dict)
                        
                        if len(captain_results_df) > 0:
                            try:
                                superfamily, family_aln_length, family_evalue = select_ship_family(captain_results_df)
                                if superfamily:
                                    family = initial_df[
                                        initial_df["familyName"] == superfamily
                                    ]["familyName"].unique()[0]
                                    if family:
                                        ship_family = dbc.Alert(
                                            [
                                                f"Your sequence is likely in Starship family: {family} (Alignment length = {family_aln_length}, evalue = {family_evalue})",
                                            ],
                                            color="warning",
                                        )
                                    else:
                                        ship_family = dbc.Alert(
                                            ["Starship family could not be determined."],
                                            color="danger",
                                        )
                            except Exception as e:
                                logger.error(f"Error selecting ship family: {str(e)}")
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
                
        return ship_family, blast_viewer
        
    except Exception as e:
        logger.error(f"Error in update_ui: {str(e)}")
        logger.exception(e)  # This will log the full traceback
        return no_update, no_update


@callback(
    [
        Output("upload-error-message", "children"),
        Output("submit-button", "disabled"),
    ],
    Input("blast-fasta-upload", "contents"),
    State("blast-fasta-upload", "filename"),
    prevent_initial_call=True,
)
def handle_fasta_upload(contents, filename):
    if contents is None:
        return "", None
    max_size = 10 * 1024 * 1024  # 10 MB in bytes

    content_type, content_string = contents.split(",")

    header, seq, fasta_length_error_message = parse_fasta_from_file(contents)

    decoded = base64.b64decode(content_string)
    file_size = len(decoded)

    if fasta_length_error_message:
        error_message = dbc.Alert(
            f"Error: {fasta_length_error_message}", color="danger"
        )
        return error_message, True
    elif file_size > max_size:
        error_message = dbc.Alert(
            f"Error: The file '{filename}' exceeds the 10 MB limit.", color="danger"
        )
        return error_message, True
    else:
        return "", False

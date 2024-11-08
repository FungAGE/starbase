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

from src.components.cache import cache
from src.utils.blast_utils import (
    check_input,
    guess_seq_type,
    write_temp_fasta,
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
from src.components.callbacks import curated_switch, create_accession_modal
from src.utils.parsing import parse_fasta, parse_fasta_from_file
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import fetch_meta_data
from src.utils.blastdb import db_list

import logging

dash.register_page(__name__)

logger = logging.getLogger(__name__)


def blast_family_button(family):
    return dbc.Button(
        family,
        color="primary",
        href=f"/wiki?page={family}",
        external_link=False,
    )


modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle(id="blast-modal-title")),
        dbc.ModalBody(id="blast-modal-content"),
    ],
    id="blast-modal",
    is_open=False,
)


layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="query-header-store"),
        dcc.Store(id="query-seq-store"),
        dcc.Store(id="query-type-store"),
        dcc.Store(id="blast-results-store"),
        dcc.Store(id="blast-results-html-store"),
        dcc.Store(id="captain-results-store"),
        dmc.Grid(
            justify="start",
            align="start",
            style={"paddingTop": "20px"},
            gutter="xl",
            children=[
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 4,
                    },
                    children=[
                        html.H3(
                            "Search protein/nucleotide sequences for Starships and Starship-associated genes.",
                            style={
                                "textAlign": "center",
                                "fontSize": "1.25rem",
                            },
                        ),
                        dcc.Textarea(
                            id="query-text",
                            placeholder="Paste FASTA sequence here...",
                            rows=5,
                            style={
                                "width": "100%",
                                "fontSize": "1rem",
                                "padding": "10px",
                            },
                        ),
                        html.H3(
                            "Or",
                            style={
                                "textAlign": "center",
                                "fontSize": "1.25rem",
                            },
                        ),
                        dcc.Upload(
                            id="blast-fasta-upload",
                            children=html.Div(
                                "Drag and drop or click to select a FASTA file.",
                                id="blast-fasta-sequence-upload",
                                style={"fontSize": "1rem"},
                            ),
                            className="upload-box",
                            multiple=False,
                            accept=".fa, .fas, .fasta, .fna",
                        ),
                        html.Div(
                            id="upload-error-message",
                            style={"color": "red", "marginTop": "1rem"},
                        ),
                        curated_switch(
                            text="Only search curated Starships", size="normal"
                        ),
                        dbc.Button(
                            "Submit BLAST",
                            id="submit-button",
                            n_clicks=0,
                            className="d-grid gap-2 col-6 mx-auto",
                            style={"fontSize": "1rem"},
                        ),
                        dcc.Loading(
                            id="family-loading",
                            type="circle",
                            children=html.Div(
                                id="ship-family",
                                style={"paddingTop": "20px"},
                            ),
                        ),
                    ],
                ),
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        # dcc.Loading(
                        #     id="blast-chord-loading",
                        #     type="circle",
                        #     children=[html.Div(id="blast-chord")],
                        # ),
                        dcc.Loading(
                            id="ship-blasterjs-iframe-loading",
                            type="circle",
                            className="dash-loading",
                            children=html.Div(id="ship-blasterjs-iframe"),
                        ),
                        html.Button("Next Result", id="next-iframe-button"),
                        modal,
                        dcc.Loading(
                            id="subject-seq-button-loading",
                            type="circle",
                            children=html.Div(id="subject-seq-button"),
                        ),
                    ],
                ),
            ],
        ),
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
        logger.info(
            f"preprocess called with n_clicks={n_clicks}, query_text_input={query_text_input}, query_file_contents={query_file_contents}"
        )

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
        # Output("blast-results-store", "data"),
        Output("blast-results-html-store", "data"),
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


@callback(
    Output("ship-family", "children"),
    [
        Input("blast-results-store", "data"),
        Input("captain-results-store", "data"),
        Input("curated-input", "value"),
    ],
    State("submit-button", "n_clicks"),
)
def update_ui(
    blast_results_dict,
    captain_results_dict,
    curated,
    n_clicks,
):
    if blast_results_dict is None and captain_results_dict is None:
        raise PreventUpdate
    if n_clicks:
        logger.info(f"Updating UI with n_clicks={n_clicks}")
        no_captain_alert = dbc.Alert(
            "No captain sequence found (e-value threshold 0.01).",
            color="warning",
        )

        try:
            ship_family = no_update

            if blast_results_dict:
                try:
                    # TODO: caching the curated dataset makes no sense. filter the full dataset based on curated flag after loading from cache.
                    initial_df = load_from_cache("meta_data")

                    if initial_df is None:
                        initial_df = fetch_meta_data(curated)
                except Exception as e:
                    logger.error(f"Error fetching metadata: {str(e)}")

                logger.info("Parsing BLAST Results for Classification")
                blast_results_df = pd.DataFrame(blast_results_dict)

                df_for_table = pd.merge(
                    initial_df,
                    blast_results_df,
                    left_on="accession_tag",
                    right_on="sseqid",
                    how="right",
                )
                # we want to remove true duplicates, while keeping other hits which may be at a different location in the same ship
                df_for_table = df_for_table.drop_duplicates(
                    subset=["accession_tag", "pident", "length"]
                )
                df_for_table = df_for_table[df_for_table["accession_tag"].notna()]
                df_for_table.fillna("", inplace=True)

                min_evalue_rows = df_for_table.loc[
                    df_for_table.groupby("qseqid")["evalue"].idxmin()
                ]
                if min_evalue_rows.empty:
                    logger.warning(
                        "min_evalue_rows is empty after grouping by qseqid and selecting min evalue."
                    )
                else:
                    print(min_evalue_rows["pident"])
                    # logger.info(f"min_evalue_rows contents: {min_evalue_rows}")

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
                                            ship_family = dbc.Alert(
                                                [
                                                    f"Your sequence is likely in Starship family: {family} (Alignment length = {family_aln_length}, evalue = {family_evalue})",
                                                ],
                                                color="warning",
                                            )
                                        else:
                                            ship_family = dbc.Alert(
                                                [
                                                    f"Starship family could not be determined."
                                                ],
                                                color="danger",
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
            return ship_family

        except Exception as e:
            logger.error(f"Error in update_ui: {str(e)}")
            return no_update


@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("ship-blast-table", "data"), State("ship-blast-table", "columns")],
)
def download_tsv(n_clicks, rows, columns):
    if n_clicks == 0:
        return None

    if not rows or not columns:
        logger.error("Error: No data available for download.")
        return None

    try:
        df = pd.DataFrame(rows, columns=[c["name"] for c in columns])

        tsv_string = df.to_csv(sep="\t", index=False)
        tsv_bytes = io.BytesIO(tsv_string.encode())
        b64 = base64.b64encode(tsv_bytes.getvalue()).decode()

        today = date.today().strftime("%Y-%m-%d")

        return dict(
            content=f"data:text/tab-separated-values;base64,{b64}",
            filename=f"starbase_blast_{today}.tsv",
        )

    except Exception as e:
        logger.error(f"Error while preparing TSV download: {e}")
        return None


@callback(
    Output("lastz-plot", "figure"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
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


@callback(
    Output("blast-modal", "is_open"),
    Output("blast-modal-content", "children"),
    Output("blast-modal-title", "children"),
    Input("ship-blast-table", "active_cell"),
    State("blast-modal", "is_open"),
    State("ship-blast-table", "data"),
)
def toggle_modal(active_cell, is_open, table_data):
    if active_cell:
        row = active_cell["row"]
        row_data = table_data[row]
        modal_content, modal_title = create_accession_modal(row_data["accession_tag"])
        return True, modal_content, modal_title
    return is_open, None, None


@callback(
    Output(
        "ship-blasterjs-iframe", "children"
    ),  # Ensure "ship-blasterjs-iframe" is a container (e.g., html.Div)
    Input("next-iframe-button", "n_clicks"),
    State("blast-results-html-store", "data"),
    prevent_initial_call=True,
)
def update_iframe(n_clicks, blast_results_html_dict):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate

    if not blast_results_html_dict:
        output = dbc.Alert("No BLAST results found.", color="danger")
    else:
        # Convert dictionary keys to a list
        keys = list(blast_results_html_dict.keys())

        # Use n_clicks to find the correct key
        index = n_clicks % len(keys)
        selected_key = keys[index]
        file_path = blast_results_html_dict[selected_key]

        # Read the HTML file contents
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                html_content = file.read()
        except FileNotFoundError:
            return dbc.Alert(f"File not found: {file_path}", color="danger")

        # Create the Iframe with inline HTML content
        output = html.Iframe(srcDoc=html_content, width="100%", height="600px")

    return output

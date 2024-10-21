import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.ERROR)

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
)
from src.components.callbacks import curated_switch
from src.utils.parsing import parse_fasta, parse_fasta_from_file
from src.components.sqlite import engine
from src.utils.blastdb import db_list


dash.register_page(__name__)


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
                        dcc.Store(id="upload-error-store"),
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
                            id="ship-blast-table-loading",
                            type="circle",
                            className="dash-loading",
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
            logging.error(e)
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
        logging.info(
            f"preprocess called with n_clicks={n_clicks}, query_text_input={query_text_input}, query_file_contents={query_file_contents}"
        )

        input_type, query_header, query_seq = check_input(
            query_text_input, query_file_contents
        )
        logging.info(
            f"check_input returned input_type={input_type}, query_header={query_header}, query_seq={query_seq}"
        )

        if input_type in ("none", "both"):
            logging.info("Invalid input type; returning None.")
            return None, None, None

        query_type = guess_seq_type(query_seq)
        logging.info(f"guess_seq_type returned query_type={query_type}")

        return query_header, query_seq, query_type

    except Exception as e:
        logging.error(f"Error in preprocess: {str(e)}")
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
            logging.error("Missing query header or sequence.")
            return None, None, None

        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        logging.info(f"Temp FASTA written: {tmp_query_fasta}")

        # Run BLAST
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name

        try:
            blast_results = run_blast(
                db_list=db_list,
                query_type=query_type,
                query_fasta=tmp_query_fasta,
                tmp_blast=tmp_blast,
                input_eval=0.01,
                threads=2,
            )
            logging.info(f"BLAST results: {blast_results.head()}")
            if blast_results is None:
                raise ValueError("BLAST returned no results!")
        except Exception as e:
            logging.error(f"BLAST error: {str(e)}")
            raise

        blast_results_dict = blast_results.to_dict("records")

        # Run HMMER
        logging.info(f"Running HMMER")

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
                logging.error("Diamond/HMMER returned no results!")
                raise
        except Exception as e:
            logging.error(f"Diamond/HMMER error: {str(e)}")
            raise

        return blast_results_dict, results_dict, subject_seq_button

    except Exception as e:
        logging.error(f"Error in fetch_captain: {str(e)}")
        return None, None, None


@callback(
    Output("subject-seq-dl-package", "data"),
    [Input("subject-seq-button", "n_clicks"), Input("subject-seq", "data")],
)
def subject_seq_download(n_clicks, filename):
    try:
        if n_clicks:
            logging.info(f"Download initiated for file: {filename}")
            return dcc.send_file(filename)
        else:
            return dash.no_update
    except Exception as e:
        logging.error(f"Error in subject_seq_download: {str(e)}")
        return dash.no_update


no_captain_alert = dbc.Alert(
    "No captain sequence found (e-value threshold 0.01).",
    color="warning",
)


@callback(
    [
        Output("ship-family", "children"),
        Output("ship-blast-table", "children"),
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
    if n_clicks:
        logging.info(f"Updating UI with n_clicks={n_clicks}")
        try:
            ship_family = no_update
            ship_table = no_update

            blast_results_df = pd.DataFrame(blast_results_dict)

            query = """
            SELECT j.starshipID, a.accession_tag, f.familyName, t.species
            FROM accessions a
            LEFT JOIN joined_ships j ON a.id = j.ship_id
            LEFT JOIN taxonomy t ON j.taxid = t.id
            LEFT JOIN family_names f ON j.ship_family_id = f.id
            WHERE j.orphan IS NULL
            """
            # user can switch back and forth between hq and all ships in the output, without having to run a new search
            if curated:
                query += "AND j.curated_status == 'curated'"

            # accessions = blast_results_df["sseqid"]

            # if accessions is not None and len(accessions) > 1:
            #     placeholders = ",".join(["?"] * len(accessions))
            #     query += f"AND accession_tag IN ({placeholders})"
            #     initial_df = pd.read_sql_query(query, engine, params=accessions)
            # else:

            initial_df = pd.read_sql_query(query, engine)

            if blast_results_dict:
                logging.info("Rendering BLAST table")
                df_for_table = pd.merge(
                    initial_df,
                    blast_results_df,
                    left_on="accession_tag",
                    right_on="sseqid",
                    how="right",
                )

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
                if min_evalue_rows["pident"][0] > 95:
                    family_name = min_evalue_rows["familyName"][0]
                    aln_len = min_evalue_rows["length"][0]
                    ev = min_evalue_rows["evalue"][0]
                    ship_family = dbc.Alert(
                        [
                            f"Your sequence is likely in Starship family: {family_name} (Alignment length = {aln_len}, evalue = {ev})",
                        ],
                        color="warning",
                    )

                else:
                    if captain_results_dict:
                        logging.info("Processing Diamond/HMMER results")
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
                                logging.error(f"Error selecting ship family: {str(e)}")
                                ship_family = html.Div(f"Error: {str(e)}")
                        else:
                            ship_family = no_captain_alert
                    else:
                        ship_family = no_captain_alert

            return ship_family, ship_table

        except Exception as e:
            logging.error(f"Error in update_ui: {str(e)}")
            return no_update, no_update


@callback(
    Output("ship-aln", "children"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
        Input("curated-input", "value"),
    ],
    State("query-type-store", "data"),
)
def blast_alignments(ship_blast_results, selected_row, curated, query_type):
    try:
        logging.info(
            f"blast_alignments called selected_row: {selected_row}, query_type: {query_type}"
        )

        if not selected_row or len(selected_row) == 0:
            return [None]

        if not ship_blast_results or len(ship_blast_results) == 0:
            logging.error(
                "No BLAST results available because ship_blast_results is empty or None."
            )
            raise

        ship_blast_results_df = pd.DataFrame(ship_blast_results)

        row_idx = selected_row[0]

        try:
            row = ship_blast_results_df.iloc[row_idx]
            qseq = str(row["qseq"])
            qseqid = str(row["qseqid"])
            sseq = str(row["sseq"])
            sseqid = str(row["sseqid"])

        except IndexError:
            logging.error(f"Error: Row index {row_idx} out of bounds.")
            raise
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)

        try:
            with open(tmp_fasta.name, "w") as f:
                f.write(f">{qseqid}\n{qseq}\n>{sseqid}\n{sseq}\n")
        except Exception as file_error:
            logging.error(f"Error writing to FASTA file: {file_error}")
            raise

        try:
            with open(tmp_fasta.name, "r") as file:
                data = file.read()
        except Exception as read_error:
            logging.error(f"Error reading FASTA file: {read_error}")
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

        return [aln]

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        raise


@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("ship-blast-table", "data"), State("ship-blast-table", "columns")],
)
def download_tsv(n_clicks, rows, columns):
    if n_clicks == 0:
        return None

    if not rows or not columns:
        logging.error("Error: No data available for download.")
        return None

    try:
        # Convert rows and columns to DataFrame
        df = pd.DataFrame(rows, columns=[c["name"] for c in columns])

        # Generate TSV string
        tsv_string = df.to_csv(sep="\t", index=False)
        tsv_bytes = io.BytesIO(tsv_string.encode())
        b64 = base64.b64encode(tsv_bytes.getvalue()).decode()

        # Get today's date
        today = date.today().strftime("%Y-%m-%d")

        # Return the data
        return dict(
            content=f"data:text/tab-separated-values;base64,{b64}",
            filename=f"starbase_blast_{today}.tsv",
        )

    except Exception as e:
        logging.error(f"Error while preparing TSV download: {e}")
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
            logging.error("Error: No blast results available for plotting.")
            return None
    except Exception as e:
        logging.error(f"Error converting blast results to DataFrame: {e}")
        return None

    # Check if selected row is valid
    if selected_row is None or not selected_row:
        logging.error("No row selected for alignment.")
        return None

    try:
        # Extract sequence info from selected row
        row = ship_blast_results_df.iloc[selected_row]
        qseq = re.sub("-", "", row["qseq"])
        qseqid = row["qseqid"]
        sseq = re.sub("-", "", row["sseq"])
        sseqid = row["sseqid"]

        # Log selected sequence information
        logging.info(f"Selected query ID: {qseqid}, subject ID: {sseqid}")

        # Write sequences to temporary FASTA file
        with open(tmp_fasta_clean.name, "w") as f:
            f.write(f">{qseqid}\n{qseq}\n>{sseqid}\n{sseq}\n")

        logging.info("Running LASTZ...")
        run_lastz(tmp_fasta_clean.name, lastz_output.name)

        lastz_df = parse_lastz_output(lastz_output.name)

        if lastz_df.empty:
            logging.error("Error: No alignment data from LASTZ.")
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
        logging.error(f"Error: Selected row index is out of bounds. Details: {e}")
        return None

    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return None


@callback(
    Output("upload-error-message", "children"),
    Output("upload-error-store", "data"),
    Output("submit-button", "disabled"),
    Input("blast-fasta-upload", "contents"),
    State("blast-fasta-upload", "filename"),
    prevent_initial_call=True,
)
def handle_fasta_upload(contents, filename):
    if contents is None:
        return "", None
    max_size = 10 * 1024 * 1024  # 10 MB in bytes

    content_type, content_string = contents.split(",")

    header, fasta_length_error_message = parse_fasta_from_file(contents)

    decoded = base64.b64decode(content_string)
    file_size = len(decoded)

    if fasta_length_error_message:
        error_message = dbc.Alert(
            f"Error: {fasta_length_error_message}", color="danger"
        )
        return error_message, error_message, True
    elif file_size > max_size:
        error_message = dbc.Alert(
            f"Error: The file '{filename}' exceeds the 10 MB limit.", color="danger"
        )
        return error_message, error_message, True
    else:
        return "", None, False

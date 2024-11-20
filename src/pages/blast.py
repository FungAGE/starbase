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


from src.components.cache import cache
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
from src.components.callbacks import curated_switch, create_accession_modal
from src.utils.seq_utils import parse_fasta, parse_fasta_from_file
from src.components.sql_manager import load_from_cache
from src.components.sql_manager import fetch_meta_data
from src.utils.blastdb import db_list

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
        dcc.Store(id="captain-results-store"),
        
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
        logger.info(f"Updating UI with n_clicks={n_clicks}")
        try:
            ship_family = no_update
            ship_table = no_update

            blast_results_df = pd.DataFrame(blast_results_dict)

            # TODO: caching the curated dataset makes no sense. filter the full dataset based on curated flag after loading from cache.
            initial_df = load_from_cache("meta_data")

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
                    print(min_evalue_rows["pident"])
                    # logger.info(f"min_evalue_rows contents: {min_evalue_rows}")

                if not min_evalue_rows.empty and "pident" in min_evalue_rows.columns:
                    if min_evalue_rows["pident"].iloc[0] > 95:
                        family_name = min_evalue_rows["familyName"].iloc[0]
                        aln_len = min_evalue_rows["length"].iloc[0]
                        ev = min_evalue_rows["evalue"].iloc[0]
                        ship_family = dmc.Alert(
                            title="Starship Family Found",
                            children=[
                                f"Your sequence is likely in Starship family: {family_name}",
                                dmc.Space(h=5),
                                dmc.Text(
                                    f"Alignment length = {aln_len}, E-value = {ev}",
                                    size="sm",
                                    c="dimmed"
                                ),
                            ],
                            color="yellow",
                            variant="light",
                            withCloseButton=False,
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
                                            ship_family = dmc.Alert(
                                                title="Starship Family Found",
                                                children=[
                                                    f"Your sequence is likely in Starship family: {family}",
                                                    dmc.Space(h=5),
                                                    dmc.Text(
                                                        f"Alignment length = {family_aln_length}, E-value = {family_evalue}",
                                                        size="sm",
                                                        c="dimmed"
                                                    ),
                                                ],
                                                color="yellow",
                                                variant="light",
                                                withCloseButton=False,
                                            )
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
            return no_update, no_update


@cache.memoize()
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
        logger.info(
            f"blast_alignments called selected_row: {selected_row}, query_type: {query_type}"
        )

        if not selected_row or len(selected_row) == 0:
            return [None]

        if not ship_blast_results or len(ship_blast_results) == 0:
            logger.error(
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

        return [aln]

    except Exception as e:
        logger.error(f"Error: {str(e)}")
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
        Output("blast-modal", "is_open"),
        Output("blast-modal-content", "children"),
        Output("blast-modal-title", "children"),
        Output("ship-blast-table", "active_cell"),
    ],
    [Input("blast-table", "active_cell")],
    [
        State("blast-modal", "is_open"),
        State("ship-blast-table", "data"),
        State("ship-blast-table", "derived_virtual_data")  # Add this
    ],
)
def toggle_modal(active_cell, is_open, table_data, filtered_data):
    if active_cell:
        # Use filtered data if available
        data_to_use = filtered_data if filtered_data is not None else table_data
        row_data = data_to_use[active_cell["row"]]
        modal_content, modal_title = create_accession_modal(row_data["accession_tag"])
        return True, modal_content, modal_title, None
    return is_open, no_update, no_update, no_update
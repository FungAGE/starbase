import warnings

warnings.filterwarnings("ignore")

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
    blast_table,
    run_lastz,
    select_ship_family,
    parse_lastz_output,
    blast_chords,
)
from src.utils.tree import plot_tree, default_highlight_families

dash.register_page(__name__)

db_list = {
    "ship": {"nucl": "database_folder/Starships/ships/fna/blastdb/concatenated.fa"},
    "gene": {
        "tyr": {
            "nucl": "database_folder/Starships/captain/tyr/fna/blastdb/concatenated.fa",
            "prot": "database_folder/Starships/captain/tyr/faa/blastdb/concatenated.faa",
            "hmm": {
                "nucl": "database_folder/Starships/captain/tyr/fna/hmm/YRsuperfamRefs.mafft.hmm",
                "prot": "database_folder/Starships/captain/tyr/faa/hmm/YRsuperfams.p1-512.hmm",
            },
        },
        "nlr": {
            "nucl": "database_folder/Starships/cargo/nlr/fna/blastdb/nlr.fa",
            "prot": "database_folder/Starships/cargo/nlr/faa/blastdb/nlr.mycoDB.faa",
        },
        "fre": {
            "nucl": "database_folder/Starships/cargo/fre/fna/blastdb/fre.fa",
            "prot": "database_folder/Starships/cargo/fre/faa/blastdb/fre.mycoDB.faa",
        },
        "plp": {
            "nucl": "database_folder/Starships/cargo/plp/fna/blastdb/plp.fa",
            "prot": "database_folder/Starships/cargo/plp/faa/blastdb/plp.mycoDB.faa",
        },
        "duf3723": {
            "nucl": "database_folder/Starships/cargo/duf3723/fna/blastdb/duf3723.fa",
            "prot": "database_folder/Starships/cargo/duf3723/faa/blastdb/duf3723.mycoDB.faa",
        },
    },
}


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
        dcc.Store(id="hmmer-results-store"),
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
                            id="fasta-upload",
                            children=html.Div(
                                "Drag and drop or click to select a FASTA file.",
                                id="fasta-sequence-upload",
                                style={"fontSize": "1rem"},
                            ),
                            className="upload-box",
                            multiple=False,
                            accept=".fa, .fas, .fasta, .fna",
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
                            type="default",
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
                        #     type="default",
                        #     children=[html.Div(id="blast-chord")],
                        # ),
                        dcc.Loading(
                            id="ship-blast-table-loading",
                            type="default",
                            children=html.Div(id="ship-blast-table"),
                        ),
                        dcc.Loading(
                            id="subject-seq-button-loading",
                            type="default",
                            children=html.Div(id="subject-seq-button"),
                        ),
                        dcc.Loading(
                            id="ship-aln-loading",
                            type="default",
                            children=html.Div(id="ship-aln"),
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@callback(
    [
        Output("query-header-store", "data"),
        Output("query-seq-store", "data"),
        Output("query-type-store", "data"),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("query-text", "value"),
        Input("fasta-upload", "contents"),
    ],
)
def preprocess(n_clicks, query_text_input, query_file_contents):
    if not n_clicks:
        raise PreventUpdate

    input_type, query_header, query_seq = check_input(
        query_text_input, query_file_contents
    )

    if input_type in ("none", "both"):
        return None, None, None

    query_type = guess_seq_type(query_seq)

    return query_header, query_seq, query_type


@callback(
    [
        Output("blast-results-store", "data"),
        Output("hmmer-results-store", "data"),
        Output("subject-seq-button", "children"),
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
    ],
)
def fetch_blast_hmmer_results(query_header, query_seq, query_type):
    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name
        tmp_hmmer = tempfile.NamedTemporaryFile(suffix=".hmmer.txt").name
        tmp_hmmer_parsed = tempfile.NamedTemporaryFile(suffix=".hmmer.parsed.txt").name

        # Run BLAST
        blast_results = run_blast(
            db_list=db_list,
            query_type=query_type,
            query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
        )
        if blast_results is None:
            raise ValueError("BLAST returned no results!")
        blast_results_dict = blast_results.to_dict("records")

        # TODO: create grouped hmm profile for nucl captains so that hit_ID returned is a captain family
        subject_seq_button = None
        subject_seq = None
        # Run HMMER
        hmmer_results, subject_seq = run_hmmer(
            db_list=db_list,
            query_type=query_type,
            input_genes="tyr",
            input_eval=0.01,
            query_fasta=tmp_query_fasta,
            tmp_hmmer=tmp_hmmer,
            tmp_hmmer_parsed=tmp_hmmer_parsed,
            threads=2,
        )

        if hmmer_results is None:
            raise ValueError("hmmsearch returned no results!")

        hmmer_results_dict = hmmer_results.to_dict("records")

        return blast_results_dict, hmmer_results_dict, subject_seq_button

    except Exception as e:
        print(f"Error: {str(e)}")
        return None, None, None


@callback(
    Output("subject-seq-dl-package", "data"),
    [Input("subject-seq-button", "n_clicks"), Input("subject-seq", "data")],
)
def subject_seq_download(n_clicks, filename):
    if n_clicks:
        return dcc.send_file(filename)
    else:
        return dash.no_update


@callback(
    [
        Output("ship-family", "children"),
        Output("ship-blast-table", "children"),
        # Output("blast-chord", "children"),
    ],
    [
        Input("blast-results-store", "data"),
        Input("hmmer-results-store", "data"),
    ],
    [
        State("submit-button", "n_clicks"),
        State("joined-ships", "data"),
        State("query-type-store", "data"),
    ],
)
def update_ui(
    blast_results_dict, hmmer_results_dict, n_clicks, cached_data, query_type
):
    ship_family = no_update
    ship_table = no_update
    # blast_chord = no_update

    if blast_results_dict is None and hmmer_results_dict is None:
        raise PreventUpdate

    if n_clicks:
        initial_df = pd.DataFrame(cached_data)
        if blast_results_dict:
            # Render BLAST table
            blast_results_df = pd.DataFrame(blast_results_dict)
            ship_table = blast_table(blast_results_df)
            # blast_chord = blast_chords(blast_results_df)

        # Process HMMER results and render family
        # if query_type == "nucl":
        #     ship_family = dbc.Alert(
        #         "hmmsearch will currently not be run for nucleotide queries",
        #         color="warning",
        #     )
        # else:
        if hmmer_results_dict:
            hmmer_results_df = pd.DataFrame(hmmer_results_dict)
            try:
                superfamily, family_aln_length, family_evalue = select_ship_family(
                    hmmer_results_df
                )
                if superfamily:
                    family = initial_df[initial_df["longFamilyID"] == superfamily][
                        "familyName"
                    ].unique()[0]
                    ship_family = dbc.Alert(
                        [
                            f"Your sequence is likely in Starship family: {family} (Alignment length = {family_aln_length}, evalue = {family_evalue})",
                        ],
                        color="warning",
                    )
            except Exception as e:
                ship_family = html.Div(f"Error: {str(e)}")
        else:
            ship_family = dbc.Alert(
                "No captain sequence found (e-value threshold 0.01).",
                color="warning",
            )

        return (
            ship_family,
            ship_table,
        )  # blast_chord


@callback(
    Output("ship-aln", "children"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
    ],
    State("query-type-store", "data"),
)
def blast_alignments(ship_blast_results, selected_row, query_type):
    try:
        # Ensure selected_row is not empty
        if not selected_row or len(selected_row) == 0:
            return [None]

        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)
        ship_blast_results_df = pd.DataFrame(ship_blast_results)

        # Get the first selected row index
        row_idx = selected_row[0]

        try:
            row = ship_blast_results_df.iloc[row_idx]
        except IndexError:
            return ["Error: Selected row index out of bounds"]

        # Write sequences to temp fasta file
        with open(tmp_fasta.name, "w") as f:
            f.write(f">{row['qseqid']}\n")
            f.write(f"{row['qseq']}\n")
            f.write(f">{row['sseqid']}\n")
            f.write(f"{row['sseq']}\n")

        with open(tmp_fasta.name, "r") as file:
            data = file.read()

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
        return [f"Error: {str(e)}"]


# Define callback to download TSV
@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("ship-blast-table", "data"), State("ship-blast-table", "columns")],
)
def download_tsv(n_clicks, rows, columns):
    if n_clicks == 0:
        return None

    df = pd.DataFrame(rows, columns=[c["name"] for c in columns])
    tsv_string = df.to_csv(sep="\t", index=False)
    tsv_bytes = io.BytesIO(tsv_string.encode())
    b64 = base64.b64encode(tsv_bytes.getvalue()).decode()

    today = date.today().strftime("%Y-%m-%d")

    return dict(
        content=f"data:text/tab-separated-values;base64,{b64}",
        filename=f"starbase_blast_{today}.tsv",
    )


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

    # Convert ship_blast_results to a Pandas DataFrame
    ship_blast_results_df = pd.DataFrame(ship_blast_results)

    if selected_row is None:
        return None
    else:
        row = ship_blast_results_df.iloc[selected_row]
        qseq = re.sub("-", "", row["qseq"].iloc[0])
        qseqid = row["qseqid"].iloc[0]
        sseq = re.sub("-", "", row["sseq"].iloc[0])
        sseqid = row["sseqid"].iloc[0]

        # Write the specific row to the file
        with open(tmp_fasta_clean.name, "w") as f:
            f.write(f">{qseqid}\n")
            f.write(f"{qseq}\n")
            f.write(f">{sseqid}\n")
            f.write(f"{sseq}\n")

        # Run LASTZ alignment
        run_lastz(tmp_fasta_clean.name, lastz_output.name)

        # Parse LASTZ output
        lastz_df = parse_lastz_output(lastz_output.name)

        # Extract individual columns and prepare the data for plotting
        x_values = []
        y_values = []

        for _, row in lastz_df.iterrows():
            x_values.append(row["qstart"])
            x_values.append(row["qend"])
            y_values.append(row["sstart"])
            y_values.append(row["send"])

        # Create the scatter plot using Plotly Graph Objects
        fig = go.Figure(data=go.Scatter(x=x_values, y=y_values, mode="lines"))

        # Add layout details
        fig.update_layout(
            title="LASTZ Alignment",
            xaxis_title=qseqid,
            yaxis_title=sseqid,
        )

        return fig

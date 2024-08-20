import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
from dash import dcc, html, callback, MATCH
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
    gene_hmmsearch,
    parse_lastz_output,
)
from src.utils.tree import plot_tree, tree_file, metadata, default_highlight_clades

dash.register_page(__name__)

db_list = {
    "ship": {"nucl": "database_folder/Starships/ships/fna/blastdb/concatenated.fa"},
    "gene": {
        "tyr": {
            "nucl": "database_folder/Starships/captain/tyr/fna/blastdb/concatenated.fa",
            "prot": "database_folder/Starships/captain/tyr/faa/blastdb/concatenated.faa",
            "hmm": "database_folder/Starships/captain/tyr/hmm/YRsuperfams.p1-512.hmm",
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

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            align="top",
            children=[
                dbc.Col(
                    style={"padding": "20px"},
                    sm=12,
                    lg=4,
                    children=[
                        dbc.Stack(
                            [
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
                            ],
                            gap=4,
                            direction="vertical",
                        )
                    ],
                ),
                dbc.Col(
                    sm=12,
                    lg=8,
                    style={"padding": "20px"},
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="default",
                            children=[
                                html.Div(id="output-container"),
                                html.Div(id="output-container-error-message"),
                            ],
                        ),
                        dcc.Loading(
                            id="loading-2",
                            type="default",
                            children=[
                                html.Div(id="ship-aln-container"),
                                html.Div(id="ship-aln-container-error-message"),
                            ],
                        ),
                        # dcc.Loading(
                        #     id="loading-3",
                        #     type="default",
                        #     children=html.Div(dcc.Graph(id="lastz-plot")),
                        # ),
                    ],
                ),
            ],
        ),
    ],
)


@callback(
    [
        Output("output-container", "children"),
        Output("output-container-error-message", "children"),
    ],
    [
        Input("submit-button", "n_clicks"),
        Input("query-text", "value"),
        Input("fasta-upload", "contents"),
        # Input("fasta-upload", "filename"),
    ],
)
def run_main(n_clicks, query_text_input, query_file_contents):
    global query_type
    tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name
    if not n_clicks:
        raise PreventUpdate

    try:
        input_type, query_header, query_seq = check_input(
            query_text_input, query_file_contents
        )
        if input_type in ("none", "both"):
            return [
                dbc.Card(
                    [
                        dbc.CardHeader(html.H4("Invalid input:")),
                        dbc.CardBody(
                            html.H4(
                                "Please provide either a query sequence in the text box, or upload a FASTA file"
                            )
                        ),
                    ],
                    color="red",
                )
            ]

        query_type = guess_seq_type(query_seq)
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)

        blast_results = run_blast(
            db_list=db_list,
            query_type=query_type,
            tmp_query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
        )
        hmmer_results = run_hmmer(
            db_list=db_list,
            query_type=query_type,
            input_genes="tyr",
            input_eval=0.01,
            query_fasta=tmp_query_fasta,
            threads=2,
        )

        ship_blast_table = blast_table(blast_results)
        # chords = blast_chords(blast_results)

        if hmmer_results is not None:
            superfamily = gene_hmmsearch(hmmer_results)
            if superfamily is not None:
                superfamily_text = html.Div(
                    [
                        html.H4(
                            f"Likely captain family: {superfamily}",
                        ),
                    ],
                )
                # superfamily_tree = dcc.Graph(
                #     id="blast-phylogeny",
                #     className="div-card",
                #     figure=plot_tree(tree_file, metadata, highlight_clades=superfamily),
                # )

        else:
            superfamily_text = """"""
            # superfamily_tree = """"""

        return (
            html.Div(
                [
                    dbc.Stack(
                        [superfamily_text, ship_blast_table],
                        gap=3,
                    )
                ],
                id="ship-blast-table-container",
            ),
            "",
        )

    except Exception as e:
        return None, f"Error: {str(e)}"


# Callback to update information about selected row
@callback(
    [
        Output("ship-aln-container", "children"),
        Output("ship-aln-container-error-message", "children"),
    ],
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
    ],
)
def blast_alignments(ship_blast_results, selected_row):
    try:
        global query_type

        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)

        ship_blast_results_df = pd.DataFrame(ship_blast_results)

        if selected_row:
            try:
                row = ship_blast_results_df.iloc[selected_row]
            except IndexError:
                row = None
            if row is not None:
                with open(tmp_fasta.name, "w") as f:
                    f.write(f">{row['qseqid'].iloc[0]}\n")
                    f.write(f"{row['qseq'].iloc[0]}\n")
                    f.write(f">{row['sseqid'].iloc[0]}\n")
                    f.write(f"{row['sseq'].iloc[0]}\n")

                with open(tmp_fasta.name, "r") as file:
                    data = file.read()
                if query_type == "nucl":
                    color = "nucleotide"
                else:
                    color = "clustal2"

                aln = dashbio.AlignmentChart(
                    id="alignment-viewer",
                    data=data,
                    height=200,
                    tilewidth=30,
                    colorscale=color,
                    # overview="slider",
                    showconsensus=False,
                    showconservation=False,
                    showgap=False,
                    showid=False,
                    ticksteps=5,
                    tickstart=0,
                )

                return aln
    except Exception as e:
        return None, f"Error: {str(e)}"


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
    Output("ship-blast-table", "selected_rows"),
    Input("ship-blast-table", "data"),
)
def update_selected_rows(data):
    if data:
        return [0]
    return []


# # Callback to open the modal
# @callback(
#     Output("modal", "is_open"),
#     Input("open-modal-link", "n_clicks"),
#     State("modal", "is_open"),
# )
# def toggle_modal(n_clicks, is_open):
#     if n_clicks is not None:
#         if n_clicks:
#             return not is_open
#         return is_open


# # Callback to close the modal
# @callback(
#     Output("modal", "is_open"),
#     Input("close-modal", "n_clicks"),
#     State("modal", "is_open"),
# )
# def close_modal(n_clicks, is_open):
#     if n_clicks:
#         return False
#     return is_open


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

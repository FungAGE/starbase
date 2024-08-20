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
    select_ship_family,
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
        dcc.Store(id="query-header-store"),
        dcc.Store(id="query-seq-store"),
        dcc.Store(id="query-type-store"),
        dcc.Store(id="blast-results-store"),
        dcc.Store(id="hmmer-results-store"),
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
                        dbc.Row(
                            [
                                dcc.Loading(
                                    id="family-loading",
                                    type="default",
                                    children=html.Div(id="ship-family"),
                                ),
                                dcc.Loading(
                                    id="tree-loading",
                                    type="default",
                                    children=html.Div(id="blast-phylogeny"),
                                ),
                                dcc.Loading(
                                    id="ship-blast-table-loading",
                                    type="default",
                                    children=html.Div(id="ship-blast-table"),
                                ),
                                dcc.Loading(
                                    id="ship-aln-loading",
                                    type="default",
                                    children=html.Div(id="ship-aln"),
                                ),
                                # dcc.Loading(
                                #     id="loading-4",
                                #     type="default",
                                #     children=html.Div(dcc.Graph(id="lastz-plot")),
                                # ),
                            ]
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
        State("query-text", "value"),
        State("fasta-upload", "contents"),
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
    ],
    [
        Input("query-header-store", "data"),
        Input("query-seq-store", "data"),
        Input("query-type-store", "data"),
    ],
)
def fetch_blast_hmmer_results(query_header, query_seq, query_type):
    if not query_header or not query_seq or not query_type:
        raise PreventUpdate

    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name

        # Run BLAST
        blast_results = run_blast(
            db_list=db_list,
            query_type=query_type,
            tmp_query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
        )
        if blast_results is None:
            raise ValueError("BLAST returned no results!")
        blast_results_dict = blast_results.to_dict("records")

        # Run HMMER
        hmmer_results = run_hmmer(
            db_list=db_list,
            query_type=query_type,
            input_genes="tyr",
            input_eval=0.01,
            query_fasta=tmp_query_fasta,
            threads=2,
        )

        if hmmer_results is not None:
            # Convert HMMER results DataFrame to dict
            hmmer_results_dict = hmmer_results.to_dict("records")
        else:
            hmmer_results_dict = None

        return blast_results_dict, hmmer_results_dict

    except Exception as e:
        return None, None


@callback(
    [
        Output("ship-family", "children"),
        Output("blast-phylogeny", "children"),
        Output("ship-blast-table", "children"),
    ],
    [
        Input("blast-results-store", "data"),
        Input("hmmer-results-store", "data"),
    ],
)
def update_ui(blast_results_dict, hmmer_results_dict):
    if blast_results_dict:
        # Render BLAST table
        blast_results_df = pd.DataFrame(blast_results_dict)
        ship_table = blast_table(blast_results_df)

        # Process HMMER results and render family and tree
        ship_family = ""
        ship_tree = ""

        if hmmer_results_dict:
            hmmer_results_df = pd.DataFrame(hmmer_results_dict)
            try:
                superfamily, family_aln_length, family_evalue = select_ship_family(
                    hmmer_results_df
                )
                if superfamily:
                    ship_family = html.Div(
                        [
                            html.H4(
                                [
                                    f"Your sequence is likely in Starship family: {superfamily}",
                                    f"(Alignment length ={family_aln_length}, evalue = {family_evalue})",
                                ]
                            )
                        ]
                    )
                    ship_tree = dcc.Graph(
                        # id="blast-phylogeny",
                        className="div-card",
                        figure=plot_tree(
                            tree_file, metadata, highlight_clades=superfamily
                        ),
                    )
            except Exception as e:
                ship_family = html.Div(f"Error: {str(e)}")
        else:
            ship_family = dbc.Alert("No captain sequence found.", color="warning")

        return ship_family, ship_tree, ship_table


@callback(
    Output("ship-blast-table", "selected_rows"),
    Input("ship-blast-table", "derived_virtual_data"),
    Input("ship-blast-table", "derived_virtual_selected_rows"),
)
def update_selected_rows(data, selected_rows):
    if not selected_rows or len(selected_rows) == 0:
        return [0]
    return selected_rows


@callback(
    [
        Output("ship-aln", "children"),
    ],
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
        Input("query-type-store", "data"),
    ],
)
def blast_alignments(ship_blast_results, selected_row, query_type):
    if ship_blast_results and selected_row and query_type:
        try:
            # Ensure selected_row is not empty
            if not selected_row or len(selected_row) == 0:
                return None

            tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)
            ship_blast_results_df = pd.DataFrame(ship_blast_results)

            # Get the first selected row index
            row_idx = selected_row[0]

            try:
                row = ship_blast_results_df.iloc[row_idx]
            except IndexError:
                return "Error: Selected row index out of bounds"

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
            return f"Error: {str(e)}"


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

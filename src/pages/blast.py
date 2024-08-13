import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.exceptions import PreventUpdate
from dash.dependencies import Output, Input, State
import dash_bio as dashbio

import io
from io import StringIO
import os
import re
import tempfile
import base64
from datetime import date
import subprocess
import json
import pandas as pd

from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbitblastnCommandline,
    NcbiblastpCommandline,
)
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import plotly.graph_objects as go
from src.components.tree import plot_tree

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

tree_file = "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile"
metadata = pd.read_csv("src/data/superfam-clades.tsv", sep="\t")

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            align="top",
            children=[
                dbc.Col(
                    style={"padding": "20px"},
                    sm=10,
                    lg=4,
                    children=[
                        dbc.Stack(
                            [
                                html.H3(
                                    [
                                        "Search protein/nucleotide sequences for Starships and Starship-associated genes."
                                    ],
                                    style={"textAlign": "center"},
                                ),
                                dcc.Textarea(
                                    id="query-text",
                                    placeholder="Paste FASTA sequence here...",
                                    rows=5,
                                    style={
                                        "width": "100%",
                                    },
                                ),
                                html.H3(
                                    ["Or"],
                                    style={"textAlign": "center"},
                                ),
                                dcc.Upload(
                                    id="fasta-upload",
                                    children=html.Div(id="fasta-sequence-upload"),
                                    style={
                                        "width": "100%",
                                        "height": "75px",
                                        "lineHeight": "75px",
                                        "borderWidth": "2px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "justify-content": "center",
                                        "textAlign": "center",
                                    },
                                    multiple=False,
                                    accept=".fa, .fas, .fasta, .fna",
                                ),
                                dbc.Button(
                                    "Submit BLAST",
                                    id="submit-button",
                                    n_clicks=0,
                                    className="d-grid gap-2 col-6 mx-auto",
                                ),
                            ],
                            gap=3,
                            direction="vertical",
                        )
                    ],
                ),
                dbc.Col(
                    sm=10,
                    lg=8,
                    style={"padding": "20px"},
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="default",
                            children=html.Div(id="output-container"),
                        ),
                        dcc.Loading(
                            id="loading-2",
                            type="default",
                            children=html.Div(id="ship-aln-container"),
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

tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name


@callback(
    Output("output-container", "children"),
    [
        Input("submit-button", "n_clicks"),
        Input("query-text", "value"),
        Input("fasta-upload", "contents"),
        # Input("fasta-upload", "filename"),
    ],
)
def run_main(n_clicks, query_text_input, query_file_contents):
    global query_type
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
            tmp_query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
        )
        hmmer_results = run_hmmer(
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
                            f"Likely captain superfamily: {superfamily}",
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

        return html.Div(
            [
                dbc.Stack(
                    [superfamily_text, ship_blast_table],
                    gap=3,
                )
            ],
            id="ship-blast-table-container",
        )

    except Exception as e:
        # Log the error and handle it appropriately
        print(f"Error occurred: {e}")
        raise PreventUpdate


def write_temp_fasta(header, sequence):
    cleaned_query_seq = SeqRecord(Seq(sequence), id=header, description="")
    tmp_query_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    SeqIO.write(cleaned_query_seq, tmp_query_fasta, "fasta")
    return tmp_query_fasta


def check_input(query_text_input, query_file_contents):
    if query_text_input in ("", None) and query_file_contents is None:
        raise ValueError("No file contents provided.")

    elif query_text_input and query_file_contents:
        return "both", None, None

    elif query_text_input:
        input_type = "text"
        header, query = parse_fasta_from_text(query_text_input)

    elif query_file_contents:
        input_type = "file"
        header, query = parse_fasta_from_file(query_file_contents)

    return input_type, header, query


def parse_fasta_from_text(text):
    try:
        queries = SeqIO.parse(StringIO(text), "fasta")
        query = next(queries)
        header, seq = query.id, query.seq
        return header, seq
    except Exception as e:
        print(f"Error parsing text: {e}")
        return None, None


def parse_fasta_from_file(file_contents):
    try:
        if not file_contents:
            raise ValueError("No file contents provided.")
        split_contents = file_contents.split(",")

        # Extract the header from the first line
        header = split_contents[0].strip()  # Remove leading/trailing whitespace

        # Join the remaining lines to form the sequence
        sequence = "".join(split_contents[1:])

        # Decode the base64 encoded sequence
        decoded_sequence = base64.b64decode(sequence).decode("utf-8")

        # Parse the decoded sequence as FASTA
        query = SeqIO.read(StringIO(decoded_sequence), "fasta")

        # Extract header and sequence from the single record
        header, seq = query.id, str(query.seq)

        return header, seq

    except Exception as e:
        print(f"Error parsing file: {e}")
        return None, None

        # nseq = 0
        # for query in queries:
        #     nseq += 1
        #     # Stop if there are more than one sequence records
        #     if nseq > 1:
        #         raise ValueError(
        #             "More than one sequence record found in the FASTA file."
        #         )

        #     # Extract header and sequence from the single record
        #     header, seq = query.id, str(query.seq)
        #     return header, seq

        # # Handle the case where no sequence record was found
        # raise ValueError("No sequence record found in the FASTA file.")


def guess_seq_type(query_seq):
    nucl_char = set("ATGC")
    prot_char = set("ARNDBCEQZGHILKMFPSTWYV")

    if query_seq is not None:
        nucl_count = sum(1 for nt in query_seq.upper() if nt in nucl_char)
        prot_count = sum(1 for aa in query_seq.upper() if aa in prot_char)
        if nucl_count == prot_count:
            query_guess = "nucl"
            print("Query is nucleotide sequence")
        else:
            query_guess = "prot"
            print("Query is protein sequence")

        return query_guess
    else:
        # Handle the case where query_seq is None
        print("Error: query_seq is None")
        return None


def run_blast(
    tmp_query_fasta=None,
    tmp_blast=None,
    input_eval=None,
    threads=None,
):
    # TODO: add another set of dirs for hq Starships and all Starships?
    # ? instead of creating an additional set of blastdbs, why not just filter by quality in the results
    # * that way, the user can switch back and forth between hq and all ships in the output, without having to run a new search

    # Determine blast program and database based on sequence type and blast type
    blastdb = db_list["ship"]["nucl"]

    if query_type == "nucl":
        blast_program = NcbiblastnCommandline
    else:
        blast_program = NcbitblastnCommandline

    if os.path.exists(blastdb) and os.path.getsize(blastdb) > 0:
        print("Performing BLAST search...")
        # Perform the BLAST search
        blast_cline = blast_program(
            query=tmp_query_fasta,
            db=blastdb,
            evalue=input_eval,
            out=tmp_blast,
            outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq",
            num_threads=threads,
        )
        stdout, stderr = blast_cline()
    else:
        raise ValueError("Issue accessing BLAST database")

    # Optionally stitch BLAST results
    # if stitch:
    #     stitched_blast_tmp = tempfile.NamedTemporaryFile(suffix=".stitch").name
    #     stitch_blast_cmd = f"python bin/BLASTstitcher.py -i {tmp_blast} -o {stitched_blast_tmp}"
    #     subprocess.run(stitch_blast_cmd, shell=True)
    #     ship_blast_out = stitched_blast_tmp

    # with open(tmp_query_fasta, "r") as file:
    #     file_contents = file.read()
    # print(file_contents)

    # with open(tmp_blast, "r") as file:
    #     file_contents = file.read()
    # print(file_contents)

    df = pd.read_csv(
        tmp_blast,
        sep="\t",
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qseq",
            "sseq",
        ],
    )
    df["qseqid"] = df["qseqid"].apply(lambda x: re.sub(r"\|.*$", "", x))
    df["sseqid"] = df["sseqid"].apply(lambda x: re.sub(r"\|.*$", "", x))

    return df


def run_hmmer(input_genes="tyr", input_eval=None, query_fasta=None, threads=None):
    hmmer_db = None
    if input_genes is not None:
        if query_type == "prot":
            hmmer_program = "hmmsearch"
            hmmer_db = db_list["gene"][input_genes]["hmm"]

            if os.path.exists(hmmer_db) and os.path.getsize(hmmer_db) > 0:
                # Run HMMER search
                tmp_hmmer = tempfile.NamedTemporaryFile(suffix=".hmmer.txt").name
                hmmer_cmd = f"{hmmer_program} -o {tmp_hmmer} --cpu {threads} --domE {input_eval} {hmmer_db} {query_fasta}"
                print("Running hmmsearch...")
                subprocess.run(hmmer_cmd, shell=True)

                # Parse HMMER output
                tmp_hmmer_parsed = tempfile.NamedTemporaryFile(
                    suffix=".hmmer.parsed.txt"
                ).name
                n_records = parse_hmmer(tmp_hmmer, tmp_hmmer_parsed)
                print(f"Number of hmmsearch records: {n_records}")

                # extract sequence from results
                # extract_hmmer(tmp_hmmer_parsed)

                # Read parsed output into DataFrame
                hmmer_results = pd.read_csv(
                    tmp_hmmer_parsed,
                    sep="\t",
                    names=[
                        "query_id",
                        "hit_IDs",
                        "aln_length",
                        "query_start",
                        "query_end",
                        "gaps",
                        "query_seq",
                        "subject_seq",
                        "evalue",
                        "bitscore",
                    ],
                )
                return hmmer_results
            else:
                raise ValueError("Issue accessing HMM database")
    if hmmer_db is None:
        return None


# Parse the HMMER results
def parse_hmmer(hmmer_output_file, parsed_file):
    print("Parsing hmmer output...")
    with open(parsed_file, "w") as tsv_file:
        tsv_file.write(
            "query_id\thit_IDs\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )
        n_records = 0
        for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
            n_records = +1
            for hit in record.hits:
                for hsp in hit.hsps:
                    query_seq = str(hsp.query.seq)
                    subject_seq = str(hsp.hit.seq)
                    aln_length = hsp.aln_span
                    query_start = hsp.query_start
                    query_end = hsp.query_end
                    gaps = str("N/A")
                    bitscore = hsp.bitscore
                    evalue = hsp.evalue
                    tsv_file.write(
                        f"{hit.id}\t{record.id}\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t{query_seq}\t{subject_seq}\t{evalue}\t{bitscore}\n"
                    )
        return n_records


def extract_hmmer(parsed_file):
    # Read the TSV file into a DataFrame
    data = pd.read_csv(parsed_file, sep="\t")

    # Get rows with the lowest e-value for each unique entry in Query
    min_evalue_rows = data.loc[data.groupby("query_id")["evalue"].idxmin()]

    # Use os.path.join to construct file paths correctly
    top_hit_out_path = os.path.join(
        os.path.dirname(parsed_file), f"{os.path.splitext(parsed_file)[0]}.besthit.txt"
    )

    with open(top_hit_out_path, "w") as top_hit_out:
        # Write the header line
        top_hit_out.write(
            "hit_IDs\tquery_id\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
        )

        # Write the rows to the file using to_csv
        min_evalue_rows.to_csv(top_hit_out, sep="\t", header=False, index=False)

        for index, row in min_evalue_rows.iterrows():
            # Create a SeqRecord
            query = row["query_id"]
            qseq = row["query_seq"]
            sequence = SeqRecord(Seq(qseq), id=query, description="")

            # Write the SeqRecord to a FASTA file
            # Use os.path.join for constructing output file paths
            output_filename = os.path.join(
                os.path.dirname(parsed_file), f"{query}_best_hsp.fa"
            )

            SeqIO.write(sequence, output_filename, "fasta")


def circos_prep(blast_output, links_output, layout_output):
    best_blast = blast_output.nsmallest(10, "evalue")
    links = [
        {
            "source": {
                "id": row["qseqid"],
                "start": int(row["qstart"]),
                "end": int(row["qend"]),
            },
            "target": {
                "id": row["sseqid"],
                "start": int(row["sstart"]),
                "end": int(row["send"]),
            },
            # Use pident as the single numeric value for simplicity
            "value": row["pident"],
        }
        for index, row in best_blast.iterrows()
    ]
    # Save to JSON file
    with open(links_output, "w") as json_file:
        json.dump(links, json_file, indent=4)

    # Determine unique sequence IDs and their lengths
    sequence_lengths = {}
    for index, row in best_blast.iterrows():
        if row["qseqid"] not in sequence_lengths:
            sequence_lengths[row["qseqid"]] = int(row["qend"])
        else:
            sequence_lengths[row["qseqid"]] = max(
                sequence_lengths[row["qseqid"]], int(row["qend"])
            )

        if row["sseqid"] not in sequence_lengths:
            sequence_lengths[row["sseqid"]] = int(row["send"])
        else:
            sequence_lengths[row["sseqid"]] = max(
                sequence_lengths[row["sseqid"]], int(row["send"])
            )

    # Define colors for each chromosome block
    colors = [
        "#FF5733",
        "#33FF57",
        "#3357FF",
        "#FF33A6",
        "#33FFF1",
        "#F1FF33",
        "#B633FF",
        "#33FFB6",
        "#FFA633",
        "#33A6FF",
    ]

    # Create Circos layout with colors
    layout = [
        {
            "id": seq_id,
            "len": length,
            "color": colors[i % len(colors)],
            "label": {"display": seq_id, "color": "#000000", "size": 14},
        }
        for i, (seq_id, length) in enumerate(sequence_lengths.items())
    ]

    # Save layout to JSON file
    with open(layout_output, "w") as json_file:
        json.dump(layout, json_file, indent=4)


def blast_chords(blast_output):
    if blast_output is not None and not blast_output.empty:
        tmp_links_json = tempfile.NamedTemporaryFile(suffix=".json").name
        tmp_layout_json = tempfile.NamedTemporaryFile(suffix=".json").name

        # Prepare data for Circos plot
        circos_prep(blast_output, tmp_links_json, tmp_layout_json)

        # Load the prepared JSON data
        try:
            with open(tmp_links_json) as f:
                circos_graph_data = json.load(f)
            with open(tmp_layout_json) as f:
                circos_layout = json.load(f)
        except Exception as e:
            print(f"Error loading JSON data: {e}")
            return html.Div(["Error loading plot data."])

        # Check if the loaded data is not empty
        if not circos_graph_data or not circos_layout:
            return html.Div(["No valid data found for the BLAST search."])

        print("Circos graph data:", circos_graph_data)
        print("Circos layout data:", circos_layout)

        # Minimal Circos plot configuration
        try:
            circos_plot = dashbio.Circos(
                layout=circos_layout,
                config={"innerRadius": 200, "outerRadius": 300},
                tracks=[
                    {
                        "type": "chords",
                        "data": circos_graph_data,
                        "config": {
                            "opacity": 0.7,
                            "color": {
                                "name": "source",
                                "field": "id",
                            },
                            "tooltipContent": {
                                "source": "source.id",
                                "target": "target.id",
                                "fields": [
                                    {
                                        "field": "value.pident",
                                        "name": "Identity (%)",
                                    },
                                    {"field": "value.evalue", "name": "E-value"},
                                    {
                                        "field": "value.bitscore",
                                        "name": "Bit Score",
                                    },
                                ],
                            },
                        },
                    },
                ],
                id="blast-chord-plot",
            )

            return html.Div(circos_plot)

        except Exception as e:
            print(f"Error creating Circos plot: {e}")
            return html.Div(["Error creating plot."])
    else:
        return html.Div(["No results found for the BLAST search."])


# TODO: link in ship classification information for subjects here
def blast_table(ship_blast_results):
    tbl = html.Div(
        [
            dash_table.DataTable(
                columns=[
                    {
                        "name": i,
                        "id": i,
                        "deletable": False,
                        "selectable": True,
                    }
                    for i in ship_blast_results.columns
                ],
                data=ship_blast_results.to_dict("records"),
                hidden_columns=["qseqid", "qseq", "sseq"],
                id="ship-blast-table",
                editable=False,
                sort_action="native",
                sort_by=[{"column_id": "pident", "direction": "desc"}],
                sort_mode="single",
                row_selectable="single",
                selected_rows=[],
                row_deletable=False,
                selected_columns=[],
                page_action="native",
                page_current=0,
                page_size=10,
                export_format="tsv",
                css=[{"selector": ".show-hide", "rule": "display: none"}],
                style_table={
                    "overflowX": "auto",
                    "maxWidth": "100%",
                },
                style_cell={
                    "minWidth": "150px",
                    "width": "150px",
                    "maxWidth": "150px",
                    "whiteSpace": "normal",
                },
            ),
            dbc.Button(
                "Download BLAST results",
                id="blast-dl-button",
                n_clicks=0,
                style={"textAlign": "center"},
            ),
            dcc.Download(id="blast-dl"),
        ]
    )
    return tbl


def gene_hmmsearch(hmmer_results):
    # Convert "evalue" column to numeric with errors='coerce'
    hmmer_results["evalue"] = pd.to_numeric(hmmer_results["evalue"], errors="coerce")

    # Drop rows where "evalue" is NaN
    hmmer_results.dropna(subset=["evalue"], inplace=True)

    # Find the index of the minimum "evalue" for each "query_id"
    idx_min_evalue = hmmer_results.groupby("query_id")["evalue"].idxmin()

    # Get corresponding "hit_IDs" for the minimum "evalue" indices
    try:
        superfamily = hmmer_results.loc[idx_min_evalue, "hit_IDs"].iloc[0]
    except IndexError:
        superfamily = None

    if superfamily is not None:
        return superfamily


# Callback to update information about selected row
@callback(
    Output("ship-aln-container", "children"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
    ],
)
def blast_alignments(ship_blast_results, selected_row):
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=True)

    # Convert ship_blast_results to a Pandas DataFrame
    ship_blast_results_df = pd.DataFrame(ship_blast_results)

    if selected_row:
        try:
            row = ship_blast_results_df.iloc[selected_row]
        except IndexError:
            row = None
        if row is not None:
            # Write the specific row to the file
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


# Define callback to download TSV
@callback(
    Output("blast-dl", "data"),
    [Input("blast-dl-button", "n_clicks")],
    [State("ship-blast-table", "data"), State("ship-blast-table", "columns")],
)
def download_tsv(n_clicks, rows, columns):
    if n_clicks == 0:
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
        return f"Error: {str(e)}"


@callback(
    Output("ship-blast-table", "selected_rows"),
    Input(
        "ship-blast-table", "data"
    ),  # This will trigger when the table data is loaded
)
def update_selected_rows(data):
    # Default to selecting the first row after sorting
    if data:
        return [0]  # Selects the first row (0-indexed)
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

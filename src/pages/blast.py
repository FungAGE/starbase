# import subprocess
import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
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
import urllib.request as urlreq
import pandas as pd

from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

dash.register_page(__name__)

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            align="center",
            children=[
                dbc.Col(
                    style={"border": "1px solid black", "padding": "20px"},
                    xs=12,
                    lg=4,
                    className="align-self-center",
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
                                    id="query-upload",
                                    children=html.Div(id="query-sequence-upload"),
                                    style={
                                        "width": "100%",
                                        "height": "75px",
                                        "lineHeight": "75px",
                                        "borderWidth": "2px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "justify-content": "center",
                                    },
                                    multiple=False,
                                    accept=".fa, .fas, .fasta, .fna",
                                ),
                                dbc.Button(
                                    html.H4("Submit BLAST"),
                                    id="submit-button",
                                    n_clicks=0,
                                    style={"textAlign": "center"},
                                    className="d-grid gap-2 col-6 mx-auto",
                                ),
                            ],
                            gap=3,
                            direction="vertical",
                        )
                    ],
                ),
                dbc.Col(
                    xs=12,
                    lg=8,
                    className="align-self-center",
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
                    ],
                ),
            ],
        ),
    ],
)

tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast").name


def parse_fasta(contents, filename):
    nseq = 1
    for sequence in contents:
        nseq += 1

    # Update the height attribute of the style based on the content height
    return [
        html.Div(
            [
                html.H6(f"File name: {filename}", style={"textAlign": "center"}),
                html.H6(f"Number of sequences: {nseq}", style={"textAlign": "center"}),
            ],
        )
    ]


@callback(
    Output("query-sequence-upload", "children"),
    [
        Input("query-upload", "contents"),
        Input("query-upload", "filename"),
    ],
)
def update_fasta_upload(seq_content, seq_filename):
    if seq_content is None:
        # Return the default style if no content is uploaded
        return [
            html.Div(html.P(["Upload a FASTA file"], style={"textAlign": "center"}))
        ]
    else:
        try:
            # "," is the delimeter for splitting content_type from content_string
            content_type, content_string = seq_content.split(",")
            query_string = str(base64.b64decode(seq_content))
            sequences = SeqIO.parse(io.StringIO(query_string), "fasta")
            children = parse_fasta(sequences, seq_filename)
            return children

        except Exception as e:
            print(e)
            return html.Div(["There was an error processing this file."])


@callback(
    Output("output-container", "children"),
    [
        Input("submit-button", "n_clicks"),
        Input("query-text", "value"),
        Input("query-upload", "contents"),
        # Input("query-upload", "filename"),
    ],
)
def run_main(n_clicks, query_text_input, query_file_contents):
    global tmp_blast
    if n_clicks > 0:
        input_type, query_header, query_seq = check_input(
            query_text_input,
            query_file_contents,
        )

        query_type = guess_seq_type(query_seq)

        tmp_query_fasta = tempfile.NamedTemporaryFile(suffix=".fa").name
        cleaned_query_seq = SeqRecord(Seq(query_seq), id=query_header, description="")
        SeqIO.write(cleaned_query_seq, tmp_query_fasta, "fasta")

        blast_results = run_blast(
            seq_type=query_type,
            db_type="ship",
            tmp_query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=0.01,
            threads=2,
        )

        hmmer_results = run_hmmer(
            seq_type=query_type,
            input_eval=0.01,
            query_fasta=tmp_query_fasta,
            threads=2,
        )

        # create results table with alignments
        ship_blast_table = blast_table(blast_results)

        if hmmer_results is not None:
            superfamily_card = hmmer_table(hmmer_results)
        else:
            superfamily_card = """"""

        # chord = blast_chords()

        return html.Div(
            [
                dbc.Stack(
                    [
                        superfamily_card,
                        ship_blast_table,
                    ],
                    gap=3,
                )
            ],
            id="ship-blast-table-container",
        )

    else:
        """"""


def check_input(query_text_input, query_file_contents):
    if query_text_input in ("", None) and query_file_contents is None:
        return html.Div(
            [
                html.Div("No input provided:", style={"color": "red"}),
                html.Div(
                    "Please provide either a query sequence in the text box, or upload a fasta file",
                    style={"color": "red"},
                ),
            ]
        )

    elif query_text_input and query_file_contents:
        return html.Div(
            [
                html.Div("Multiple inputs:", style={"color": "red"}),
                html.Div(
                    "Please provide either a query sequence in the text box, or upload a fasta file, not both",
                    style={"color": "red"},
                ),
            ]
        )

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
        # Count characters for deciding type later
        nucl_count = sum(query_seq.count(nuc) for nuc in nucl_char)
        prot_count = sum(query_seq.count(aa) for aa in prot_char)
        # Guess if sequence is nucleotide or protein
        # if prot_count > nucl_count:
        if prot_count >= (0.1 * len(query_seq)):
            query_type = "prot"
            print("Query is protein sequence")
        else:
            query_type = "nucl"
            print("Query is nucleotide sequence")

        return query_type
    else:
        # Handle the case where query_seq is None
        print("Error: query_seq is None")
        return None


def run_blast(
    seq_type=None,
    db_type=None,
    tmp_query_fasta=None,
    tmp_blast=None,
    input_eval=None,
    threads=None,
):

    print("starting BLAST")

    # Define paths to databases
    db_list = {
        "ship": {"nucl": "database_folder/Starships/ships/fna/blastdb/concatenated.fa"},
        "gene": {
            "tyr": {
                "prot": "database_folder/Starships/captain/tyr/faa/blastdb/concatenated.faa"
            },
        },
    }

    # Determine blast program and database based on sequence type and blast type
    if db_type == "ship":
        blastdb = db_list["ship"]["nucl"]

    if seq_type == "nucl":
        blast_program = NcbiblastnCommandline
    else:
        blast_program = NcbitblastnCommandline

    # if db_type == "gene":
    #     if seq_type == "nucl":
    #         blast_program = NcbiblastnCommandline
    #         blastdb = db_list["gene"][gene_type]["nucl"]
    #     else:
    #         blast_program = NcbiblastpCommandline
    #         blastdb = db_list["gene"][gene_type]["prot"]

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


def run_hmmer(seq_type=None, input_eval=None, query_fasta=None, threads=None):
    hmmer_program = "hmmsearch"
    hmmer_db = "database_folder/Starships/captain/tyr/hmm/YRsuperfams.p1-512.hmm"

    # Run HMMER search
    tmp_hmmer = tempfile.NamedTemporaryFile(suffix=".hmmer.txt").name
    hmmer_cmd = f"{hmmer_program} -o {tmp_hmmer} --cpu {threads} --domE {input_eval} {hmmer_db} {query_fasta}"
    print("Running hmmersearch...")
    subprocess.run(hmmer_cmd, shell=True)

    # Parse HMMER output
    tmp_hmmer_parsed = tempfile.NamedTemporaryFile(suffix=".hmmer.parsed.txt").name
    n_records = parse_hmmer(tmp_hmmer, tmp_hmmer_parsed)
    print(f"Number of hmmersearch records: {n_records}")

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


def blast_chords():
    print("Making circos plot...")
    # circos_data = df[["qseqid", "sseqid", "qstart", "qend"]]
    # print(df)

    data = (
        urlreq.urlopen("https://git.io/circos_graph_data.json").read().decode("utf-8")
    )

    circos_graph_data = json.loads(data)

    layout_config = {
        "labels": {
            "size": 10,
            "color": "#4d4d4d",
        },
        "ticks": {
            "color": "#4d4d4d",
            "labelColor": "#4d4d4d",
            "spacing": 10000000,
            "labelSuffix": "Mb",
            "labelDenominator": 1000000,
            "labelSize": 10,
        },
    }

    chords_config = {"color": "RdYlBu", "opacity": 0.8}

    circos_plot = dashbio.Circos(
        layout=circos_graph_data["GRCh37"],
        config=layout_config,
        tracks=[
            {
                "type": "CHORDS",
                "data": circos_graph_data["chords"],
                "config": chords_config,
            }
        ],
    )

    return circos_plot


# TODO: link in ship classification information for subjects here
def blast_table(ship_blast_results):
    tbl_dat = ship_blast_results
    tbl = html.Div(
        [
            dbc.Button(
                "Download BLAST results",
                id="blast-dl-button",
                n_clicks=0,
                style={"textAlign": "center"},
            ),
            dcc.Download(id="blast-dl"),
            html.Br(),
            dash_table.DataTable(
                columns=[
                    {
                        "name": i,
                        "id": i,
                        "deletable": False,
                        "selectable": True,
                    }
                    for i in tbl_dat.columns
                ],
                data=tbl_dat.to_dict("records"),
                hidden_columns=["qseqid", "qseq", "sseq"],
                id="ship-blast-table",
                editable=False,
                sort_action="native",
                sort_mode="multi",
                row_selectable="single",
                selected_rows=[0],  # Select the first row by default
                row_deletable=False,
                selected_columns=[],
                page_action="native",
                page_current=0,
                page_size=10,
                export_format="tsv",
                css=[{"selector": ".show-hide", "rule": "display: none"}],
            ),
        ]
    )
    return tbl


def hmmer_table(hmmer_results):
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

    if superfamily is None:
        return None
    else:
        superfamily_card = html.Div(
            children=[
                dbc.Card(
                    [
                        dbc.CardHeader(
                            html.H4(f"Likely captain superfamily: {superfamily}")
                        ),
                    ],
                    color="secondary",
                    inverse=True,
                ),
                html.Hr(),
            ],
            style={"width": "33%"},
        )
        return superfamily_card


# Callback to update information about selected row
@callback(
    Output("ship-aln-container", "children"),
    [
        Input("ship-blast-table", "derived_virtual_data"),
        Input("ship-blast-table", "derived_virtual_selected_rows"),
    ],
)
def blast_alignments(ship_blast_results, selected_row):

    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False)

    # Convert ship_blast_results to a Pandas DataFrame
    ship_blast_results_df = pd.DataFrame(ship_blast_results)

    if selected_row is not None:
        row_sel = selected_row
    else:
        row_sel = 0

    try:
        row = ship_blast_results_df.iloc[row_sel]
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

        aln = dashbio.AlignmentChart(
            id="alignment-viewer",
            data=data,
            height=200,
            tilewidth=30,
            # overview="slider",
            showconsensus=False,
            showconservation=False,
            showgap=False,
            showid=False,
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

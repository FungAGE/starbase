# import subprocess
import tempfile
import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State


dash.register_page(__name__)

layout = html.Div(
    [
        html.H1("BLAST/hmmersearch starbase"),
        html.H2(
            "Search protein/nucleotide sequences for Starships and Starship-associated genes."
        ),
        html.Br(),
        dbc.Row(
            [
                dbc.CardHeader("Upload FASTA Sequence:"),
                dbc.CardBody(
                    [
                        dbc.Card(
                            [
                                html.Div(
                                    [
                                        html.Table(
                                            style={"width": "100%"},
                                            children=[
                                                html.Tr(
                                                    [
                                                        html.Td(
                                                            style={"width": "50%"},
                                                            children=[
                                                                dcc.Textarea(
                                                                    id="query_text",
                                                                    placeholder="Paste FASTA sequence here...",
                                                                    rows=5,
                                                                    style={
                                                                        "width": "75%"
                                                                    },
                                                                )
                                                            ],
                                                        ),
                                                        html.Td(
                                                            style={"width": "50%"},
                                                            children=[
                                                                dcc.Upload(
                                                                    id="query_file",
                                                                    children=html.Button(
                                                                        "Upload FASTA File"
                                                                    ),
                                                                )
                                                            ],
                                                        ),
                                                    ]
                                                )
                                            ],
                                        )
                                    ]
                                ),
                                html.Div(id="output-container"),
                            ]
                        )
                    ]
                ),
            ]
        ),
        html.Br(),
        dbc.Button("Submit", id="submit-button", n_clicks=0),
    ]
)


@callback(
    [Output("submit-button", "children"), Output("output-container", "children")],
    [Input("submit-button", "n_clicks")],
    [
        State("query_text", "value"),
        State("query_file", "contents"),
        State("query_file", "filename"),
    ],
)
def check_input(n_clicks, query_text_input, query_file_contents, query_filename):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    if n_clicks is None:
        """"""
    else:

        if (
            query_text_input == "" and query_text_input is None
        ) and query_file_contents is None:
            return html.Div(
                [
                    html.Div("No input provided:", style={"color": "red"}),
                    html.Div(
                        "Please provide either a query sequence in the text box, or upload a fasta file",
                        style={"color": "red"},
                    ),
                ]
            )
        elif (
            query_text_input != "" and query_text_input is not None
        ) and query_file_contents is not None:
            return html.Div(
                [
                    html.Div("Multiple inputs:", style={"color": "red"}),
                    html.Div(
                        "Please provide either a query sequence in the text box, or upload a fasta file, not both",
                        style={"color": "red"},
                    ),
                ]
            )
        elif query_text_input is not None and query_text_input != "":
            input_type = "text"
            query = query_text_input
        elif query_file_contents is not None:
            input_type = "file"
            try:
                if query_filename.endswith({".fa", ".fasta"}):
                    query = query_file_contents
            except Exception as e:
                print(e)
                return html.Div(
                    "There was an error processing this file.", style={"color": "red"}
                )
        else:
            return html.Div(
                [
                    html.Div("Error with reading input:", style={"color": "red"}),
                    html.Div(
                        "This was not supposed to happen... Try again?",
                        style={"color": "red"},
                    ),
                ]
            )

        query_header = parse_fasta_header(query)

        query_list = clean_lines(query)

        query_seq = SeqRecord(
            Seq(query_list["cleaned_query"]), id=query_header, description=""
        )

        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa").name
        SeqIO.write(query_seq, tmp_fasta, "fasta")

        if n_clicks > 0:
            blast_table = run_blast(
                seq_type=query_list["query_type"],
                blast_type="blastn",
                tmp_fasta=tmp_fasta,
                input_eval=0.01,
                threads=2,
            )
            blast_table_ui = html.Div(
                [
                    html.Div("Parsed Input Details:"),
                    html.Div(f"Temporary file path: {tmp_fasta}"),
                    html.Div(f"Query type: {query_list['query_type']}"),
                    html.Div(f"Query header: {query_header}"),
                    html.Div(f"Input type: {input_type}"),
                    html.Div(
                        [
                            dash_table.DataTable(
                                id="blast-table",
                                data=blast_table,
                                editable=False,
                                filter_action="native",
                                sort_action="native",
                                sort_mode="multi",
                                column_selectable="single",
                                row_selectable="multi",
                                row_deletable=True,
                                selected_columns=[],
                                selected_rows=[],
                                page_action="native",
                                page_current=0,
                                page_size=10,
                            )
                        ]
                    ),
                    html.Div(id="table-container"),
                ]
            )

            return (
                dbc.Spinner(html.Span("Processing...")),
                "Button clicked!",
                blast_table_ui,
            )

        else:
            return "Submit", ""


def parse_fasta_header(query):
    header_grep = query.count(">")

    if header_grep > 1:
        return html.Div(
            [
                html.Div(
                    "Multi-fastas not supported at this time", style={"color": "red"}
                ),
                html.Div(
                    "Please provide one sequence at a time.", style={"color": "red"}
                ),
            ]
        )
    elif header_grep == 0:
        query_header = "QUERY"
    elif header_grep == 1:
        query_header = query.split("\n")[0][1:]
    else:
        return html.Div(
            [html.Div("Error reading header from input", style={"color": "red"})]
        )

    return query_header


def clean_lines(query):
    import re

    cleaned_seq = ""
    nucl_char = set("ATGC")
    prot_char = set("ARNDBCEQZGHILKMFPSTWYV")

    for line in query:
        cleaned_seq += re.sub("[^A-Za-z]", "", str(line))

    # Count characters for deciding type later
    nucl_count = sum(cleaned_seq.count(base) for base in nucl_char)
    prot_count = sum(cleaned_seq.count(aa) for aa in prot_char)

    # Guess if sequence is nucleotide or protein
    if prot_count > nucl_count:
        # if (prot_count >= (0.1 * len(cleaned_seq))):
        query_type = "prot"
        print("Query is protein sequence")
    else:
        query_type = "nucl"
        print("Query is nucleotide sequence")

    return {"query": query, "query_type": query_type, "cleaned_query": cleaned_seq}


def run_blast(
    seq_type=None, blast_type=None, tmp_fasta=None, input_eval=None, threads=None
):

    from Bio.Blast.Applications import (
        NcbiblastnCommandline,
        NcbiblastpCommandline,
        NcbitblastnCommandline,
    )
    import tempfile
    import pandas as pd

    # Define paths to databases
    db_list = {
        "ship": {"nucl": "project-vol/Starships/ships/fna/blastdb/concatenated.fa"},
        "gene": {
            "tyr": {
                "prot": "project-vol/Starships/captain/tyr/faa/blastdb/concatenated.faa"
            },
            "fre": {
                "prot": "project-vol/Starships/cargo/fre/faa/blastdb/fre.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/fre/fna/blastdb/fre.fa",
            },
            "nlr": {
                "prot": "project-vol/Starships/cargo/nlr/faa/blastdb/nlr.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/nlr/fna/blastdb/nlr.fa",
            },
            "DUF3723": {
                "prot": "project-vol/Starships/cargo/duf3723/faa/blastdb/duf3723.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/duf3723/fna/blastdb/duf3723.fa",
            },
            "plp": {
                "prot": "project-vol/Starships/cargo/plp/faa/blastdb/plp.mycoDB.faa",
                "nucl": "project-vol/Starships/cargo/plp/fna/blastdb/plp.fa",
            },
        },
    }

    print("Running BLAST...")
    # Determine blast program and database based on sequence type and blast type
    if seq_type == "nucl":
        if blast_type == "ship":
            blast_program = NcbitblastnCommandline
            blastdb = db_list["ship"]["nucl"]
        else:
            blast_program = NcbiblastnCommandline
            blastdb = db_list["gene"][blast_type]["nucl"]
    else:
        if blast_type == "ship":
            blast_program = NcbitblastnCommandline
            blastdb = db_list["ship"]["nucl"]
        else:
            blast_program = NcbiblastpCommandline
            blastdb = db_list["gene"][blast_type]["prot"]

    if len(blastdb) != 1:
        raise ValueError("Issue accessing BLAST database")

    # Perform the BLAST search
    blast_tmp = tempfile.NamedTemporaryFile(suffix=".blast").name
    blast_cline = blast_program(
        query=tmp_fasta,
        db=blastdb,
        evalue=input_eval,
        out=blast_tmp,
        outfmt=6,
        num_threads=threads,
    )
    stdout, stderr = blast_cline()

    # Read the BLAST results into a DataFrame
    blast_results = pd.read_csv(
        blast_tmp,
        sep="\t",
        names=[
            "query_id",
            "hit_IDs",
            "pident",
            "aln_length",
            "mismatches",
            "gaps",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "evalue",
            "bitscore",
            "query_seq",
            "subject_seq",
        ],
    )

    # Optionally stitch BLAST results
    # if stitch:
    #     stitched_blast_tmp = tempfile.NamedTemporaryFile(suffix=".stitch").name
    #     stitch_blast_cmd = f"python bin/BLASTstitcher.py -i {blast_tmp} -o {stitched_blast_tmp}"
    #     subprocess.run(stitch_blast_cmd, shell=True)
    #     ship_blast_out = stitched_blast_tmp

    return html.Div(
        [
            html.H3("Table for metadata of all Starships in starbase"),
            dash_table.DataTable(
                id="blast_table",
                columns=[
                    {"name": i, "id": i, "deletable": False, "selectable": True}
                    for i in blast_results.columns
                ],
                data=blast_results.to_dict("records"),
                editable=True,
                filter_action="native",
                sort_action="native",
                sort_mode="multi",
                column_selectable="single",
                row_selectable="multi",
                row_deletable=True,
                selected_columns=[],
                selected_rows=[],
                page_action="native",
                page_current=0,
                page_size=10,
            ),
            html.Div(id="blast-table-container"),
        ],
        style={"display": "inline-block", "width": "50%"},
    )

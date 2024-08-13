import dash
from dash import dcc, callback, html
from dash.dependencies import Output, Input
import base64
import pandas as pd
from Bio import SeqIO
import io
from io import StringIO


def parse_gff(contents, filename):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    gff = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep="\t")
    # cols = [
    #     "seqid",
    #     "source",
    #     "type",
    #     "start",
    #     "end",
    #     "score",
    #     "strand",
    #     "phase",
    #     "attributes",
    # ]

    # annotations = pd.DataFrame(gff)
    nanno = len(gff)

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of annotations: {nanno}"),
            ]
        )
    ]


def parse_fasta(contents, filename):
    sequences = SeqIO.parse(io.StringIO(contents), "fasta")

    records = []
    nseq = 0
    for sequence in sequences:
        records.append({"ID": sequence.id, "Sequence": str(sequence.seq)})
        nseq += 1

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of sequences: {nseq}"),
            ],
        )
    ]


def dl_package(app):
    @app.callback(Output("dl-package", "data"), [Input("dl-button", "n_clicks")])
    def generate_download(n_clicks):
        if n_clicks is None:
            return dash.no_update
        return dcc.send_file(
            "database_folder/Starships/ships/fna/blastdb/concatenated.fa"
        )


def update_fasta_upload(app):
    @app.callback(
        Output("fasta-sequence-upload", "children"),
        [
            Input("fasta-upload", "contents"),
            Input("fasta-upload", "filename"),
        ],
    )
    def update_fasta_details(seq_content, seq_filename):
        if seq_content is None:
            return [
                html.Div(
                    html.P(
                        ["Select a FASTA file to upload"],
                        style={
                            "justify-content": "center",
                            "textAlign": "center",
                        },
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
                print(e)
                return html.Div(["There was an error processing this file."])


def update_gff_upload(app):
    @app.callback(
        Output("output-gff-upload", "children"),
        [
            Input("upload-gff", "contents"),
            Input("upload-gff", "filename"),
        ],
    )
    def update_gff_details(anno_content, anno_filename):
        if anno_content is None:
            return [
                html.Div(["Select a GFF file to upload"]),
            ]
        else:
            try:
                children = parse_gff(anno_content, anno_filename)
                return children

            except Exception as e:
                print(e)
                return html.Div(["There was an error processing this file."])

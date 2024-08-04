import base64

import dash
import dash_bootstrap_components as dbc
from dash.dependencies import Output, Input
from dash import dcc, html, callback
import sqlite3

from Bio import SeqIO

import time
import datetime
import io

import pandas as pd

dash.register_page(__name__)

# Define form layout
form = html.Div(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            className="align-self-center",
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            html.H2(
                                                [
                                                    "Submission of multiple Starships to ",
                                                    html.Span(
                                                        "starbase",
                                                        className="logo-text",
                                                    ),
                                                ],
                                            )
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.Div(
                                                    [
                                                        html.P(
                                                            [
                                                                "Unfortunately, we can only handle submission for one Starship at a time. If you have a batch of Starships that you'd like to submit, please send the submission via ",
                                                                html.A(
                                                                    "email.",
                                                                    href="mailto:adrian.e.forsythe@gmail.com",
                                                                ),
                                                            ],
                                                            style={
                                                                "fontSize": "1vw",
                                                            },
                                                        ),
                                                    ]
                                                )
                                            ],
                                        ),
                                    ]
                                )
                            ],
                        )
                    ],
                    class_name="mb-3",
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            className="align-self-center",
                            children=[
                                html.H2(
                                    [
                                        "Submit individual Starships to ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
                                ),
                                html.H5(
                                    [
                                        "Fields in ",
                                        html.Span(
                                            "red",
                                            style={"color": "red"},
                                        ),
                                        " = manditory.",
                                    ],
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            className="align-self-center",
                            children=[
                                html.H5(
                                    ["Upload Starship sequence:"],
                                ),
                                dcc.Upload(
                                    id="upload-sequence",
                                    children=html.Div(id="output-sequence-upload"),
                                    style={
                                        "color": "red",
                                        "width": "50%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                        "margin": "10px",
                                    },
                                    multiple=False,
                                    accept=".fa, .fas, .fasta, .fna",
                                ),
                                dcc.Loading(
                                    id="loading-1",
                                    type="default",
                                    children=html.Div(id="loading-output-1"),
                                ),
                                html.H5(
                                    [
                                        "Upload gene annotations associated with Starship sequence (GFF[3] format):"
                                    ],
                                ),
                                dcc.Upload(
                                    id="upload-annotations",
                                    children=html.Div(id="output-annotation-upload"),
                                    accept=".gff, .gff3",
                                    multiple=False,
                                    style={
                                        "width": "50%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "5px",
                                        "textAlign": "center",
                                        "margin": "10px",
                                    },
                                ),
                                dcc.Loading(
                                    id="loading-2",
                                    type="default",
                                    children=html.Div(id="loading-output-2"),
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            className="align-self-center",
                            children=[
                                html.H5(
                                    ["Starship metadata:"],
                                ),
                                html.P(
                                    [
                                        "Name of curator: ",
                                        dcc.Input(
                                            id="uploader",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                                html.P(
                                    [
                                        "How were Starships annotated? (i.e. starfish): ",
                                        dcc.Input(
                                            id="evidence",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                                html.P(
                                    [
                                        "Enter genus name: ",
                                        dcc.Input(
                                            id="genus",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                                html.P(
                                    [
                                        "Enter species name: ",
                                        dcc.Input(
                                            id="species",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            className="align-self-center",
                            children=[
                                html.H5("Coordinates of Starship in host genome:"),
                                html.P(
                                    [
                                        "Host genome contig/scaffold/chromosome ID: ",
                                        dcc.Input(
                                            id="hostchr",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                                html.P(
                                    [
                                        "Start coordinate of Starship: ",
                                        dcc.Input(
                                            id="shipstart",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                                html.P(
                                    [
                                        "End coordinate for Starship: ",
                                        dcc.Input(
                                            id="shipend",
                                            type="text",
                                            style={"width": "15%"},
                                            required=True,
                                        ),
                                    ],
                                    style={"font-size": "1vw"},
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            className="align-self-center",
                            lg=8,
                            sm=12,
                            children=[
                                html.H5("Additional information:"),
                                dcc.Textarea(
                                    id="comment",
                                    placeholder="Any comments about the Starship features, annotations, or host genome?",
                                    style={
                                        "height": "100px",
                                        "width": "50%",
                                    },
                                    required=False,
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            className="align-self-center",
                            lg=8,
                            sm=12,
                            children=[
                                dbc.Button(
                                    html.H5(
                                        ["Submit"],
                                        style={
                                            "align-items": "center",
                                            "justify-content": "center",
                                            "textAlign": "center",
                                        },
                                    ),
                                    id="submit-ship",
                                    n_clicks=0,
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        )
    ],
)

@callback(Output("loading-output-1", "children"), Input("upload-sequence", "value"))
def input_triggers_spinner(value):
    time.sleep(1)
    return value


@callback(Output("loading-output-2", "children"), Input("loading-input-2", "value"))
def input_triggers_nested(value):
    time.sleep(1)
    return value


def parse_fasta(contents, filename):
    if contents:
        # Assume that the user uploaded a FASTA file
        sequences = SeqIO.parse(io.StringIO(contents), "fasta")

        records = []
        nseq = 1
        for sequence in sequences:
            records.append({"ID": sequence.id, "Sequence": str(sequence.seq)})
            nseq += 1

        # Update the height attribute of the style based on the content height
        return [
            html.Div(
                [
                    html.H6(f"File name: {filename}"),
                    html.H6(f"Number of sequences: {nseq}"),
                ],
            )
        ]

    else:
        # Return the default style if no content is uploaded
        return [
            html.Div(
                [
                    "Drag and Drop or ",
                    html.A("Select a File"),
                ],
            )
        ]


def parse_gff(contents, filename):
    if contents:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
        # Assume that the user uploaded a TSV file
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

        return html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of annotations: {nanno}"),
            ]
        )

    else:
        return [
            html.Div(
                [
                    "Drag and Drop or ",
                    html.A("Select a File"),
                ]
            ),
        ]


@callback(
    Output("output-sequence-upload", "children"),
    [
        Input("upload-sequence", "contents"),
        Input("upload-sequence", "filename"),
    ],
)
def update_fasta_upload(seq_content, seq_filename):
    try:
        children = parse_fasta(seq_content, seq_filename)
        return children

    except Exception as e:
        print(e)
        return html.Div(["There was an error processing this file."])


@callback(
    Output("output-annotation-upload", "children"),
    [
        Input("upload-annotations", "contents"),
        Input("upload-annotations", "filename"),
    ],
)
def update_anno_upload(anno_content, anno_filename):
    try:
        children = parse_gff(anno_content, anno_filename)
        return children

    except Exception as e:
        print(e)
        return html.Div(["There was an error processing this file."])


@callback(
    Output("output-data-upload", "children"),
    [
        Input("upload-sequence", "contents"),
        Input("upload-sequence", "filename"),
        Input("upload-sequence", "last_modified"),
        Input("upload-annotations", "contents"),
        Input("upload-annotations", "filename"),
        Input("upload-annotations", "last_modified"),
        Input("uploader", "value"),
        Input("evidence", "value"),
        Input("genus", "value"),
        Input("species", "value"),
        Input("hostchr", "value"),
        Input("shipstart", "value"),
        Input("shipend", "value"),
        Input("comment", "value"),
        Input("submit-ship", "n_clicks"),
    ],
)
def submit_ship(
    seq_content,
    seq_filename,
    seq_date,
    anno_content,
    anno_filename,
    anno_date,
    uploader,
    evidence,
    genus,
    species,
    hostchr,
    shipstart,
    shipend,
    comment,
    n_clicks,
):
    if n_clicks > 0 and seq_content is not None:
        try:
            # Create SQLite database connection
            conn = sqlite3.connect("database_folder/starbase.sqlite")
            c = conn.cursor()

            # Check if the table already exists
            c.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='submissions_new'"
            )
            existing_table = c.fetchone()

        except sqlite3.Error as error:
            print("Error connecting to database.", error)

        if not existing_table:
            # Create table for storing form submissions
            c.execute(
                """CREATE TABLE submissions_new (
                seq_contents TEXT NOT NULL,
                seq_filename TEXT NOT NULL,
                seq_date TEXT NOT NULL,
                anno_contents TEXT,
                anno_filename TEXT,
                anno_date TEXT,
                uploader TEXT NOT NULL,
                evidence TEXT NOT NULL,
                genus TEXT NOT NULL,
                species TEXT NOT NULL,
                hostchr TEXT NOT NULL,
                shipstart INTEGER NOT NULL,
                shipend INTEGER NOT NULL,
                comment TEXT,
                id	INTEGER NOT NULL,
                PRIMARY KEY(id AUTOINCREMENT)
            )"""
            )
            conn.commit()

        try:
            if seq_content is None or len(seq_content) == 0:
                return "No fasta file uploaded"
            else:
                content_type, content_string = seq_content.split(",")
                seq_decoded = base64.b64decode(content_string).decode("utf-8")
                seq_datetime_obj = datetime.datetime.fromtimestamp(seq_date)
                # update_fasta_output(seq_decoded, seq_filename, seq_datetime_obj)

            # ? put annotations into separate SQL submission table?
            if anno_content is not None and len(anno_content) > 0:
                content_type, content_string = anno_content.split(",")
                decoded = base64.b64decode(content_string).decode("utf-8")
                anno_datetime_obj = datetime.datetime.fromtimestamp(anno_date)
                # update_gff_output(decoded, anno_filename, anno_datetime_obj)
            else:
                anno_contents = ""
                anno_filename = ""
                anno_date = ""

            if comment is None:
                comment = ""

            c.execute(
                "INSERT INTO submissions_new (seq_contents,seq_filename,seq_date,anno_contents,anno_filename,anno_date,uploader,evidence,genus,species,hostchr,shipstart,shipend,comment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    seq_content,
                    seq_filename,
                    seq_datetime_obj,
                    anno_contents,
                    anno_filename,
                    anno_datetime_obj,
                    uploader,
                    evidence,
                    genus,
                    species,
                    hostchr,
                    shipstart,
                    shipend,
                    comment,
                ),
            )
            conn.commit()

            return html.Div(
                [html.H5(f"Successfully submitted '{seq_filename}' to starbase")]
            )

        except sqlite3.Error as error:
            print("Failed to insert record into SQLite table:", error)
        finally:
            c.close()
            conn.close()
            print("SQLite connection is closed.")

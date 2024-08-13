import base64

import dash
import dash_bootstrap_components as dbc
from dash.dependencies import Output, Input, State
from dash import dcc, html, callback
import sqlite3

from Bio import SeqIO

import time
import datetime
import io

import pandas as pd

dash.register_page(__name__)


layout = html.Div(
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
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            html.H3(
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
                                                                "fontSize": "1rem",
                                                            },
                                                        ),
                                                    ]
                                                )
                                            ],
                                        ),
                                    ],
                                    className="w-100 w-lg-75",
                                )
                            ],
                        )
                    ],
                    class_name="mb-3",
                    style={"paddingTop": "20px"},
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            children=[
                                html.H3(
                                    [
                                        "Submit individual Starships to ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
                                    className="text-center",
                                ),
                                html.H4(
                                    [
                                        "Fields in ",
                                        html.Span(
                                            "red",
                                            style={"color": "red"},
                                        ),
                                        " = manditory.",
                                    ],
                                    className="text-center",
                                ),
                            ],
                        )
                    ],
                    style={"padding": "10px"},
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            children=[
                                html.H4(
                                    ["Upload Starship sequence:"],
                                ),
                                dcc.Upload(
                                    id="fasta-upload",
                                    children=html.Div(id="fasta-sequence-upload"),
                                    style={
                                        "color": "red",
                                        "width": "100%",
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
                                html.H4(
                                    [
                                        "Upload gene annotations associated with Starship sequence (GFF[3] format):"
                                    ],
                                ),
                                dcc.Upload(
                                    id="upload-gff",
                                    children=html.Div(id="output-gff-upload"),
                                    accept=".gff, .gff3, .tsv",
                                    multiple=False,
                                    style={
                                        "width": "100%",
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
                            children=[
                                html.H4(
                                    ["Starship Metadata"],
                                ),
                                html.P(
                                    [
                                        "Email of curator: ",
                                        dcc.Input(
                                            id="uploader",
                                            type="email",
                                            style={"width": "100%"},
                                            className="form-control",
                                            placeholder="Enter email",
                                            required=True,
                                        ),
                                    ],
                                ),
                                html.P(
                                    [
                                        "How were Starships annotated? (i.e. starfish): ",
                                        dcc.Input(
                                            id="evidence",
                                            type="text",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                                html.P(
                                    [
                                        "Enter genus name: ",
                                        dcc.Input(
                                            id="genus",
                                            type="text",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                                html.P(
                                    [
                                        "Enter species name: ",
                                        dcc.Input(
                                            id="species",
                                            type="text",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                            ],
                        )
                    ],
                    style={"padding": "10px"},
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            children=[
                                html.H4("Coordinates of Starship in host genome:"),
                                html.P(
                                    [
                                        "Host genome contig/scaffold/chromosome ID: ",
                                        dcc.Input(
                                            id="hostchr",
                                            type="text",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                                html.P(
                                    [
                                        "Start coordinate of Starship: ",
                                        dcc.Input(
                                            id="shipstart",
                                            type="number",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                                html.P(
                                    [
                                        "End coordinate for Starship: ",
                                        dcc.Input(
                                            id="shipend",
                                            type="number",
                                            style={"width": "100%"},
                                            className="form-control",
                                            required=True,
                                        ),
                                    ],
                                ),
                            ],
                        )
                    ],
                    style={"padding": "10px"},
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            children=[
                                html.H4("Additional information:"),
                                dcc.Textarea(
                                    id="comment",
                                    placeholder="Any comments about the Starship features, annotations, or host genome?",
                                    style={
                                        "height": "100px",
                                        "width": "100%",
                                    },
                                    required=False,
                                ),
                            ],
                        )
                    ],
                    style={"padding": "10px"},
                ),
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            children=[
                                dbc.Button(
                                    html.H4(
                                        ["Submit"],
                                        className="text-center",
                                    ),
                                    id="submit-ship",
                                    n_clicks=0,
                                    className="d-grid gap-2 col-6 mx-auto",
                                ),
                                dbc.Modal(
                                    [
                                        dbc.ModalHeader(
                                            dbc.ModalTitle("New Submission")
                                        ),
                                        dbc.ModalBody(
                                            html.Div(id="output-data-upload")
                                        ),
                                        dbc.ModalFooter(
                                            dbc.Button(
                                                "Close",
                                                id="close",
                                                className="ms-auto",
                                                n_clicks=0,
                                            )
                                        ),
                                    ],
                                    id="modal",
                                    is_open=False,
                                ),
                            ],
                        ),
                    ],
                    style={"padding": "10px"},
                ),
            ],
        )
    ],
)


@callback(
    [Output("modal", "is_open"), Output("output-data-upload", "children")],
    [
        Input("upload-fasta", "contents"),
        Input("upload-fasta", "filename"),
        Input("upload-fasta", "last_modified"),
        Input("upload-gff", "contents"),
        Input("upload-gff", "filename"),
        Input("upload-gff", "last_modified"),
        Input("uploader", "value"),
        Input("evidence", "value"),
        Input("genus", "value"),
        Input("species", "value"),
        Input("hostchr", "value"),
        Input("shipstart", "value"),
        Input("shipend", "value"),
        Input("comment", "value"),
        Input("submit-ship", "n_clicks"),
        Input("close", "n_clicks"),
    ],
    [State("modal", "is_open")],
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
    close_modal,
    is_open,
):
    modal = not is_open
    message = None
    if n_clicks > 0:
        if seq_content is None or len(seq_content) == 0:
            return "No fasta file uploaded"
        else:
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
                    anno_datetime_obj = ""

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

                modal = True
                if close_modal:
                    modal = False

                message = html.H5(
                    f"Successfully submitted '{seq_filename}' to starbase"
                )

                return modal, message

            except sqlite3.Error as error:
                print("Failed to insert record into SQLite table:", error)

            finally:
                c.close()
                conn.close()
                print("SQLite connection is closed.")

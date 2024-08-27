import warnings

warnings.filterwarnings("ignore")

import base64

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback
import sqlite3

import datetime

dash.register_page(__name__)


layout = dmc.Container(
    fluid=True,
    children=[
        dbc.Form(
            [
                dmc.Grid(
                    justify="start",
                    align="center",
                    children=[
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
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
                                    className="auto-resize-600",
                                )
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                dmc.Title(
                                    [
                                        "Submit individual Starships to ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
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
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                html.H4(
                                    ["Upload Starship sequence:"],
                                ),
                                dcc.Upload(
                                    id="fasta-upload",
                                    children=html.Div(id="fasta-sequence-upload"),
                                    className="upload-box text-red auto-resize-600 text-center",
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
                                    className="upload-box text-red auto-resize-600 text-center",
                                ),
                                dcc.Loading(
                                    id="loading-2",
                                    type="default",
                                    children=html.Div(id="loading-output-2"),
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                html.H4(
                                    ["Starship Metadata"],
                                ),
                                dbc.Label(
                                    "Email of curator: ",
                                    html_for="uploader",
                                ),
                                dcc.Input(
                                    id="uploader",
                                    type="email",
                                    className="form-control auto-resize-600",
                                    placeholder="Enter email",
                                    required=True,
                                ),
                                dbc.Label("How were Starships annotated?"),
                                dcc.Input(
                                    id="evidence",
                                    type="text",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="i.e. starfish",
                                ),
                                dbc.Label("Enter genus name:"),
                                dcc.Input(
                                    id="genus",
                                    type="text",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="Alternaria",
                                ),
                                dbc.Label("Enter species name:"),
                                dcc.Input(
                                    id="species",
                                    type="text",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="alternata",
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                html.H4("Coordinates of Starship in host genome:"),
                                dbc.Label("Host genome contig/scaffold/chromosome ID:"),
                                dcc.Input(
                                    id="hostchr",
                                    type="text",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="'chr1', or GenBank Accession",
                                ),
                                dbc.Label(
                                    "Start coordinate of Starship (relative to contig/scaffold/chromosome):"
                                ),
                                dcc.Input(
                                    id="shipstart",
                                    type="number",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="i.e. 1200",
                                ),
                                dbc.Label(
                                    "End coordinate for Starship (relative to contig/scaffold/chromosome):"
                                ),
                                dcc.Input(
                                    id="shipend",
                                    type="number",
                                    className="form-control auto-resize-600",
                                    required=True,
                                    placeholder="i.e. 20500",
                                ),
                                dbc.Label("Starship found on strand:"),
                                dbc.RadioItems(
                                    id="strand-radios",
                                    options=[
                                        {
                                            "label": "Positive strand",
                                            "value": 1,
                                        },
                                        {
                                            "label": "Negative strand",
                                            "value": 2,
                                        },
                                    ],
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                dbc.Label("Additional information:"),
                                html.Br(),
                                dcc.Textarea(
                                    id="comment",
                                    placeholder="Any comments about the Starship features, annotations, or host genome?",
                                    style={
                                        "height": "100px",
                                    },
                                    required=False,
                                    className="auto-resize-600",
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={
                                "lg": 6,
                                "sm": 12,
                            },
                            offset={
                                "lg": 3,
                                "sm": 0,
                            },
                            children=[
                                dbc.Button(
                                    html.H4(
                                        ["Submit"],
                                        className="text-center auto-resize-600",
                                    ),
                                    id="submit-ship",
                                    n_clicks=0,
                                    # className="d-grid gap-2 col-6 mx-auto",
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
            ]
        )
    ],
)


@callback(
    [Output("modal", "is_open"), Output("output-data-upload", "children")],
    [
        Input("fasta-upload", "contents"),
        Input("fasta-upload", "filename"),
        Input("fasta-upload", "last_modified"),
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
        Input("strand-radios", "value"),
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
    strand_radio,
    comment,
    n_clicks,
    close_modal,
    is_open,
):
    modal = is_open  # Keep the modal state as it is unless toggled
    message = """"""

    if n_clicks and n_clicks > 0:
        if strand_radio == 1:
            shipstrand = "+"
        else:
            shipstrand = "-"
        if not seq_content:
            return (
                modal,
                "No fasta file uploaded",
            )  # Return the error message if no file

        try:
            # Create SQLite database connection
            conn = sqlite3.connect("database_folder/starbase.sqlite")
            c = conn.cursor()

            # Check if the table already exists
            c.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name='submissions_new'"
            )
            existing_table = c.fetchone()

            if not existing_table:
                # Create the table if it doesn't exist
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
                    shipstrand TEXT,
                    comment TEXT,
                    id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT
                )"""
                )
                conn.commit()

            # Process sequence and annotation contents
            content_type, content_string = seq_content.split(",")
            seq_decoded = base64.b64decode(content_string).decode("utf-8")
            seq_datetime_obj = datetime.datetime.fromtimestamp(seq_date)

            anno_contents = ""
            anno_filename = ""
            anno_datetime_obj = ""

            if anno_content:
                content_type, content_string = anno_content.split(",")
                anno_contents = base64.b64decode(content_string).decode("utf-8")
                anno_datetime_obj = datetime.datetime.fromtimestamp(anno_date)

            if not comment:
                comment = ""

            # Insert submission data into the table
            c.execute(
                """INSERT INTO submissions_new (seq_contents, seq_filename, seq_date, anno_contents,
                anno_filename, anno_date, uploader, evidence, genus, species, hostchr, shipstart,
                shipend, shipstrand, comment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
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
                    shipstrand,
                    comment,
                ),
            )
            conn.commit()

            modal = not is_open if not close_modal else False
            message = html.H5(f"Successfully submitted '{seq_filename}' to starbase")

        except sqlite3.Error as error:
            print("Failed to insert record into SQLite table:", error)
            message = html.H5("Failed to submit data. Please try again.")

        finally:
            if conn:
                c.close()
                conn.close()

    return modal, message

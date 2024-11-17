import base64
import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback

import datetime
from src.utils.parsing import parse_fasta, parse_gff
import logging

logger = logging.getLogger(__name__)

dash.register_page(__name__)


from src.components.sql_engine import submissions_engine, submissions_connected

layout = dmc.Container(
    size="md",
    children=[
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title(
                    [
                        "Submit Starships to ",
                        html.Span("starbase", className="logo-text"),
                    ],
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    "Fields marked with * are required",
                    c="red",
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),

        # Form Content
        dmc.Paper(
            children=dbc.Form([
                dmc.Stack([
                    # File Upload Section
                    dmc.Stack([
                        dmc.Title("Upload Files", order=2, mb="md"),
                        
                        # FASTA Upload
                        dmc.Paper(
                            children=[
                                dmc.Text("Starship Sequence *", fw=500, mb="sm"),
                                dcc.Upload(
                                    id="submit-fasta-upload",
                                    children=html.Div(
                                        id="submit-fasta-sequence-upload",
                                        style={"textAlign": "center", "padding": "20px"}
                                    ),
                                    className="upload-box",
                                    multiple=False,
                                    accept=".fa, .fas, .fasta, .fna",
                                    max_size=10000000,
                                ),
                                dcc.Loading(
                                    id="loading-1",
                                    type="circle",
                                    children=html.Div(id="loading-output-1"),
                                ),
                            ],
                            p="md",
                            radius="md",
                            withBorder=False,
                        ),
                        
                        # GFF Upload
                        dmc.Paper(
                            children=[
                                dmc.Text("Gene Annotations (GFF3)", fw=500, mb="sm"),
                                dcc.Upload(
                                    id="submit-upload-gff",
                                    children=html.Div(
                                        id="submit-output-gff-upload",
                                        style={"textAlign": "center", "padding": "20px"}
                                    ),
                                    className="upload-box",
                                    accept=".gff, .gff3, .tsv",
                                    multiple=False,
                                    max_size=10000000,
                                ),
                                dcc.Loading(
                                    id="loading-2",
                                    type="circle",
                                    children=html.Div(id="loading-output-2"),
                                ),
                            ],
                            p="md",
                            radius="md",
                            withBorder=False,
                        ),
                    ], gap="lg"),

                    # Metadata Section
                    dmc.Stack([
                        dmc.Title("Metadata", order=2, mb="md"),
                        
                        # Curator Info
                        dmc.TextInput(
                            id="uploader",
                            label="Email of curator *",
                            placeholder="Enter email",
                            required=True,
                        ),
                        
                        dmc.TextInput(
                            id="evidence",
                            label="How were Starships annotated? *",
                            placeholder="i.e. starfish",
                            required=True,
                        ),
                        
                        # Taxonomy Info
                        dmc.Group([
                            dmc.TextInput(
                                id="genus",
                                label="Genus *",
                                placeholder="Alternaria",
                                required=True,
                                style={"flex": 1},
                            ),
                            dmc.TextInput(
                                id="species",
                                label="Species *",
                                placeholder="alternata",
                                required=True,
                                style={"flex": 1},
                            ),
                        ]),
                        
                        # Location Info
                        dmc.TextInput(
                            id="hostchr",
                            label="Host genome contig/scaffold/chromosome ID *",
                            placeholder="'chr1', or GenBank Accession",
                            required=True,
                        ),
                        
                        dmc.Group([
                            dmc.NumberInput(
                                id="shipstart",
                                label="Start coordinate *",
                                placeholder="1200",
                                required=True,
                                style={"flex": 1},
                            ),
                            dmc.NumberInput(
                                id="shipend",
                                label="End coordinate *",
                                placeholder="20500",
                                required=True,
                                style={"flex": 1},
                            ),
                        ]),
                        
                        dmc.RadioGroup(
                            id="strand-radios",
                            label="Strand",
                            value=1,
                            children=[
                                dmc.Radio(label="Positive strand", value=1),
                                dmc.Space(h=10),
                                dmc.Radio(label="Negative strand", value=2),
                            ],
                        ),
                        
                        # Comments
                        dmc.Textarea(
                            id="comment",
                            label="Additional information",
                            placeholder="Any comments about the Starship features, annotations, or host genome?",
                            minRows=3,
                        ),
                    ], gap="md"),

                    # Submit Button
                    dmc.Center(
                        dmc.Button(
                            "Submit Starship",
                            id="submit-ship",
                            size="lg",
                            variant="gradient",
                            gradient={"from": "indigo", "to": "cyan"},
                        ),
                    ),
                ], gap="xl"),
            ]),
            p="xl",
            radius="md",
            withBorder=True,
        ),
        
        # Modal
        dbc.Modal([
            dbc.ModalHeader(dbc.ModalTitle("New Submission")),
            dbc.ModalBody(html.Div(id="output-data-upload")),
            dbc.ModalFooter(
                dbc.Button("Close", id="close", className="ms-auto", n_clicks=0)
            ),
        ], id="submit-modal", is_open=False),
    ],
    style={
        "margin": "0 auto",
        "padding": "2rem",
    },
)


# Function to insert a new submission into the database
def insert_submission(
    submissions_engine,
    seq_contents,
    seq_filename,
    seq_date,
    anno_contents,
    anno_filename,
    anno_date,
    uploader,
    evidence,
    genus,
    species,
    hostchr,
    shipstart,
    shipend,
    shipstrand,
    comment,
):

    content_type, content_string = seq_contents.split(",")
    seq_decoded = base64.b64decode(content_string).decode("utf-8")
    seq_datetime_obj = datetime.datetime.fromtimestamp(seq_date).strftime(
        "%Y-%m-%d %H:%M:%S"
    )

    anno_contents = ""
    anno_datetime_obj = ""

    if anno_contents:
        content_type, content_string = anno_contents.split(",")
        anno_contents = base64.b64decode(content_string).decode("utf-8")
        anno_datetime_obj = datetime.datetime.fromtimestamp(anno_date).strftime(
            "%Y-%m-%d %H:%M:%S"
        )

    if not comment:
        comment = ""

    with submissions_engine.connect() as connection:
        query = """
            INSERT INTO submissions (
                seq_contents, seq_filename, seq_date, anno_contents,
                anno_filename, anno_date, uploader, evidence,
                genus, species, hostchr, shipstart, shipend,
                shipstrand, comment
            )
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
        connection.execute(
            query,
            (
                seq_contents,
                seq_filename,
                seq_date,
                anno_contents,
                anno_filename,
                anno_date,
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


@callback(
    [Output("submit-modal", "is_open"), Output("output-data-upload", "children")],
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
        Input("submit-fasta-upload", "last_modified"),
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
        Input("submit-upload-gff", "last_modified"),
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
    [State("submit-modal", "is_open")],
)
def submit_ship(
    seq_contents,
    seq_filename,
    seq_date,
    anno_contents,
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
        if not seq_contents:
            return (
                modal,
                "No fasta file uploaded",
            )  # Return the error message if no file

        insert_submission(
            submissions_engine,
            seq_contents,
            seq_filename,
            seq_date,
            anno_contents,
            anno_filename,
            anno_date,
            uploader,
            evidence,
            genus,
            species,
            hostchr,
            shipstart,
            shipend,
            shipstrand,
            comment,
        )

        modal = not is_open if not close_modal else False
        message = html.H5(f"Successfully submitted '{seq_filename}' to starbase")

    return modal, message


@callback(
    Output("submit-fasta-sequence-upload", "children"),
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
    ],
)
def update_fasta_details(seq_contents, seq_filename):
    if seq_contents is None:
        return [
            html.Div(
                html.P(
                    ["Select a FASTA file to upload"],
                )
            )
        ]
    else:
        try:
            # "," is the delimeter for splitting content_type from content_string
            content_type, content_string = seq_contents.split(",")
            query_string = base64.b64decode(content_string).decode("utf-8")
            children = parse_fasta(query_string, seq_filename)
            return children

        except Exception as e:
            logger.error(e)
            return html.Div(["There was an error processing this file."])


@callback(
    Output("submit-output-gff-upload", "children"),
    [
        Input("submit-upload-gff", "contents"),
        Input("submit-upload-gff", "filename"),
    ],
)
def update_gff_details(anno_contents, anno_filename):
    if anno_contents is None:
        return [
            html.Div(["Select a GFF file to upload"]),
        ]
    else:
        try:
            children = parse_gff(anno_contents, anno_filename)
            return children

        except Exception as e:
            logger.error(e)
            return html.Div(["There was an error processing this file."])

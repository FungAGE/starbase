import base64
import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash.dependencies import Output, Input, State
from dash import dcc, html, callback

import datetime
from src.utils.seq_utils import parse_fasta, parse_gff

from src.database.sql_engine import get_submissions_session
from sqlalchemy.exc import SQLAlchemyError
from sqlalchemy import text
from src.components.callbacks import create_file_upload
from src.utils.classification_utils import assign_accession
from src.database.sql_manager import fetch_ships

from src.config.logging import get_logger

logger = get_logger(__name__)

dash.register_page(__name__)


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
            children=dbc.Form(
                [
                    dmc.Stack(
                        [
                            # File Upload Section
                            dmc.Stack(
                                [
                                    dmc.Title("Upload Files", order=2, mb="md"),
                                    # FASTA Upload
                                    dmc.Paper(
                                        children=[
                                            dmc.Text(
                                                [
                                                    "Starship Sequence ",
                                                    html.Span(
                                                        "*", style={"color": "red"}
                                                    ),
                                                ],
                                                fw=500,
                                                mb="sm",
                                            ),
                                            create_file_upload(
                                                upload_id="submit-fasta-upload",
                                                output_id="submit-fasta-sequence-upload",
                                                accept_types=[
                                                    ".fa",
                                                    ".fas",
                                                    ".fasta",
                                                    ".fna",
                                                ],
                                                placeholder_text="Select a FASTA file to upload",
                                            ),
                                            dcc.Loading(
                                                id="loading-1",
                                                type="circle",
                                                children=html.Div(
                                                    id="loading-output-1"
                                                ),
                                            ),
                                        ],
                                        p="md",
                                        radius="md",
                                        withBorder=False,
                                    ),
                                    # GFF Upload
                                    dmc.Paper(
                                        children=[
                                            dmc.Text(
                                                "Gene Annotations (GFF3)",
                                                fw=500,
                                                mb="sm",
                                            ),
                                            create_file_upload(
                                                upload_id="submit-upload-gff",
                                                output_id="submit-output-gff-upload",
                                                accept_types=[".gff", ".gff3", ".tsv"],
                                                placeholder_text="Select a GFF file to upload",
                                            ),
                                            dcc.Loading(
                                                id="loading-2",
                                                type="circle",
                                                children=html.Div(
                                                    id="loading-output-2"
                                                ),
                                            ),
                                        ],
                                        p="md",
                                        radius="md",
                                        withBorder=False,
                                    ),
                                ],
                                gap="lg",
                            ),
                            # Metadata Section
                            dmc.Stack(
                                [
                                    dmc.Title("Metadata", order=2, mb="md"),
                                    # Curator Info
                                    dmc.TextInput(
                                        id="uploader",
                                        label="Email of curator",
                                        placeholder="Enter email",
                                        required=True,
                                    ),
                                    dmc.TextInput(
                                        id="evidence",
                                        label="How were Starships annotated?",
                                        placeholder="i.e. starfish",
                                        required=True,
                                    ),
                                    # Taxonomy Info
                                    dmc.Group(
                                        [
                                            dmc.TextInput(
                                                id="genus",
                                                label="Genus",
                                                placeholder="Alternaria",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                            dmc.TextInput(
                                                id="species",
                                                label="Species",
                                                placeholder="alternata",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                        ]
                                    ),
                                    # Location Info
                                    dmc.TextInput(
                                        id="hostchr",
                                        label="Host genome contig/scaffold/chromosome ID",
                                        placeholder="'chr1', or GenBank Accession",
                                        required=True,
                                    ),
                                    dmc.Group(
                                        [
                                            dmc.NumberInput(
                                                id="shipstart",
                                                label="Start coordinate",
                                                placeholder="1200",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                            dmc.NumberInput(
                                                id="shipend",
                                                label="End coordinate",
                                                placeholder="20500",
                                                required=True,
                                                style={"flex": 1},
                                            ),
                                        ]
                                    ),
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
                                ],
                                gap="md",
                            ),
                            # Submit Button
                            dmc.Center(
                                dmc.Button(
                                    "Submit Starship",
                                    id="submit-ship",
                                    size="lg",
                                    variant="gradient",
                                    gradient={"from": "indigo", "to": "cyan"},
                                    loading=False,
                                ),
                            ),
                        ],
                        gap="xl",
                    ),
                ]
            ),
            p="xl",
            radius="md",
            withBorder=True,
        ),
        # Modal
        dbc.Modal(
            [
                dbc.ModalHeader(
                    dbc.ModalTitle("Submission Success", className="text-success"),
                    close_button=True,
                ),
                dbc.ModalBody(
                    [
                        html.Div(
                            [
                                html.I(
                                    className="fas fa-check-circle text-success fa-3x mb-3"
                                ),
                                html.Div(id="output-data-upload", className="mt-3"),
                                dmc.Text(
                                    "Your Starship has been successfully added to the database.",
                                    className="text-muted mt-2",
                                ),
                            ],
                            className="text-center",
                        )
                    ]
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close",
                        id="close",
                        className="ms-auto",
                        color="primary",
                        n_clicks=0,
                    )
                ),
            ],
            id="submit-modal",
            is_open=False,
            centered=True,
        ),
    ],
    style={
        "margin": "0 auto",
        "padding": "2rem",
    },
)


# Function to insert a new submission into the database
def insert_submission(
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
    accession=None,
    needs_review=False,
):
    try:
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

        with get_submissions_session() as session:
            query = """
                INSERT INTO submissions (
                    seq_contents, seq_filename, seq_date, anno_contents,
                    anno_filename, anno_date, uploader, evidence,
                    genus, species, hostchr, shipstart, shipend,
                    shipstrand, comment, accession_tag, needs_review
                )
                VALUES (
                    :seq_contents, :seq_filename, :seq_date, :anno_contents, 
                    :anno_filename, :anno_date, :uploader, :evidence, 
                    :genus, :species, :hostchr, :shipstart, :shipend, 
                    :shipstrand, :comment, :accession_tag, :needs_review
                )
            """
            session.execute(
                text(query),
                {
                    "seq_contents": seq_decoded,
                    "seq_filename": seq_filename,
                    "seq_date": seq_datetime_obj,
                    "anno_contents": anno_contents,
                    "anno_filename": anno_filename,
                    "anno_date": anno_datetime_obj,
                    "uploader": uploader,
                    "evidence": evidence,
                    "genus": genus,
                    "species": species,
                    "hostchr": hostchr,
                    "shipstart": shipstart,
                    "shipend": shipend,
                    "shipstrand": shipstrand,
                    "comment": comment,
                    "accession_tag": accession,
                    "needs_review": needs_review,
                },
            )
            session.commit()
            logger.debug(f"Successfully inserted submission for {seq_filename}")
            return True

    except SQLAlchemyError as e:
        logger.error(f"Database error during submission: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error processing submission: {str(e)}")
        raise


def create_fasta_display(records, filename):
    """Create HTML components to display FASTA file info."""
    return html.Div(
        [
            html.H6(f"File name: {filename}"),
            html.H6(f"Number of sequences: {len(records)}"),
        ]
    )


@callback(
    [
        Output("submit-modal", "is_open"),
        Output("output-data-upload", "children"),
        Output("submit-ship", "loading"),
    ],
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
    modal = is_open
    message = ""
    loading = False
    accession = None
    needs_review = True

    if n_clicks and n_clicks > 0:
        loading = True
        if strand_radio == 1:
            shipstrand = "+"
        else:
            shipstrand = "-"
        if not seq_contents:
            return (
                modal,
                "No fasta file uploaded",
                loading,
            )  # Return the error message if no file
        else:
            # Decode sequence content
            content_type, content_string = seq_contents.split(",")
            seq_decoded = base64.b64decode(content_string).decode("utf-8")

            # Parse FASTA to get sequence
            # Assuming single sequence in FASTA
            sequence = parse_fasta(seq_decoded, seq_filename)[0]["sequence"]

            # Get existing ships for comparison
            existing_ships = fetch_ships(curated=True, with_sequence=True)

            # Assign accession and check if review needed
            accession, needs_review = assign_accession(sequence, existing_ships)

        insert_submission(
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
            accession=accession,
            needs_review=needs_review,
        )

        modal = not is_open if not close_modal else False
        message = html.Div(
            [
                html.H4("Successfully submitted!", className="mb-3"),
                dmc.Text(f"Assigned accession: {accession}", className="text-muted"),
                dmc.Text(
                    f"Review status: {'Needs review' if needs_review else 'Auto-approved'}",
                    className="text-muted",
                ),
                dmc.Text(f"Filename: {seq_filename}", className="text-muted"),
                dmc.Text(f"Uploaded by: {uploader}", className="text-muted"),
            ]
        )

        loading = False

    return modal, message, loading


@callback(
    Output("submit-fasta-sequence-upload", "children"),
    [
        Input("submit-fasta-upload", "contents"),
        Input("submit-fasta-upload", "filename"),
    ],
)
def update_fasta_details(seq_contents, seq_filename):
    if seq_contents is None:
        return html.Div(html.P(["Select a FASTA file to upload"]))

    try:
        content_type, content_string = seq_contents.split(",")
        query_string = base64.b64decode(content_string).decode("utf-8")
        records = parse_fasta(query_string, seq_filename)
        return create_fasta_display(records, seq_filename)
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

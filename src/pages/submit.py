import dash
import dash_bootstrap_components as dbc
from dash import dcc, html, callback
from dash.dependencies import Output, Input
import sqlite3
from dash.exceptions import PreventUpdate

dash.register_page(__name__)

# Create SQLite database connection
conn = sqlite3.connect("project-vol/starbase.sqlite")
c = conn.cursor()

# Create table for storing form submissions
# c.execute('''CREATE TABLE submissions (
# 	fna	TEXT NOT NULL,
# 	gff3	TEXT,
# 	genus	TEXT NOT NULL,
# 	species	INTEGER NOT NULL,
# 	evidence	TEXT,
# 	uploader	TEXT NOT NULL,
# 	comment	TEXT,
# 	id	INTEGER NOT NULL,
# 	PRIMARY KEY(id AUTOINCREMENT)
# )''')
# conn.commit()

# Define form layout
form = html.Div(
    [
        html.H1("Submit Starships to starbase"),
        dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            "Submission of multiple Starships to starbase"
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.Div(
                                                    [
                                                        html.P(
                                                            [
                                                                "Unfortunately, we can only handle submission for one Starship at a time. If you have a batch of Starships that you'd like to submit, please send the submission via ",
                                                                html.A(
                                                                    "email",
                                                                    href="mailto:adrian.e.forsythe@gmail.com",
                                                                ),
                                                            ]
                                                        ),
                                                    ]
                                                )
                                            ]
                                        ),
                                    ]
                                )
                            ],
                            width=8,
                        )
                    ],
                    className="mb-3",
                ),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            "Submit individual Starships to starbase"
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.H4("File Upload:"),
                                                dcc.Upload(
                                                    id="uploadsequence",
                                                    children="Upload ship sequence",
                                                    accept=".fa, .fna, .fasta",
                                                    multiple=False,
                                                    style={
                                                        "width": "75%",
                                                        "height": "60px",
                                                        "lineHeight": "60px",
                                                        "borderWidth": "1px",
                                                        "borderStyle": "dashed",
                                                        "borderRadius": "5px",
                                                        "textAlign": "center",
                                                        "margin": "10px",
                                                    },
                                                ),
                                                dcc.Upload(
                                                    id="uploadannotations",
                                                    children="Upload gene annotations associated with Starship sequence (GFF[3] or BED format)",
                                                    accept=".gff, .gff3, .bed",
                                                    multiple=False,
                                                    style={
                                                        "width": "75%",
                                                        "height": "60px",
                                                        "lineHeight": "60px",
                                                        "borderWidth": "1px",
                                                        "borderStyle": "dashed",
                                                        "borderRadius": "5px",
                                                        "textAlign": "center",
                                                        "margin": "10px",
                                                    },
                                                ),
                                                html.Br(),
                                                html.H4("Starship metadata:"),
                                                dcc.Input(
                                                    id="uploader",
                                                    type="text",
                                                    placeholder="Name of curator",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                dcc.Input(
                                                    id="evidence",
                                                    type="text",
                                                    placeholder="How were Starships annotated? (i.e. starfish)",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                dcc.Input(
                                                    id="genus",
                                                    type="text",
                                                    placeholder="Enter genus name",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                dcc.Input(
                                                    id="species",
                                                    type="text",
                                                    placeholder="Enter species name",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                html.Br(),
                                                html.Br(),
                                                html.H4(
                                                    "Coordinates of Starship in host genome:"
                                                ),
                                                dcc.Input(
                                                    id="hostchr",
                                                    type="text",
                                                    placeholder="Host genome contig/scaffold/chromosome ID",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                dcc.Input(
                                                    id="shipstart",
                                                    type="text",
                                                    placeholder="Start coordinate of Starship",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                dcc.Input(
                                                    id="shipend",
                                                    type="text",
                                                    placeholder="End coordinate for Starship",
                                                    style={"width": "75%"},
                                                    required=True,
                                                ),
                                                html.Br(),
                                                html.Br(),
                                                html.H4("Additional information:"),
                                                dcc.Textarea(
                                                    id="comment",
                                                    placeholder="Any comments about the Starship features, annotations, or host genome?",
                                                    style={
                                                        "height": "100px",
                                                        "width": "75%",
                                                    },
                                                    required=False,
                                                ),
                                                html.Br(),
                                                html.Br(),
                                                html.Button(
                                                    "Submit",
                                                    id="submitship",
                                                    n_clicks=0,
                                                ),
                                                html.Div(id="outputstate"),
                                            ]
                                        ),
                                    ]
                                )
                            ]
                        )
                    ]
                ),
            ],
            className="mt-4",
        ),
    ],
)

# Define app layout
layout = html.Div([dbc.Container(form, className="mt-4")])


# Callback to handle form submission
@callback(
    Output("outputstate", "children"),
    [Input("submitship", "n_clicks")],
    [Input("uploadsequence", "uploadsequence")],
    [Input("uploadannotations", "uploadannotations")],
    [Input("uploader", "uploader")],
    [Input("evidence", "evidence")],
    [Input("genus", "genus")],
    [Input("species", "species")],
    [Input("hostchr", "hostchr")],
    [Input("shipstart", "shipstart")],
    [Input("shipend", "shipend")],
    [Input("comment", "comment")],
)
def update_output(
    n_clicks,
    uploadsequence,
    uploadannotations,
    uploader,
    evidence,
    genus,
    species,
    hostchr,
    shipstart,
    shipend,
    comment,
):
    if n_clicks == 0:
        raise PreventUpdate

    if uploadsequence is None or len(uploadsequence) == 0:
        return "No fasta file uploaded"

    if uploadannotations is None or len(uploadannotations) == 0:
        return "No gff file uploaded"

    # Save file and form data to SQLite database
    for filename in uploadsequence:
        file_path = f"uploads/{filename}"  # Assuming uploads folder exists
        c.execute(
            "INSERT INTO submissions (upload-sequence, upload-annotations, uploader, evidence, genus, species, host-chr, ship-start, ship-end, comment) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (
                uploadsequence,
                uploadannotations,
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

    return "Successfully submitted to starbase"

import warnings

import dash
from dash import dcc, html
import dash_bootstrap_components as dbc

warnings.filterwarnings("ignore")

dash.register_page(__name__)

layout = html.Div(
    [
        html.Div(
            [
                html.H1(
                    [
                        "Annotate Starships with ",
                        html.Span(
                            "starfish",
                            className="logo-text",
                        ),
                    ]
                ),
                html.P(
                    [
                        html.Span(
                            "starfish",
                            className="logo-text",
                        ),
                        " is a modular toolkit for giant mobile element annotation. Built primarily for annotating Starship elements in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal Starship or a bacterial integrative and conjugative element: a 'captain' gene with zero or more 'cargo genes downstream of its 3' end. It is particularly well suited for annotating low-copy number elements in a content independent manner.",
                    ]
                ),
                html.Img(
                    src="assets/logos/STARFISH_LOGO.png",
                    width="50%",
                ),
            ],
            style={"padding": "20px"},
        ),
        html.Div(
            [
                html.H1(
                    [
                        "Run ",
                        html.Span(
                            "starfish",
                            className="logo-text",
                        ),
                        " on your genome",
                    ]
                ),
                dcc.Upload(
                    id="input_fasta",
                    children=html.Div(["Drag and Drop or ", html.A("Select Files")]),
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
                    multiple=False,
                ),
                dcc.Upload(
                    id="input_gff",
                    children=html.Div(["Drag and Drop or ", html.A("Select Files")]),
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
                    multiple=False,
                ),
                dbc.Button(
                    [
                        "Launch ",
                        html.Span(
                            "starfish",
                            className="logo-text",
                        ),
                    ],
                    id="starfish",
                    n_clicks=0,
                    style={"margin": "10px"},
                ),
            ],
            style={"padding": "20px"},
        ),
    ],
    style={
        "width": "40%",
        "margin": "20px",
        "border": "1px solid #ddd",
        "borderRadius": "5px",
    },
)

# # Define server logic
# @callback(
#     Output("output-data-upload", "children"),
#     Input("starfish", "n_clicks"),
#     State("input_fasta", "contents"),
#     State("input_gff", "contents"),
# )
# def launch_starfish(n_clicks, fasta_content, gff_content):
#     if n_clicks == 0:
#         raise dash.exceptions.PreventUpdate
#     else:
#         if fasta_content is None or gff_content is None:
#             return html.Div(
#                 "Provide both genome and annotations in the required format."
#             )

#         # Process the files
#         # Example: Convert base64 string to bytes and save to file
#         fasta_bytes = base64.b64decode(fasta_content.split(",")[1])
#         with open("genome.fasta", "wb") as f:
#             f.write(fasta_bytes)

#         gff_bytes = base64.b64decode(gff_content.split(",")[1])
#         with open("annotations.gff", "wb") as f:
#             f.write(gff_bytes)

#         # Execute starfish script or function here

#         return html.Div("Starfish launched successfully!")

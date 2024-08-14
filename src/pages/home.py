import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from src.components.paperTable import table

dash.register_page(__name__, title="Home", name="Home", path="/")

working = {
    "wiki": "Catalogue/Wiki of Starship Metadata",
    "submit": "Submission of new Starship sequences",
    "blast": "BLAST/HMMER searches",
}
working_buttons = [
    dbc.Button(
        value,
        href=f"/{key}",
        external_link=False,
        color="primary",
        className="mx-auto",
    )
    for key, value in working.items()
]

not_working = [
    "Synteny/Genome Browser",
    html.P(
        [
            html.Span(
                "starfish",
                className="logo-text",
            ),
            " webserver",
        ],
    ),
]
not_working_ul = html.Ul(
    [html.Li(item) for item in not_working],
)

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            style={"paddingTop": "20px", "paddingBottom": "20px"},
            children=[
                dbc.Col(
                    lg=8,
                    sm=12,
                    children=[
                        dbc.Stack(
                            [
                                html.H1(
                                    [
                                        html.Span(
                                            "starbase: ",
                                            className="logo-text",
                                        ),
                                        "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                                    ],
                                    className="text-center",
                                ),
                                dbc.Stack(
                                    [
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H4(
                                                        [
                                                            "What can I currently use ",
                                                            html.Span(
                                                                "starbase",
                                                                className="logo-text",
                                                            ),
                                                            " for?",
                                                        ],
                                                    )
                                                ),
                                                dbc.CardBody(
                                                    dbc.Stack(
                                                        working_buttons,
                                                        direction="horizontal",
                                                        gap=3,
                                                        className="justify-content-center",
                                                    )
                                                ),
                                            ],
                                            className="w-100 mb-3",
                                            style={
                                                "minHeight": "150px",
                                            },
                                        ),
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H4(
                                                        [
                                                            "Functions of ",
                                                            html.Span(
                                                                "starbase",
                                                                className="logo-text",
                                                            ),
                                                            " under active development:",
                                                        ],
                                                    )
                                                ),
                                                dbc.CardBody([not_working_ul]),
                                            ],
                                            className="w-100",
                                            style={
                                                "minHeight": "150px",
                                            },
                                        ),
                                    ],
                                    direction="vertical",
                                    gap=3,
                                ),
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            html.H4(
                                                ["What is a Starship?"],
                                            ),
                                        ),
                                        dbc.CardBody(
                                            [
                                                dbc.Row(
                                                    [
                                                        dbc.Col(
                                                            lg=6,
                                                            sm=12,
                                                            children=[
                                                                html.P(
                                                                    [
                                                                        "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes). They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in ",
                                                                        html.Span(
                                                                            "Paecilomyces",
                                                                            style={
                                                                                "font-style": "italic"
                                                                            },
                                                                        ),
                                                                        " cheese making in ",
                                                                        html.Span(
                                                                            "Penicillium",
                                                                            style={
                                                                                "font-style": "italic",
                                                                            },
                                                                        ),
                                                                        ", and enable the transfer of formaldehyde resistance in ",
                                                                        html.Span(
                                                                            "Aspergillus nidulans",
                                                                            style={
                                                                                "font-style": "italic",
                                                                            },
                                                                        ),
                                                                        " and ",
                                                                        html.Span(
                                                                            "Penicillium chrysogenum.",
                                                                            style={
                                                                                "font-style": "italic",
                                                                            },
                                                                        ),
                                                                    ],
                                                                ),
                                                            ],
                                                        ),
                                                        dbc.Col(
                                                            lg=6,
                                                            sm=12,
                                                            children=[
                                                                dbc.CardImg(
                                                                    src="assets/images/starship-model.png",
                                                                    style={
                                                                        "backgroundColor": "white",
                                                                        "margin": "auto",
                                                                        "display": "block",
                                                                        "height": "100%",
                                                                        "object-fit": "cover",
                                                                    },
                                                                    className="text-center",
                                                                )
                                                            ],
                                                        ),
                                                    ]
                                                ),
                                            ]
                                        ),
                                    ],
                                    color="primary",
                                    inverse=True,
                                    className="mt-3",
                                    style={
                                        "width": "100%",
                                    },
                                ),
                                table,
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [html.H4("Data Availability")],
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.P(
                                                    [
                                                        "We have been maintaining ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                    ],
                                                ),
                                                html.Div(
                                                    style={"textAlign": "center"},
                                                    children=[
                                                        dbc.Button(
                                                            html.P(
                                                                [
                                                                    "Download the latest version of ",
                                                                    html.Span(
                                                                        "starbase",
                                                                        className="logo-text",
                                                                    ),
                                                                    ".",
                                                                ],
                                                            ),
                                                            id="dl-button",
                                                            color="primary",
                                                            className="d-grid gap-2 col-3 mx-auto",
                                                        ),
                                                        dcc.Download(id="dl-package"),
                                                    ],
                                                ),
                                            ]
                                        ),
                                    ],
                                    className="mt-3",
                                ),
                            ],
                            gap=4,
                            direction="vertical",
                        )
                    ],
                ),
            ],
        ),
    ],
)

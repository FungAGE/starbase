import dash
import dash_bootstrap_components as dbc
from dash import html

dash.register_page(__name__, title="Home", name="Home", path="/")


def mod_home_ui():
    working = [
        "Catalogue/Wiki of Starship Metadata",
        "Submission of new Starship sequences",
        "BLAST/HMMER searches",
    ]
    working_ul = html.Ul(
        [html.Li(item) for item in working],
        style={
            "fontSize": "0.6vw",
        },
    )

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
        style={
            "fontSize": "0.6vw",
        },
    )

    return html.Div(
        [
            dbc.Container(
                fluid=True,
                children=[
                    dbc.Col(
                        width={"size": 6, "offset": 3},
                        children=[
                            dbc.Row(
                                html.H1(
                                    [
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                            style={
                                                "fontSize": "7vw",
                                            },
                                        ),
                                        html.Br(),
                                        "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                                    ],
                                    style={
                                        "fontSize": "2vw",
                                        "textAlign": "center",
                                    },
                                )
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
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
                                                            style={
                                                                "fontSize": "1vw",
                                                            },
                                                        )
                                                    ),
                                                    dbc.CardBody([working_ul]),
                                                ]
                                            )
                                        ]
                                    ),
                                    dbc.Col(
                                        [
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
                                                            style={
                                                                "fontSize": "1vw",
                                                            },
                                                        )
                                                    ),
                                                    dbc.CardBody([not_working_ul]),
                                                ]
                                            )
                                        ]
                                    ),
                                ]
                            ),
                            dbc.Row(
                                dbc.Col(
                                    html.Img(
                                        src="assets/images/starbase-map.png",
                                        width="100%",
                                    )
                                )
                            ),
                            dbc.Row(
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [html.H4("Data Availability")],
                                            style={"fontSize": "2vw"},
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
                                                    style={"fontSize": "0.6vw"},
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
                                                                style={
                                                                    "fontSize": "0.6vw"
                                                                },
                                                            ),
                                                            id="dl_package",
                                                            color="primary",
                                                            class_name="mr-1",
                                                        )
                                                    ],
                                                ),
                                            ]
                                        ),
                                    ]
                                )
                            ),
                        ],
                        md=5,
                        align="center",
                    ),
                ],
            ),
        ]
    )


layout = mod_home_ui

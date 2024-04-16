import dash
import dash_bootstrap_components as dbc
from dash import html

dash.register_page(__name__, title="Home", name="Home", path="/")


def mod_home_ui():
    working = [
        "Catalogue/Wiki of Starship Metadata",
        "Submission of new Starship sequences",
    ]
    working_ul = html.Ul([html.Li(item) for item in working])

    not_working = [
        "BLAST/HMMER searches",
        "Synteny/Genome Browser",
        "starfish webserver",
    ]
    not_working_ul = html.Ul([html.Li(item) for item in not_working])

    return html.Div(
        [
            html.Div(
                style={
                    "display": "flex",
                    "justify-content": "center",
                    "align-items": "center",
                },
                children=[
                    html.Table(
                        style={"width": "50%"},
                        children=[
                            html.Tr(
                                [
                                    html.Td(
                                        style={"width": "50%"},
                                        children=[
                                            html.H1(
                                                [
                                                    html.Span(
                                                        "starbase",
                                                        className="logo-text",
                                                    ),
                                                    ": A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                                                ]
                                            )
                                        ],
                                    )
                                ]
                            ),
                            html.Br(),
                            html.Tr(
                                [
                                    html.Td(
                                        style={"width": "50%"},
                                        children=[
                                            html.Div(
                                                [
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
                                                                                    ]
                                                                                )
                                                                            ),
                                                                            dbc.CardBody(
                                                                                [
                                                                                    working_ul
                                                                                ]
                                                                            ),
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
                                                                                    ]
                                                                                )
                                                                            ),
                                                                            dbc.CardBody(
                                                                                [
                                                                                    not_working_ul
                                                                                ]
                                                                            ),
                                                                        ]
                                                                    )
                                                                ]
                                                            ),
                                                        ]
                                                    )
                                                ]
                                            )
                                        ],
                                    )
                                ]
                            ),
                            html.Tr(
                                [
                                    html.Td(
                                        style={"width": "100%", "text-align": "center"},
                                        children=[
                                            html.Img(
                                                src="assets/images/starbase-map.png",
                                                width="100%",
                                            )
                                        ],
                                    )
                                ]
                            ),
                        ],
                    )
                ],
            ),
            html.Br(),
            html.Div(
                style={
                    "display": "flex",
                    "justify-content": "center",
                    "align-items": "center",
                },
                children=[
                    html.Table(
                        style={"width": "50%"},
                        children=[
                            html.Tr(
                                [
                                    html.Td(
                                        style={"width": "50%"},
                                        children=[
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader(
                                                        [html.H4("Data Availability")]
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
                                                                ]
                                                            ),
                                                            html.Div(
                                                                style={
                                                                    "textAlign": "center"
                                                                },
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
                                                                            ]
                                                                        ),
                                                                        id="dl_package",
                                                                        color="primary",
                                                                        className="mr-1",
                                                                    )
                                                                ],
                                                            ),
                                                        ]
                                                    ),
                                                ]
                                            ),
                                        ],
                                    )
                                ]
                            )
                        ],
                    )
                ],
            ),
        ]
    )


layout = mod_home_ui

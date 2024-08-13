import dash
import dash_bootstrap_components as dbc
from dash import dcc, html

dash.register_page(__name__)

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            children=[
                dbc.Stack(
                    [
                        dbc.Row(
                            justify="center",
                            align="center",
                            children=[
                                dbc.Col(
                                    lg=8,
                                    sm=12,
                                    children=[
                                        html.P(
                                            [
                                                html.Span(
                                                    "starbase",
                                                    className="logo-text",
                                                ),
                                                " was developed by the ",
                                                html.A(
                                                    "FungAGE lab",
                                                    href="https://fungage.github.io/",
                                                ),
                                                " in collaboration with the Gluck-Thaler lab.",
                                            ],
                                            className="text-center",
                                            style={"font-size": "1rem"},
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
                                    align="center",
                                    children=[
                                        dbc.Stack(
                                            [
                                                dbc.Card(
                                                    [
                                                        dbc.Row(
                                                            [
                                                                dbc.Col(
                                                                    dbc.CardImg(
                                                                        src="assets/images/aaron.png",
                                                                        className="img-fluid rounded-start",
                                                                    ),
                                                                    width={
                                                                        "size": 4,
                                                                        "order": 1,
                                                                    },
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                "Aaron Vogan",
                                                                                className="card-title",
                                                                            ),
                                                                            html.P(
                                                                                "FungAGE group leader",
                                                                                className="card-text",
                                                                            ),
                                                                            dbc.Button(
                                                                                html.I(
                                                                                    className="bi bi-envelope"
                                                                                ),
                                                                                href="mailto:aaron.vogan@ebc.uu.se",
                                                                                color="teal",
                                                                                size="lg",
                                                                                className="me-2",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    width={
                                                                        "size": 8,
                                                                        "order": 2,
                                                                    },
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={
                                                        "maxWidth": "100%",
                                                    },
                                                ),
                                                dbc.Card(
                                                    [
                                                        dbc.Row(
                                                            [
                                                                dbc.Col(
                                                                    dbc.CardImg(
                                                                        src="assets/images/adrian.png",
                                                                        className="img-fluid rounded-start",
                                                                    ),
                                                                    width={
                                                                        "size": 4,
                                                                        "order": 1,
                                                                    },
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                "Adrian Forsythe",
                                                                                className="card-title",
                                                                            ),
                                                                            html.P(
                                                                                [
                                                                                    html.Span(
                                                                                        "starbase",
                                                                                        className="logo-text",
                                                                                    ),
                                                                                    " lead developer",
                                                                                ],
                                                                                className="card-text",
                                                                            ),
                                                                            dbc.Button(
                                                                                html.I(
                                                                                    className="bi bi-envelope"
                                                                                ),
                                                                                href="mailto:adrian.e.forsythe@gmail.com",
                                                                                color="teal",
                                                                                size="lg",
                                                                                className="me-2",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    width={
                                                                        "size": 8,
                                                                        "order": 2,
                                                                    },
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={
                                                        "maxWidth": "100%",
                                                    },
                                                ),
                                                dbc.Card(
                                                    [
                                                        dbc.Row(
                                                            [
                                                                dbc.Col(
                                                                    dbc.CardImg(
                                                                        src="assets/images/emile.png",
                                                                        className="img-fluid rounded-start",
                                                                    ),
                                                                    width={
                                                                        "size": 4,
                                                                        "order": 1,
                                                                    },
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                "Emile Gluck-Thaler",
                                                                                className="card-title",
                                                                            ),
                                                                            html.P(
                                                                                [
                                                                                    html.Span(
                                                                                        "starfish",
                                                                                        className="logo-text",
                                                                                    ),
                                                                                    " lead developer, Gluck-Thaler lab group leader",
                                                                                ],
                                                                                className="card-text",
                                                                            ),
                                                                            dbc.Button(
                                                                                html.I(
                                                                                    className="bi bi-envelope"
                                                                                ),
                                                                                href="mailto:emilegluckthaler@gmail.com",
                                                                                color="teal",
                                                                                size="lg",
                                                                                className="me-2",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    width={
                                                                        "size": 8,
                                                                        "order": 2,
                                                                    },
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={
                                                        "maxWidth": "100%",
                                                    },
                                                ),
                                            ],
                                            gap=3,
                                        )
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
                                        html.P(
                                            [
                                                "The source code for ",
                                                html.Span(
                                                    "starbase",
                                                    className="logo-text",
                                                ),
                                                " webserver will soon be available on GitHub",
                                            ],
                                            className="text-center",
                                            style={"font-size": "1rem"},
                                        ),
                                        html.Div(
                                            html.Img(
                                                src="assets/images/starbase-map.png",
                                                className="img-fluid",
                                            )
                                        ),
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H2(
                                                        "Data Availability",
                                                        className="text-center",
                                                    )
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
                                                                " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the meantime, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                            ],
                                                            style={
                                                                "font-size": "0.875rem"
                                                            },
                                                        ),
                                                        html.Div(
                                                            [
                                                                dbc.Button(
                                                                    "Download the latest version of starbase.",
                                                                    id="dl-button",
                                                                    color="primary",
                                                                    className="mt-2",
                                                                ),
                                                                dcc.Download(
                                                                    id="dl-package"
                                                                ),
                                                            ],
                                                            className="text-center",
                                                            style={
                                                                "font-size": "0.875rem",
                                                            },
                                                        ),
                                                    ]
                                                ),
                                            ],
                                            className="mt-4",
                                        ),
                                    ],
                                )
                            ],
                        ),
                    ],
                    gap=4,
                    direction="vertical",
                )
            ],
            style={"padding": "20px"},
        ),
    ],
)

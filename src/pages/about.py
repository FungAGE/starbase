import dash
import dash_bootstrap_components as dbc
from dash import html


dash.register_page(__name__)

layout = html.Div(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.Stack(
                    [
                        dbc.Row(),
                        dbc.Row(
                            justify="center",
                            align="center",
                            children=[
                                dbc.Col(
                                    lg=8,
                                    sm=12,
                                    className="align-self-center",
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
                                                " in collaboration with the Gluck-Thaler lab. The sourcecode for ",
                                                html.Span(
                                                    "starbase",
                                                    className="logo-text",
                                                ),
                                                " webserver will soon be available on GitHub",
                                            ],
                                            className="mr-1",
                                        )
                                    ],
                                ),
                            ],
                        ),
                        dbc.Row(
                            justify="center",
                            align="center",
                            children=[
                                dbc.Col(
                                    lg=8,
                                    sm=12,
                                    className="align-self-center",
                                    children=[
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H2(
                                                        ["Data Availability"],
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
                                                                " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                            ],
                                                            className="mr-1",
                                                        ),
                                                        html.Div(
                                                            [
                                                                dbc.Button(
                                                                    html.P(
                                                                        [
                                                                            "Download the latest version of ",
                                                                            html.Span(
                                                                                "starbase",
                                                                                className="logo-text",
                                                                            ),
                                                                        ],
                                                                    ),
                                                                    id="dl_package",
                                                                    color="primary",
                                                                    className="mr-1",
                                                                ),
                                                            ],
                                                            style={
                                                                "textAlign": "center"
                                                            },
                                                        ),
                                                    ]
                                                ),
                                            ],
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
                                    className="align-self-center",
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
                                                                    className="col-md-4",
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                [
                                                                                    "Aaron Vogan",
                                                                                ],
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
                                                                                className="mr-1",
                                                                                size="lg",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    className="col-md-8",
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={"maxWidth": "540px"},
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
                                                                    className="col-md-4",
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                [
                                                                                    "Adrian Forsythe",
                                                                                ],
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
                                                                                className="mr-1",
                                                                                size="lg",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    className="col-md-8",
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={"maxWidth": "540px"},
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
                                                                    className="col-md-4",
                                                                ),
                                                                dbc.Col(
                                                                    dbc.CardBody(
                                                                        [
                                                                            html.H3(
                                                                                [
                                                                                    "Emile Gluck-Thaler",
                                                                                ],
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
                                                                                className="mr-1",
                                                                                size="lg",
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    className="col-md-8",
                                                                ),
                                                            ],
                                                            className="g-0 d-flex align-items-center",
                                                        )
                                                    ],
                                                    className="mb-3",
                                                    style={"maxWidth": "540px"},
                                                ),
                                            ],
                                            gap=3,
                                        )
                                    ],
                                )
                            ],
                        ),
                    ],
                    gap=3,
                    direction="vertical",
                )
            ],
        ),
    ],
)

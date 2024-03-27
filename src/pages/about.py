import dash
import dash_bootstrap_components as dbc
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State

import dash_ag_grid as dag

dash.register_page(__name__)

layout = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader(html.H5("About starbase")),
                    dbc.CardBody([
                        html.Div([html.P(["starbase was developed by the ",html.A("FungAGE lab", href="https://fungage.github.io/"), " and the code for the starbase webserver will soon be available on our GitHub"])]),
                        dbc.Card([
                            dbc.CardHeader(html.H6("Data Availability")),
                            dbc.CardBody([
                                html.P("We have been maintaining starbase data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)"),
                                dbc.Button("Download the latest version of starbase.", id="dl_package", color="primary", className="mr-1")
                            ])
                        ])
                    ])
                ])
            ])
        ]),
        html.Br(),
        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.Row([
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
                                        html.H4("Aaron Vogan", className="card-title"),
                                        html.P("FungAGE group leader",
                                            className="card-text",
                                        ),
                                        dbc.Button(html.I(className="bi bi-envelope"), href="mailto:aaron.vogan@ebc.uu.se", color="teal", className="mr-1")
                                    ]
                                ),
                                className="col-md-8",
                            ),
                        ],
                        className="g-0 d-flex align-items-center",
                    )
                ],
                className="mb-3",
                style={"maxWidth": "540px"})
            ])
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Card([
                        dbc.Row([
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
                                            html.H4("Adrian Forsythe", className="card-title"),
                                            html.P("starbase lead developer",
                                                className="card-text",
                                            ),
                                            dbc.Button(html.I(className="bi bi-envelope"), href="mailto:adrian.e.forsythe@gmail.com", color="teal", className="mr-1")
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
                )
            ])
        ]),
        dbc.Row([
            dbc.Col([        
                dbc.Card([
                        dbc.Row([
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
                                            html.H4("Emile Gluck-Thaler", className="card-title"),
                                            html.P("starfish lead developer",
                                                className="card-text",
                                            ),
                                            dbc.Button(html.I(className="bi bi-envelope"), href="mailto:emilegluckthaler@gmail.com", color="teal", className="mr-1")
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
                )
            ])
        ])
    ])
])
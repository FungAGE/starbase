import dash_bootstrap_components as dbc
from dash import html

from src.pages import (
    HOME_URL,
    WIKI_URL,
    EXPLORE_URL,
    BLAST_URL,
    # MUMMER_URL,
    SUBMIT_URL,
    ABOUT_URL,
)

navbar = dbc.Navbar(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                dbc.Collapse(
                    [
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(
                                    ["Home"],
                                    className="nav-item-text",
                                ),
                                href=HOME_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(
                                    [
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                        " Wiki",
                                    ],
                                    className="nav-item-text",
                                ),
                                href=WIKI_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(
                                    [
                                        "Explore ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
                                    className="nav-item-text",
                                ),
                                href=EXPLORE_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(
                                    [
                                        "BLAST/hmmersearch ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
                                    className="nav-item-text",
                                ),
                                href=BLAST_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                        # dbc.NavItem(
                        #     dbc.NavLink(
                        #         html.P(
                        #             [
                        #                 "Mummer",
                        #             ],
                        #             className="nav-item-text",
                        #         ),
                        #         href=MUMMER_URL,
                        #         active="exact",
                        #         className="nav-item-link",
                        #     )
                        # ),
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(
                                    [
                                        "Submit to ",
                                        html.Span(
                                            "starbase",
                                            className="logo-text",
                                        ),
                                    ],
                                    className="nav-item-text",
                                ),
                                href=SUBMIT_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                        dbc.NavItem(
                            dbc.NavLink(
                                html.P(["About"], className="nav-item-text"),
                                href=ABOUT_URL,
                                active="exact",
                                className="nav-item-link",
                            )
                        ),
                    ],
                    id="navbar-collapse",
                    is_open=False,
                    navbar=True,
                    className="justify-content-start",
                ),
            ],
        ),
    ]
)


def sidebar():
    return html.Div(
        [
            html.Img(
                src="assets/logos/favicon.svg",
                style={"height": "12rem", "width": "12rem"},
            ),
            html.Hr(),
            dbc.Nav(
                [
                    dbc.NavLink("Home", href=HOME_URL, active="exact"),
                    dbc.NavLink("Wiki", href=WIKI_URL, active="exact"),
                    dbc.NavLink(
                        html.P(
                            ["Explore ", html.Span("starbase", className="logo-text")]
                        ),
                        href=EXPLORE_URL,
                        active="exact",
                    ),
                    dbc.NavLink("BLAST", href=BLAST_URL, active="exact"),
                    # dbc.NavLink("Genome Browser", href=IGV_URL, active="exact"),
                    # dbc.NavLink("Starfish", href=STARFISH_URL, active="exact"),
                    dbc.NavLink(
                        html.P(
                            ["Submit to ", html.Span("starbase", className="logo-text")]
                        ),
                        href=SUBMIT_URL,
                        active="exact",
                    ),
                    dbc.NavLink("About", href=ABOUT_URL, active="exact"),
                ],
                vertical=True,
                pills=True,
            ),
        ],
        style=SIDEBAR_STYLE,
    )


# def test_blast_ui():
#     html.Iframe(src='assets/blaster.html')
# layout = test_blast_ui

# def blasterjs_ui():
#     html.Div([
#         html.Div([
#             # dcc.Upload(
#             #     id='blastinput',
#             #     children=html.Div([
#             #         'Drag and Drop or ',
#             #         html.A('Select Files')
#             #     ]),
#             #     multiple=False
#             # ),
#             html.Div(id='blast-multiple-alignments'),
#             html.Div(id='blast-alignments-table'),
#             html.Div(id='blast-single-alignment')
#         ]),
#         html.Script(src='lib/html2canvas.min.js', type='text/javascript'),
#         html.Script(src='lib/blaster.min.js', type='text/javascript'),
#         html.Script("""
#             var blasterjs = require("biojs-vis-blasterjs")
#             var instance = new blasterjs({
#                 string: "assets/blast.out",
#                 multipleAlignments: "blast-multiple-alignments",
#                 alignmentsTable: "blast-alignments-table",
#                 singleAlignment: "blast-single-alignment",
#             })
#         """, type='text/javascript')
#     ])

# layout = blasterjs_ui

import dash_bootstrap_components as dbc
from dash import html
from dash import callback
from dash.dependencies import Output, Input, State

from src.pages import (
    HOME_URL,
    WIKI_URL,
    EXPLORE_URL,
    BLAST_URL,
    # IGV_URL,
    SUBMIT_URL,
    ABOUT_URL,
)

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

navbar = dbc.Navbar(
    dbc.Container(
        [
            dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
            dbc.Collapse(
                dbc.Row(
                    align="center",
                    className="g-0",
                    children=[
                        dbc.Col(
                            width={"size": 12, "offset": 2},
                            children=[
                                dbc.NavbarSimple(
                                    fluid=True,
                                    children=[
                                        dbc.NavItem(
                                            dbc.NavLink(
                                                "Home", href=HOME_URL, active="exact"
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
                                                    ]
                                                ),
                                                href=WIKI_URL,
                                                active="exact",
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
                                                    ]
                                                ),
                                                href=EXPLORE_URL,
                                                active="exact",
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
                                                    ]
                                                ),
                                                href=BLAST_URL,
                                                active="exact",
                                            )
                                        ),
                                        # dbc.NavItem(dbc.NavLink("Genome Browser", href=IGV_URL, active="exact")),
                                        # dbc.NavItem(dbc.NavLink("Starfish", href=STARFISH_URL, active="exact")),
                                        dbc.NavItem(
                                            dbc.NavLink(
                                                html.P(
                                                    [
                                                        "Submit to ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                    ]
                                                ),
                                                href=SUBMIT_URL,
                                                active="exact",
                                            )
                                        ),
                                        dbc.NavItem(
                                            dbc.NavLink(
                                                "About", href=ABOUT_URL, active="exact"
                                            )
                                        ),
                                    ],
                                    # vertical=True,
                                    # pills=True,
                                    style={
                                        "fontSize": "1vw",
                                    },
                                )
                            ],
                        )
                    ],
                ),
                id="navbar-collapse",
                is_open=False,
                navbar=True,
            ),
        ]
    ),
)


def sidebar():
    return html.Div(
        [navbar],
        # style=SIDEBAR_STYLE,
    )


# add callback for toggling the collapse on small screens
@callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open


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

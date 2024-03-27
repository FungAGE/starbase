import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash import dcc, html, callback
from flask import session

from src.pages import HOME_URL, WIKI_URL, EXPLORE_URL, BLAST_URL, ABOUT_URL, IGV_URL, SUBMIT_URL
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

def sidebar():
    return html.Div([
        html.Img(src="assets/logos/favicon.svg", style={'height': '12rem', 'width': '12rem'}),
        html.Hr(),
        dbc.Nav(
            [
                dbc.NavLink("Home", href=HOME_URL, active="exact"),
                dbc.NavLink("Wiki", href=WIKI_URL, active="exact"),
                dbc.NavLink("Explore starbase", href=EXPLORE_URL, active="exact"),
                dbc.NavLink("BLAST", href=BLAST_URL, active="exact"),
                dbc.NavLink("Genome Browser", href=IGV_URL, active="exact"),
                dbc.NavLink("Submit to starbase", href=SUBMIT_URL, active="exact"),
                dbc.NavLink("About", href=ABOUT_URL, active="exact")
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

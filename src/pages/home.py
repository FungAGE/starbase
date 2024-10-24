import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input

import pandas as pd

from src.components.mariadb import engine

from src.components.tables import make_paper_table
from src.components.callbacks import download_ships_card

import logging

logger = logging.getLogger(__name__)

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
        class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
    )
    for key, value in working.items()
]

not_working = [
    html.Div(
        ["Synteny/Genome Browser"],
    ),
    html.Div(
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
    [
        html.Li(
            item,
            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
        )
        for item in not_working
    ],
)

layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dmc.Container(
            fluid=True,
            children=[
                dmc.Center(
                    children=[
                        dmc.Title(
                            [
                                html.Span(
                                    "starbase: ",
                                    className="logo-text",
                                ),
                                "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                            ],
                            className="text-center",
                            style={"paddingTop": "20px"},
                            # className="text-center text-custom text-custom-xl",
                        ),
                    ],
                ),
                dmc.Grid(
                    justify="center",
                    align="center",
                    style={"paddingTop": "20px"},
                    gutter="xl",
                    children=[
                        dmc.GridCol(
                            span=12,
                            children=[
                                dmc.Center(
                                    [
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    [
                                                        html.Div(
                                                            ["What is a Starship?"],
                                                            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                        ),
                                                    ],
                                                    className="card-header-custom",
                                                ),
                                                dbc.CardBody(
                                                    [
                                                        dmc.Grid(
                                                            [
                                                                dmc.GridCol(
                                                                    span="content",
                                                                    children=[
                                                                        html.Div(
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
                                                                            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl align-items-center",
                                                                            style={
                                                                                "justify-content": "center"
                                                                            },
                                                                        ),
                                                                    ],
                                                                ),
                                                                dmc.GridCol(
                                                                    span="content",
                                                                    children=[
                                                                        dmc.Image(
                                                                            src="assets/images/starship-model.png",
                                                                            style={
                                                                                "backgroundColor": "white",
                                                                                "maxWidth": "1000px",
                                                                            },
                                                                        )
                                                                    ],
                                                                ),
                                                            ]
                                                        ),
                                                    ],
                                                ),
                                            ],
                                            color="primary",
                                            inverse=True,
                                            className="auto-resize-900",
                                        ),
                                    ]
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span="content",
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [
                                                html.Div(
                                                    [
                                                        "What can I currently use ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                        " for?",
                                                    ],
                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                )
                                            ],
                                            style={
                                                "justify-content": "center",
                                            },
                                            className="card-header-custom",
                                        ),
                                        dbc.CardBody(
                                            [
                                                dbc.Stack(
                                                    working_buttons,
                                                    # direction="horizontal",
                                                    gap=3,
                                                    className="justify-content-center",
                                                )
                                            ],
                                            className="d-flex align-items-center",
                                        ),
                                    ],
                                    className="w-100 mb-3",
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span="content",
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [
                                                html.Div(
                                                    [
                                                        "Functions of ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                        " under active development:",
                                                    ],
                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                )
                                            ],
                                            style={
                                                "justify-content": "center",
                                            },
                                            className="card-header-custom",
                                        ),
                                        dbc.CardBody(
                                            [not_working_ul],
                                            className="d-flex align-items-center",
                                        ),
                                    ],
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span="content",
                            children=[
                                dmc.Center([make_paper_table(engine)]),
                            ],
                        ),
                    ],
                ),
                dmc.Grid(
                    justify="center",
                    align="center",
                    style={"paddingTop": "20px"},
                    grow=True,
                    children=[
                        dmc.GridCol(span=12, children=dmc.Center(download_ships_card))
                    ],
                ),
            ],
        ),
    ]
)

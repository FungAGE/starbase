import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

from dash import html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import os


download_ships_button = dbc.Button(
    html.Div(
        [
            "Download Starships from the latest version of ",
            html.Span(
                "starbase",
                className="logo-text",
            ),
            ".",
        ],
    ),
    href="/download",
    color="primary",
    class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
)

download_ships_card = dbc.Card(
    [
        dbc.CardHeader(
            [
                html.Div(
                    "Data Availability",
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                )
            ],
            className="card-header-custom",
        ),
        dbc.CardBody(
            [
                dmc.Text(
                    [
                        "We have been maintaining ",
                        html.Span(
                            "starbase",
                            className="logo-text",
                        ),
                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export",
                    ],
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                    style={"paddingBottom": "20px"},
                ),
                dmc.Center(download_ships_button),
            ],
        ),
    ],
    className="auto-resize-750",
    # style={"height": "350px"},
    #
)

download_starbase_button = html.Div(
    [
        dbc.Button(
            html.P(
                [
                    "Download Starships from the latest version of ",
                    html.Span(
                        "starbase",
                        className="logo-text",
                    ),
                    ".",
                ]
            ),
            id="dl-button",
            color="primary",
            className="mt-2",
        ),
    ],
    className="text-center",
    style={
        "font-size": "0.875rem",
    },
)


def curated_switch(text, size="normal"):
    if size == "normal":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "justify-content": "center",
        }
    if size == "large":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "transform": "scale(1.5)",
            "justify-content": "center",
        }
    switch = dbc.Row(
        justify="center",
        align="start",
        style={"paddingTop": "20px"},
        children=[
            dbc.Col(
                lg=6,
                sm=8,
                children=[
                    dbc.Switch(
                        id="curated-input",
                        label=text,
                        value=False,
                        style=style,
                    ),
                ],
            )
        ],
    )
    return switch

from dash import html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import logging
import pandas as pd

from src.components.cache_manager import load_from_cache
from src.components.sql_queries import fetch_meta_data

logger = logging.getLogger(__name__)

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


def create_accession_modal(accession):    
    logger.debug(f"Looking for accession: {accession}")
    
    # Load DataFrame from cache
    initial_df = load_from_cache("meta_data")
    if initial_df is None or initial_df.empty:
        initial_df = fetch_meta_data()
        
    # Filter for the specific accession
    modal_data = initial_df[initial_df["accession_tag"] == accession]
    logger.debug(f"Found {len(modal_data)} rows for accession {accession}")
    
    if modal_data.empty:
        return html.Div([
            html.P(f"No data found for accession: {accession}"),
            html.P("This might be because:"),
            html.Ul([
                html.Li("The accession is not in the curated dataset"),
                html.Li("The accession data is missing from the meta_data cache"),
                html.Li(f"The accession format doesn't match ({accession})")
            ])
        ]), f"Accession: {accession}"
        
    # Create modal content using the first row of matched data
    try:
        modal_title = html.H2(f"Ship Accession: {accession}")
        modal_content = html.Div(
            [
                html.Div(
                    [
                        html.Strong("starshipID: "),
                        html.Span(modal_data["starshipID"].iloc[0]),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("curated_status: "),
                        html.Span(modal_data["curated_status"].iloc[0]),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("Number of genomes with ship present:"),
                        html.Span(len(modal_data)),
                    ]
                ),
                html.Hr(),
                # TODO: update genome accessions field
                html.Div(
                    [
                        html.Strong("Assembly Accession(s): "),
                        html.Span(""),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("Genome Source: "),
                        html.Span(modal_data["genomeSource"].iloc[0]),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("Citation: "),
                        html.Span(modal_data["citation"].iloc[0]),
                    ]
                ),
                html.Div(
                    [html.Strong("contigID: "), html.Span(modal_data["contigID"].iloc[0])]
                ),
                html.Div(
                    [
                        html.Strong("elementBegin: "),
                        html.Span(modal_data["elementBegin"].iloc[0]),
                    ]
                ),
                html.Div(
                    [
                        html.Strong("elementEnd: "),
                        html.Span(modal_data["elementEnd"].iloc[0]),
                    ]
                ),
                html.Div([html.Strong("size: "), html.Span(modal_data["size"].iloc[0])]),
                html.Hr(),
                html.Div([html.Strong("order: "), html.Span(modal_data["order"].iloc[0])]),
                html.Div(
                    [html.Strong("family: "), html.Span(modal_data["family"].iloc[0])]
                ),
                html.Div(
                    [html.Strong("species: "), html.Span(modal_data["species"].iloc[0])]
                ),
                html.Div(
                    [html.Strong("strain: "), html.Span(modal_data["strain"].iloc[0])]
                ),
                html.Div(
                    [
                        html.Strong("NCBI Taxonomy ID: "),
                        html.A(
                            (
                                int(modal_data["taxID"].iloc[0])
                                if isinstance(modal_data["taxID"].iloc[0], (float, int))
                                else modal_data["taxID"].iloc[0]
                            ),
                            href=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={int(modal_data['taxID'].iloc[0]) if isinstance(modal_data['taxID'].iloc[0], (float, int)) else modal_data['taxID'].iloc[0]}",
                        ),
                    ]
                ),
            ]
        )
        return modal_content, modal_title
    except Exception as e:
        logger.error(f"Error creating modal content: {str(e)}")
        return html.Div(f"Error creating modal: {str(e)}"), "Error"


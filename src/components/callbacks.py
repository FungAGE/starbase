from dash import html, Output, Input, State, callback, no_update
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import logging
import pandas as pd
import traceback

from src.components.sql_engine import sql_connected
from src.components.sql_manager import load_from_cache
from src.components.sql_manager import fetch_meta_data

logger = logging.getLogger(__name__)


download_ships_button = dmc.Anchor(
    dmc.Button(
        [
            dmc.Group(
                [
                    DashIconify(icon="mdi:download"),
                    dmc.Stack(
                        [
                            dmc.Text(
                                "Download Starships",
                                style={"display": "block"}
                            ),
                            dmc.Text(
                                [
                                    "from the latest version of ",
                                    html.Span(
                                        "starbase",
                                        className="logo-text",
                                    ),
                                ],
                                style={"display": "block"}
                            ),
                        ],
                        gap=0,
                        style={
                            "whiteSpace": "normal",
                            "textAlign": "center",
                            "lineHeight": "1.2",
                        }
                    ),
                ],
                gap="xs",
                style={"flexWrap": "nowrap", "justifyContent": "center"}
            ),
        ],
        id="navigate-to-download-btn",
        variant="gradient",
        gradient={"from": "indigo", "to": "cyan"},
        size="lg",
        radius="md",
        fullWidth=True,
        disabled=not sql_connected,
        styles={
            "root": {
                "minHeight": "auto",
                "height": "auto",
                "whiteSpace": "normal",
                "padding": "1rem",
            },
            "inner": {
                "justifyContent": "center",
            }
        }
    ),
    href="/download",
    style={
        "textDecoration": "none",
        "width": "100%",
        "display": "block",
    }
)

download_ships_card = dmc.Paper([
    dmc.Title("Data Availability",order=2,mb="md"),
    dmc.Text(
                        [
                            "We have been maintaining ",
                            html.Span(
                                "starbase",
                                className="logo-text",
                            ),
                            " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export",
                        ],
                        size="lg", c="dimmed",
                        style={"paddingBottom": "20px"},
                    ),
                    dmc.Center(download_ships_button),
                ], p="xl", radius="md", shadow="sm", withBorder=True
            )

def curated_switch(text="Only search curated Starships", size="sm"):
    """Create a switch component for toggling curated-only searches."""
    return dmc.Switch(
        id="curated-input",
        label=text,
        size=size,
        checked=True
    )

def dereplicated_switch(text="Only search dereplicated Starships", size="sm"):
    """Create a switch component for toggling dereplicated-only searches."""
    return dmc.Switch(
        id="dereplicated-input",
        label=text,
        size=size,
        checked=True
    )

def create_accession_modal(accession):    
    try:
        # Always fetch fresh meta_data
        initial_df = fetch_meta_data()
        
        # Clean the accession format in both the input and the dataframe
        accession = str(accession).strip("[]").split("/")[-1].strip()
        initial_df["accession_tag"] = initial_df["accession_tag"].astype(str).apply(
            lambda x: x.strip("[]").split("/")[-1].strip()
        )
        
        # Find the matching row(s)
        modal_data = initial_df[initial_df["accession_tag"] == accession]
        
        if modal_data.empty:
            return html.Div([
                html.P(f"No data found for accession: {accession}"),
                html.P("Cache status:"),
                html.Ul([
                    html.Li(f"Total records in cache: {len(initial_df)}"),
                    html.Li(f"Sample accessions: {', '.join(initial_df['accession_tag'].head().tolist())}"),
                    html.Li(f"Searched for: {accession}")
                ])
            ]), f"Accession: {accession}"
        
        # Create modal content using the first row of matched data
        modal_title = html.H2(f"Ship Accession: {accession}")
        modal_content = html.Div([
            html.Div([
                html.Strong("starshipID: "),
                html.Span(modal_data["starshipID"].iloc[0]),
            ]),
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
        logger.error(f"Error in create_accession_modal: {str(e)}")
        logger.error(traceback.format_exc())
        raise

def create_modal_callback(table_id, modal_id, content_id, title_id, column_check=None):
    @callback(
        Output(modal_id, "opened"),
        Output(content_id, "children"),
        Output(title_id, "children"),
        Output(table_id, "active_cell"),
        [Input(table_id, "active_cell")],
        [
            State(modal_id, "opened"),
            State(table_id, "data"),
            State(table_id, "derived_virtual_data")
        ],
    )
    def toggle_modal(active_cell, is_open, table_data, filtered_data):
        try:
            if not active_cell:
                return False, no_update, no_update, no_update
                
            # Debug the table data
            data_to_use = filtered_data if filtered_data is not None else table_data
            logger.debug(f"Table accessions: {[row['accession_tag'] for row in data_to_use[:5]]}")
            
            # Debug the cache data
            initial_df = load_from_cache("meta_data")
            if initial_df is not None:
                logger.debug(f"Cache accessions: {initial_df['accession_tag'].head().tolist()}")
            else:
                logger.debug("No data in meta_data cache")
                initial_df = fetch_meta_data()
                logger.debug(f"Freshly fetched accessions: {initial_df['accession_tag'].head().tolist()}")
            
            # Get the row data
            data_to_use = filtered_data if filtered_data is not None else table_data
            row_data = data_to_use[active_cell["row"]]
            
            # Clean and standardize the accession tag
            accession = str(row_data["accession_tag"]).strip("[]").split("/")[-1].strip()
            logger.debug(f"Looking for accession in cache: {accession}")
            logger.debug(f"Available accessions in cache: {initial_df['accession_tag'].unique()[:5]}")
            
            modal_content, modal_title = create_accession_modal(accession)
            return True, modal_content, modal_title, None
            
        except Exception as e:
            logger.error(f"Error in toggle_modal: {str(e)}")
            logger.error(traceback.format_exc())
            # Return a more user-friendly error modal instead of raising
            error_content = html.Div([
                html.P("Error loading modal content"),
                html.P(f"Details: {str(e)}"),
            ])
            return True, error_content, "Error", None


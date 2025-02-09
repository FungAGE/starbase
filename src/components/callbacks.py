from dash import html, Output, Input, State, callback, no_update, dcc, clientside_callback, ClientsideFunction
import dash_mantine_components as dmc
from dash_iconify import DashIconify
from typing import List  

import logging
import traceback
import json

from src.config.cache import cache
from src.database.sql_manager import fetch_meta_data

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
                html.Div([
                    html.Strong("Genome Details:"),
                    dmc.Table(
                        striped=True,
                        highlightOnHover=True,
                        children=[
                            html.Thead(
                                html.Tr([
                                    html.Th("Field"),
                                    *[html.Th(f"{modal_data['assembly_accession'].iloc[i] if 'assembly_accession' in modal_data else f'Genome {i+1}'}")
                                      for i in range(len(modal_data))]
                                ])
                            ),
                            html.Tbody([
                                html.Tr([
                                    html.Td("Genome Source"),
                                    *[html.Td(modal_data["genomeSource"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                                html.Tr([
                                    html.Td("Citation"),
                                    *[html.Td(modal_data["citation"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                                html.Tr([
                                    html.Td("ContigID"),
                                    *[html.Td(modal_data["contigID"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                                html.Tr([
                                    html.Td("Element Begin"),
                                    *[html.Td(modal_data["elementBegin"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                                html.Tr([
                                    html.Td("Element End"),
                                    *[html.Td(modal_data["elementEnd"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                                html.Tr([
                                    html.Td("Size"),
                                    *[html.Td(modal_data["size"].iloc[i])
                                      for i in range(len(modal_data))]
                                ]),
                            ])
                        ]
                    ) if len(modal_data) > 1 else html.Div([
                        html.Div([html.Strong("Assembly Accession: "), 
                                html.Span(modal_data["assembly_accession"].iloc[0] if "assembly_accession" in modal_data else "")]),
                        html.Div([html.Strong("Genome Source: "), 
                                html.Span(modal_data["genomeSource"].iloc[0])]),
                        html.Div([html.Strong("Citation: "), 
                                html.Span(modal_data["citation"].iloc[0])]),
                        html.Div([html.Strong("ContigID: "), 
                                html.Span(modal_data["contigID"].iloc[0])]),
                        html.Div([html.Strong("Element Begin: "), 
                                html.Span(modal_data["elementBegin"].iloc[0])]),
                        html.Div([html.Strong("Element End: "), 
                                html.Span(modal_data["elementEnd"].iloc[0])]),
                        html.Div([html.Strong("Size: "), 
                                html.Span(modal_data["size"].iloc[0])]),
                    ])
                ]),
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
    # First, add the clientside callback to handle cell clicks
    clientside_callback(
        ClientsideFunction(
            namespace='clientside',
            function_name='handleCellClick'
        ),
        Output(f"{table_id}-click-data", "children"),
        Input(table_id, "cellClicked"),
        prevent_initial_call=True
    )

    @callback(
        Output(modal_id, "opened"),
        Output(content_id, "children"),
        Output(title_id, "children"),
        Input(f"{table_id}-click-data", "children"),
        prevent_initial_call=True
    )
    def toggle_modal(click_data):
        try:
            if not click_data:
                return False, no_update, no_update
            
            data = json.loads(click_data)
            if data["column"] == "accession_tag":
                # Clean and standardize the accession tag
                accession = str(data["accession"]).strip("[]").split("/")[-1].strip()
                logger.debug(f"Looking for accession in cache: {accession}")
                
                modal_content, modal_title = create_accession_modal(accession)
                return True, modal_content, modal_title
            
            return False, no_update, no_update
            
        except Exception as e:
            logger.error(f"Error in toggle_modal: {str(e)}")
            logger.error(traceback.format_exc())
            error_content = html.Div([
                html.P("Error loading modal content"),
                html.P(f"Details: {str(e)}"),
            ])
            return True, error_content, "Error"


def create_file_upload(
    upload_id: str,
    output_id: str,
    accept_types: List[str],
    placeholder_text: str = "Drag and drop or click to select a file",
    icon: str = "mdi:file-upload",
    **kwargs
) -> dmc.Stack:
    return dmc.Stack([
        dmc.Center(
            DashIconify(
                icon=icon,
                width=40,
                height=40,
                color="#228be6"
            )
        ),
        dcc.Upload(
            id=upload_id,
            children=dmc.Stack([
                html.Div(
                    id=output_id,
                    children=placeholder_text,
                ),
                dmc.Text(
                    f"Accepted formats: {', '.join(accept_types)}",
                    size="sm",
                    c="dimmed"
                ),
            ], align="center", gap="xs"),
            className="upload-box",
            **kwargs
        ),
        dmc.Progress(
            id=f"{upload_id}-progress",
            value=0,
            animated=True,
            style={"display": "none"}
        ),
    ], gap="md")

def create_feedback_button():
    return dmc.Affix(
        dmc.Anchor(
            dmc.Button(
                children=[
                    dmc.Group([
                        DashIconify(icon="octicon:mark-github-16", width=20),
                        dmc.Stack([
                            "Found an issue?",
                            dmc.Text(
                                "Report it on GitHub",
                                size="xs",
                                c="gray.3"
                            )
                        ], gap=2)
                    ], gap=8, style={"minWidth": "180px"})
                ],
                variant="filled",
                color="gray",
                radius="md",
                size="lg",
                style={
                    "padding": "16px 24px",
                    "boxShadow": "0 2px 10px rgba(0, 0, 0, 0.1)",
                    "cursor": "pointer",
                    "width": "100%",
                    "height": "auto",
                }
            ),
            href="https://github.com/FungAGE/starbase/issues",
            target="_blank",
            style={"textDecoration": "none", "display": "block"}
        ),
        position={"bottom": 30, "right": 30},
        zIndex=1000,
    )
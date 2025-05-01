from dash import html, Output, Input, State, callback, no_update, dcc, callback_context
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
from typing import List

import logging
import traceback
import pandas as pd

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
                            dmc.Text("Download Starships", style={"display": "block"}),
                            dmc.Text(
                                [
                                    "from the latest version of ",
                                    html.Span(
                                        "starbase",
                                        className="logo-text",
                                    ),
                                ],
                                style={"display": "block"},
                            ),
                        ],
                        gap=0,
                        style={
                            "whiteSpace": "normal",
                            "textAlign": "center",
                            "lineHeight": "1.2",
                        },
                    ),
                ],
                gap="xs",
                style={"flexWrap": "nowrap", "justifyContent": "center"},
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
            },
        },
    ),
    href="/download",
    style={
        "textDecoration": "none",
        "width": "100%",
        "display": "block",
    },
)

download_ships_card = dmc.Paper(
    [
        dmc.Title("Data Availability", order=2, mb="md"),
        dmc.Text(
            [
                "We have been maintaining ",
                html.Span(
                    "starbase",
                    className="logo-text",
                ),
                " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export",
            ],
            size="lg",
            c="dimmed",
            style={"paddingBottom": "20px"},
        ),
        dmc.Center(download_ships_button),
    ],
    p="xl",
    radius="md",
    shadow="sm",
    withBorder=True,
    h="100%",
)


def curated_switch(text="Only search curated Starships", size="sm"):
    """Create a switch component for toggling curated-only searches."""
    return dmc.Switch(id="curated-input", label=text, size=size, checked=True)


def dereplicated_switch(text="Only search dereplicated Starships", size="sm"):
    """Create a switch component for toggling dereplicated-only searches."""
    return dmc.Switch(id="dereplicated-input", label=text, size=size, checked=True)


def create_accession_modal(accession):
    try:
        initial_df = fetch_meta_data()

        accession = str(accession).strip("[]").split("/")[-1].strip()
        initial_df["accession_tag"] = (
            initial_df["accession_tag"]
            .astype(str)
            .apply(lambda x: x.strip("[]").split("/")[-1].strip())
        )

        modal_data = initial_df[initial_df["accession_tag"] == accession]

        if modal_data.empty:
            return dmc.Stack(
                [
                    dmc.Alert(
                        title="No Data Found",
                        color="red",
                        children=[
                            f"No data found for accession: {accession}",
                            dmc.Space(h=10),
                            "Cache Status:",
                            dmc.List(
                                [
                                    dmc.ListItem(
                                        f"Total records in cache: {len(initial_df)}"
                                    ),
                                    dmc.ListItem(
                                        f"Sample accessions: {', '.join(initial_df['accession_tag'].head().tolist())}"
                                    ),
                                    dmc.ListItem(f"Searched for: {accession}"),
                                ]
                            ),
                        ],
                    )
                ]
            ), f"Accession: {accession}"

        # Basic ship information section
        ship_info = dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Ship Details", order=4, mb=10),
                dmc.SimpleGrid(
                    cols={"base": 1, "sm": 2},
                    spacing="lg",
                    children=[
                        dmc.Group(
                            [
                                dmc.Text("starshipID:", fw=700),
                                dmc.Text(modal_data["starshipID"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Curation Status:", fw=700),
                                dmc.Badge(
                                    modal_data["curated_status"].iloc[0],
                                    color="green"
                                    if modal_data["curated_status"].iloc[0] == "curated"
                                    else "yellow",
                                ),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Starship Family:", fw=700),
                                dmc.Text(modal_data["familyName"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Genomes Present:", fw=700),
                                dmc.Badge(str(len(modal_data)), color="blue"),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Starship Navis:", fw=700),
                                dmc.Text(modal_data["starship_navis"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Starship Haplotype:", fw=700),
                                dmc.Text(modal_data["starship_haplotype"].iloc[0]),
                            ]
                        ),
                    ],
                ),
            ],
        )

        # Taxonomy section
        taxonomy_info = dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Taxonomy", order=4, mb=10),
                dmc.SimpleGrid(
                    cols={"base": 1, "sm": 2},
                    spacing="lg",
                    children=[
                        dmc.Group(
                            [
                                dmc.Text("Order:", fw=700),
                                dmc.Text(modal_data["order"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Family:", fw=700),
                                dmc.Text(modal_data["family"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Species:", fw=700),
                                dmc.Text(modal_data["name"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Strain:", fw=700),
                                dmc.Text(modal_data["strain"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("NCBI Taxonomy ID:", fw=700),
                                dmc.Anchor(
                                    str(
                                        int(modal_data["taxID"].iloc[0])
                                        if isinstance(
                                            modal_data["taxID"].iloc[0], (float, int)
                                        )
                                        else modal_data["taxID"].iloc[0]
                                    ),
                                    href=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={int(modal_data['taxID'].iloc[0]) if isinstance(modal_data['taxID'].iloc[0], (float, int)) else modal_data['taxID'].iloc[0]}",
                                    target="_blank",
                                ),
                            ]
                        ),
                    ],
                ),
            ],
        )

        # Genome details section
        genome_details = dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Genome Details", order=4, mb=10),
                dmc.ScrollArea(
                    children=[
                        dmc.Table(
                            striped=True,
                            highlightOnHover=True,
                            withColumnBorders=True,
                            children=[
                                html.Thead(
                                    html.Tr(
                                        [
                                            html.Th(
                                                "Field",
                                                style={"backgroundColor": "#f8f9fa"},
                                            ),
                                            *[
                                                html.Th(
                                                    f"{modal_data['assembly_accession'].iloc[i] if 'assembly_accession' in modal_data else f'Genome {i + 1}'}",
                                                    style={
                                                        "backgroundColor": "#f8f9fa"
                                                    },
                                                )
                                                for i in range(len(modal_data))
                                            ],
                                        ]
                                    )
                                ),
                                html.Tbody(
                                    [
                                        html.Tr(
                                            [
                                                html.Td("Genome Source"),
                                                *[
                                                    html.Td(
                                                        modal_data["genomeSource"].iloc[
                                                            i
                                                        ]
                                                    )
                                                    for i in range(len(modal_data))
                                                ],
                                            ]
                                        ),
                                        html.Tr(
                                            [
                                                html.Td("ContigID"),
                                                *[
                                                    html.Td(
                                                        modal_data["contigID"].iloc[i]
                                                    )
                                                    for i in range(len(modal_data))
                                                ],
                                            ]
                                        ),
                                        html.Tr(
                                            [
                                                html.Td("Element Position"),
                                                *[
                                                    html.Td(
                                                        f"{int(modal_data['elementBegin'].iloc[i])} - {int(modal_data['elementEnd'].iloc[i])}"
                                                        if pd.notna(
                                                            modal_data[
                                                                "elementBegin"
                                                            ].iloc[i]
                                                        )
                                                        and pd.notna(
                                                            modal_data[
                                                                "elementEnd"
                                                            ].iloc[i]
                                                        )
                                                        else ""
                                                    )
                                                    for i in range(len(modal_data))
                                                ],
                                            ]
                                        ),
                                        html.Tr(
                                            [
                                                html.Td("Size"),
                                                *[
                                                    html.Td(
                                                        f"{int(modal_data['size'].iloc[i])} bp"
                                                    )
                                                    for i in range(len(modal_data))
                                                ],
                                            ]
                                        ),
                                    ]
                                ),
                            ],
                        )
                    ]
                )
                if len(modal_data) > 1
                else dmc.SimpleGrid(
                    cols={"base": 1, "sm": 2},
                    spacing="lg",
                    children=[
                        dmc.Group(
                            [
                                dmc.Text("Assembly Accession:", fw=700),
                                dmc.Text(
                                    modal_data["assembly_accession"].iloc[0]
                                    if "assembly_accession" in modal_data
                                    else ""
                                ),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Genome Source:", fw=700),
                                dmc.Text(modal_data["genomeSource"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("ContigID:", fw=700),
                                dmc.Text(modal_data["contigID"].iloc[0]),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Element Position:", fw=700),
                                dmc.Text(
                                    f"{int(modal_data['elementBegin'].iloc[0])} - {int(modal_data['elementEnd'].iloc[0])}"
                                    if pd.notna(modal_data["elementBegin"].iloc[0])
                                    and pd.notna(modal_data["elementEnd"].iloc[0])
                                    else ""
                                ),
                            ]
                        ),
                        dmc.Group(
                            [
                                dmc.Text("Size:", fw=700),
                                dmc.Text(f"{int(modal_data['size'].iloc[0])} bp"),
                            ]
                        ),
                    ],
                ),
            ],
        )

        modal_content = dmc.Stack([ship_info, taxonomy_info, genome_details], gap="md")

        modal_title = dmc.Title(f"Ship Accession: {accession}", order=2)
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
        [
            Input(table_id, "cellClicked"),  # AG Grid
            Input(table_id, "active_cell"),
        ],  # Dash DataTable
        [
            State(table_id, "derived_virtual_data"),
            State(table_id, "page_current"),
            State(table_id, "page_size"),
        ],
        prevent_initial_call=True,
    )
    def toggle_modal(cell_clicked, active_cell, table_data, page_current, page_size):
        try:
            ctx = callback_context
            triggered_id = ctx.triggered[0]["prop_id"]

            if not (cell_clicked or active_cell):
                return False, no_update, no_update

            accession = None
            if f"{table_id}.cellClicked" in triggered_id and cell_clicked:
                # AG Grid format
                if cell_clicked["colId"] == "accession_tag":
                    accession = str(cell_clicked["value"])
            elif f"{table_id}.active_cell" in triggered_id and active_cell:
                # Dash DataTable format
                if active_cell["column_id"] == "accession_tag":
                    # Calculate the actual row index based on pagination
                    actual_row_idx = (page_current or 0) * page_size + active_cell[
                        "row"
                    ]
                    if table_data and actual_row_idx < len(table_data):
                        accession = str(table_data[actual_row_idx]["accession_tag"])

            if accession:
                # Clean and standardize the accession tag
                accession = accession.strip("[]").split("/")[-1].strip()
                logger.debug(f"Looking for accession in cache: {accession}")

                modal_content, modal_title = create_accession_modal(accession)
                return True, modal_content, modal_title

            return False, no_update, no_update

        except Exception as e:
            logger.error(f"Error in toggle_modal: {str(e)}")
            logger.error(traceback.format_exc())
            error_content = html.Div(
                [
                    html.P("Error loading modal content"),
                    html.P(f"Details: {str(e)}"),
                ]
            )
            return True, error_content, "Error"


def create_file_upload(
    upload_id: str,
    output_id: str,
    accept_types: List[str],
    placeholder_text: str = "Drag and drop or click to select a file",
    icon: str = "mdi:file-upload",
    **kwargs,
) -> dmc.Stack:
    return dmc.Stack(
        [
            dmc.Center(DashIconify(icon=icon, width=40, height=40, color="#228be6")),
            dcc.Upload(
                id=upload_id,
                children=dmc.Stack(
                    [
                        html.Div(
                            id=output_id,
                            children=placeholder_text,
                        ),
                        dmc.Text(
                            f"Accepted formats: {', '.join(accept_types)}",
                            size="sm",
                            c="dimmed",
                        ),
                    ],
                    align="center",
                    gap="xs",
                ),
                className="upload-box",
                **kwargs,
            ),
            dmc.Progress(
                id=f"{upload_id}-progress",
                value=0,
                animated=True,
                style={"display": "none"},
            ),
        ],
        gap="md",
    )


def create_feedback_button():
    return dmc.Notification(
        title="Found an issue?",
        id="feedback-notify",
        action="show",
        # icon=DashIconify(icon="octicon:mark-github-16", width=20),
        autoClose=20000,
        color="gray",
        radius="md",
        message=[
            dmc.Anchor(
                dmc.Button(
                    "Report it on GitHub",
                    variant="light",
                    color="blue",
                    size="sm",
                    leftSection=DashIconify(icon="octicon:mark-github-16", width=20),
                ),
                href="https://github.com/FungAGE/starbase/issues",
                target="_blank",
            )
        ],
        className="d-none d-md-block",
        style={
            "width": "250px",
            "backgroundColor": "var(--mantine-color-gray-1)",
        },
    )


def make_progress_bar(message, id_prefix, color):
    """Create a progress bar with given prefix and color"""
    return dmc.Group(
        [
            dmc.Text(message, size="sm"),
            dbc.Progress(
                id=f"{id_prefix}-progress",
                value=0,
                color=color,
                animated=False,
                striped=False,
                style={"width": "100%", "marginBottom": "5px"},
            ),
            html.Div(id=f"{id_prefix}-progress-spacer"),  # Spacer div
        ]
    )

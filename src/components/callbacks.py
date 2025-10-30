from dash import html, Output, Input, State, callback, no_update, dcc, callback_context
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash_iconify import DashIconify
from typing import List

import functools
import traceback
import pandas as pd

from src.database.sql_manager import fetch_meta_data, get_quality_tags

from src.config.logging import get_logger

logger = get_logger(__name__)

def handle_callback_error(callback_func):
    """
    Decorator to handle callback errors gracefully

    Args:
        callback_func: The callback function to wrap

    Returns:
        Wrapped function with error handling
    """

    @functools.wraps(callback_func)
    def wrapper(*args, **kwargs):
        try:
            return callback_func(*args, **kwargs)
        except PreventUpdate:
            raise
        except Exception as e:
            # Log the error with full traceback
            logger.error(f"Callback error in {callback_func.__name__}: {str(e)}")
            logger.error(f"Inputs: args={args}, kwargs={kwargs}")
            logger.error(traceback.format_exc())

            # Return detailed error alert
            return dmc.Alert(
                title="Error Loading Data",
                children=[
                    "We encountered a problem processing your request.",
                    dmc.Space(h=10),
                    dmc.Text(f"Error: {str(e)}", size="sm"),
                    dmc.Space(h=10),
                    dmc.Code(str(traceback.format_exc()), block=True),
                ],
                color="red",
                variant="filled",
            )

    return wrapper


def safe_get_value(df, column, index=0, default="N/A", format_func=None):
    """
    Safely extract a value from a DataFrame, handling NA/Null values.

    Args:
        df: DataFrame to extract from
        column: Column name
        index: Row index (default 0)
        default: Default value if NA/Null (default "N/A")
        format_func: Optional function to format the value

    Returns:
        Formatted value or default
    """
    try:
        if column not in df.columns:
            return default

        value = df[column].iloc[index]

        # Check for various forms of null/empty values
        if (
            pd.isna(value)
            or value is None
            or value == ""
            or str(value).lower() in ["nan", "none", "null"]
        ):
            return default

        # Apply formatting function if provided
        if format_func and callable(format_func):
            try:
                return format_func(value)
            except (ValueError, TypeError):
                return default

        return str(value)

    except (IndexError, KeyError):
        return default


def safe_get_numeric(df, column, index=0, default="N/A"):
    """Helper for numeric values"""
    return safe_get_value(df, column, index, default, lambda x: str(int(float(x))))


def safe_get_position(df, begin_col, end_col, index=0, default="N/A"):
    """Helper for position ranges"""
    begin = safe_get_value(df, begin_col, index, None, lambda x: int(float(x)))
    end = safe_get_value(df, end_col, index, None, lambda x: int(float(x)))

    if begin != "N/A" and end != "N/A" and begin is not None and end is not None:
        return f"{begin} - {end}"
    return default


def create_table_cells(df, column, format_func=None, default="N/A"):
    """
    Create table cells for all rows in a DataFrame column, handling NA values.

    Args:
        df: DataFrame
        column: Column name
        format_func: Optional formatting function
        default: Default value for NA/Null

    Returns:
        List of html.Td elements
    """
    cells = []
    for i in range(len(df)):
        value = safe_get_value(df, column, i, default, format_func)
        cells.append(html.Td(value))
    return cells


def create_position_cells(df, begin_col, end_col, default="N/A"):
    """Create table cells for position ranges"""
    cells = []
    for i in range(len(df)):
        position = safe_get_position(df, begin_col, end_col, i, default)
        cells.append(html.Td(position))
    return cells


def create_genome_headers(df):
    """
    Create table headers for genome columns, using assembly_accession if available,
    otherwise fallback to genome numbers.
    """
    headers = []
    for i in range(len(df)):
        # Try to get assembly accession as header
        header_value = safe_get_value(df, "assembly_accession", i, default="")

        # If no assembly accession, use genome number
        if not header_value or header_value == "N/A":
            header_value = f"Genome {i + 1}"

        headers.append(html.Th(header_value, style={"backgroundColor": "#f8f9fa"}))
    return headers


def create_genome_cards(df):
    """
    Create a list of card components, one per genome, to avoid wide tables
    that overflow small viewports.
    """
    cards = []
    for i in range(len(df)):
        header_value = safe_get_value(df, "assembly_accession", i, default="")
        if not header_value or header_value == "N/A":
            header_value = f"Genome {i + 1}"

        genome_source = safe_get_value(df, "genomeSource", i, default="N/A")
        contig_id = safe_get_value(df, "contigID", i, default="N/A")
        position = safe_get_position(df, "elementBegin", "elementEnd", i, default="N/A")
        length_bp = safe_get_value(
            df,
            "elementLength",
            i,
            default="N/A",
            format_func=lambda x: f"{int(float(x))} bp",
        )

        cards.append(
            dmc.Paper(
                p="md",
                withBorder=True,
                radius="sm",
                children=[
                    dmc.Text(header_value, fw=700, mb=6),
                    dmc.Stack(
                        gap="xs",
                        children=[
                            dmc.Group([dmc.Text("Genome Source:", fw=700), dmc.Text(genome_source)]),
                            dmc.Group([dmc.Text("ContigID:", fw=700), dmc.Text(contig_id)]),
                            dmc.Group([dmc.Text("Element Position:", fw=700), dmc.Text(position)]),
                            dmc.Group([dmc.Text("Size:", fw=700), dmc.Text(length_bp)]),
                        ],
                    ),
                ],
            )
        )

    return cards


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


def create_quality_tag_badges(quality_tags):
    """
    Create badge components for quality tags.
    
    Args:
        quality_tags (list): List of tag strings in format "tag_type" or "tag_type:tag_value"
    
    Returns:
        list: List of dmc.Badge components
    """
    if not quality_tags:
        return []
    
    # Define colors for different tag types
    tag_colors = {
        "incomplete": "orange",
        "fragmented": "red",
        "partial": "yellow",
        "nested": "grape",
        "verified": "teal",
        "high_quality": "green",
        "low_quality": "red",
        "default": "gray"
    }
    
    badges = []
    for tag in quality_tags:
        # Parse tag_type and tag_value if present
        if ":" in tag:
            tag_type, tag_value = tag.split(":", 1)
            display_text = f"{tag_type}: {tag_value}"
        else:
            tag_type = tag
            display_text = tag_type
        
        # Get color for this tag type
        color = tag_colors.get(tag_type.lower(), tag_colors["default"])
        
        badges.append(
            dmc.Badge(
                display_text,
                color=color,
                variant="light",
                size="xs",
            )
        )
    
    return badges


def dereplicated_switch(text="Only search dereplicated Starships", size="sm"):
    """Create a switch component for toggling dereplicated-only searches."""
    return dmc.Switch(id="dereplicated-input", label=text, size=size, checked=True)

def create_accession_modal_data(accession):
    """Create structured data for accession modal instead of Dash components."""
    try:
        accession = str(accession).strip("[]").split("/")[-1].strip()
        base_accession = accession.split(".")[0] if "." in accession else accession

        modal_data = fetch_meta_data(accession_tags=[base_accession])

        if modal_data.empty:
            return {
                "title": f"Accession: {accession}",
                "error": f"No data found for accession: {accession}"
            }

        # Validate modal_data
        if not isinstance(modal_data, pd.DataFrame) or modal_data.empty:
            return {
                "title": f"Accession: {accession}",
                "error": "Invalid modal data received"
            }

        # Check for required columns
        required_columns = ["accession_tag", "familyName"]
        missing_columns = [
            col for col in required_columns if col not in modal_data.columns
        ]
        if missing_columns:
            return {
                "title": f"Accession: {accession}",
                "error": f"Missing required columns: {missing_columns}"
            }

        # HACK: applying a fix for extra rows in the starship_features table, only take the first begin/end coordinates for each ship_id/accession_id
        # ! this might cause some issues if coordinates are not updated for all rows for a ship_id/accession_id pair, updated only if begin/end coordinates are the same
        # TODO: split features table or move coordinate information to separate table or another existing table
        modal_data = modal_data.groupby("accession_tag").first().reset_index()

        # Fetch quality tags separately using joined_ship_id
        joined_ship_id = safe_get_numeric(modal_data, "joined_ship_id")
        quality_tags = []
        accepted_quality_tags = ["missing_direct_repeats","missing_tir","missing_boundaries","missing_genome_context","unannotated","missing_empty_site"]
        if joined_ship_id:
            try:
                quality_tags_data = get_quality_tags(joined_ship_id)
                # Format tags as "tag_type:tag_value" or just "tag_type" if no value
                for tag in quality_tags_data:
                    if tag.get("tag_value"):
                        tag_value = f"{tag['tag_type']}:{tag['tag_value']}"
                    else:
                        tag_value = tag['tag_type']
                    if tag_value in accepted_quality_tags:
                        quality_tags.append(tag_value)
            except Exception as e:
                logger.warning(f"Error fetching quality tags for joined_ship_id {joined_ship_id}: {e}")

        # Create structured data
        result = {
            "title": f"Starship Accession: {accession}",
            "familyName": safe_get_value(modal_data, "familyName"),
            "taxonomic_family": safe_get_value(modal_data, "family"),
            "genomes_present": str(len(modal_data)),
            "navis_name": safe_get_value(modal_data, "navis_name"),
            "haplotype_name": safe_get_value(modal_data, "haplotype_name"),
            "order": safe_get_value(modal_data, "order"),
            "species_name": safe_get_value(modal_data, "name"),
            "tax_id": safe_get_numeric(modal_data, "taxID"),
            "assembly_accession": safe_get_value(modal_data, "assembly_accession", default=""),
            "genome_source": safe_get_value(modal_data, "genomeSource", default=""),
            "contig_id": safe_get_value(modal_data, "contigID", default=""),
            "element_length": safe_get_numeric(modal_data, "elementLength", default=""),
            "element_position": safe_get_position(modal_data, "elementBegin", "elementEnd"),
            "curated_status": safe_get_value(modal_data, "curated_status", default="unknown"),
            "quality_tags": quality_tags
        }

        return result

    except Exception as e:
        logger.error(f"Error in create_accession_modal_data: {str(e)}")
        logger.error(traceback.format_exc())
        raise

def create_accession_modal(modal_data):
    accession = modal_data["accession_tag"]
    familyName = modal_data["familyName"]
    genomes_present = modal_data["genomes_present"]
    navis_name = modal_data["navis_name"]
    haplotype_name = modal_data["haplotype_name"]
    order = modal_data["order"]
    taxonomic_family = modal_data["taxonomic_family"]
    species_name = modal_data["species_name"]
    tax_id = modal_data["tax_id"]
    assembly_accession = modal_data["assembly_accession"]
    genome_source = modal_data["genome_source"]
    contig_id = modal_data["contig_id"]
    element_length = modal_data["element_length"]
    element_position = modal_data["element_position"]
    curated_status = modal_data["curated_status"]
    quality_tags = modal_data["quality_tags"]
    badge_color = modal_data["badge_color"]
    try:
        modal_title = dmc.Title(f"Starship Accession: {accession}", order=2)

        # Basic ship information section
        ship_info = dmc.Paper(
                        p="md",
                        withBorder=True,
                        children=[
                            dmc.SimpleGrid(
                                cols={"base": 1, "sm": 2},
                                spacing="lg",
                                children=[
                                        dmc.Group(
                                            [
                                                dmc.Text("Starship Family:", fw=700),
                                                dmc.Text(familyName),
                                            ]
                                        ),
                                        dmc.Group(
                                            [
                                                dmc.Text("Genomes Present:", fw=700),
                                                dmc.Badge(str(genomes_present), color="blue"),
                                            ]
                                        ),
                                        dmc.Group(
                                            [
                                                dmc.Text("Starship Navis:", fw=700),
                                                dmc.Text(navis_name),
                                            ]
                                        ),
                                        dmc.Group(
                                            [
                                                dmc.Text("Starship Haplotype:", fw=700),
                                                dmc.Text(haplotype_name),
                                            ]
                                        ),
                                    ]
                                    )
                        ])
        modal_data.append(ship_info)

        # Add quality tags if present
        if quality_tags:
            quality_tag_badges = create_quality_tag_badges(quality_tags)

            curation_info = dmc.Paper(
                p="md",
                withBorder=True,
                children=[
                            dmc.Stack([
                                dmc.Group(
                                    [
                                        dmc.Text("Curation Status:", fw=700),
                                        dmc.Badge(
                                            curated_status,
                                            color=badge_color,
                                        ),
                                    ]
                                ),
                                dmc.Stack([
                                    dmc.Text("Quality Tags:", fw=700),
                                    dmc.Flex(quality_tag_badges, wrap="wrap", gap="xs"),
                                ])
                            ])
                        ]
            )
            modal_data.append(curation_info)
        
        # Taxonomy section
        taxonomy_info = dmc.Paper(
                    p="md",
                    withBorder=True,
                    children=[
                        dmc.SimpleGrid(
                            cols={"base": 1, "sm": 2},
                            spacing="lg",
                            children=[
                                dmc.Group(
                                    [
                                        dmc.Text("Order:", fw=700),
                                        dmc.Text(order),
                                    ]
                                ),
                                dmc.Group(
                                    [
                                        dmc.Text("Family:", fw=700),
                                        dmc.Text(taxonomic_family),
                                    ]
                                ),
                                dmc.Group(
                                    [
                                        dmc.Text("Species:", fw=700),
                                        dmc.Text(species_name),
                                    ]
                                ),
                                dmc.Group(
                                    [
                                        dmc.Text("NCBI Taxonomy ID:", fw=700),
                                        (
                                            dmc.Anchor(
                                                tax_id,
                                                href=(
                                                    f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}"
                                                    if tax_id != "N/A"
                                                    else "#"
                                                ),
                                                target="_blank" if tax_id != "N/A" else None,
                                            )
                                            if tax_id != "N/A"
                                            else dmc.Text(tax_id)
                                        ),
                                    ]
                                ),
                            ],
                        ),
                    ],
                )
        modal_data.append(taxonomy_info)

        # Genome details section
        genome_details = dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Genome Details", order=4, mb=10),
                (
                    dmc.SimpleGrid(
                        cols={"base": 1, "sm": 2, "md": 3},
                        spacing="lg",
                        children=create_genome_cards(modal_data),
                    )
                    if len(modal_data) > 1
                    else dmc.SimpleGrid(
                        cols={"base": 1, "sm": 2},
                        spacing="lg",
                        children=[
                            dmc.Group(
                                [
                                    dmc.Text("Assembly Accession:", fw=700),
                                    dmc.Text(assembly_accession or "N/A"),
                                ]
                            ),
                            dmc.Group(
                                [
                                    dmc.Text("Genome Source:", fw=700),
                                    dmc.Text(genome_source or "N/A"),
                                ]
                            ),
                            dmc.Group(
                                [
                                    dmc.Text("ContigID:", fw=700),
                                    dmc.Text(contig_id or "N/A"),
                                ]
                            ),
                            dmc.Group(
                                [
                                    dmc.Text("Element Position:", fw=700),
                                    dmc.Text(element_position),
                                ]
                            ),
                            dmc.Group(
                                [
                                    dmc.Text("Size:", fw=700),
                                    dmc.Text(
                                        f"{element_length} bp"
                                        if element_length != "N/A"
                                        else "N/A"
                                    ),
                                ]
                            ),
                        ],
                    )
                ),
            ],
        )
        modal_data.append(genome_details)

        return dmc.Stack(modal_data, gap="md"), modal_title

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
                # AG Grid format - handle both accession_tag and accession_display
                if cell_clicked["colId"] in ["accession_tag", "accession_display"]:
                    accession = str(cell_clicked["value"])
            elif f"{table_id}.active_cell" in triggered_id and active_cell:
                # Dash DataTable format - handle both accession_tag and accession_display
                if active_cell["column_id"] in ["accession_tag", "accession_display"]:
                    # Calculate the actual row index based on pagination
                    actual_row_idx = (page_current or 0) * page_size + active_cell[
                        "row"
                    ]
                    if table_data and actual_row_idx < len(table_data):
                        accession = str(
                            table_data[actual_row_idx][active_cell["column_id"]]
                        )

            if accession:
                # Clean and standardize the accession tag
                accession = accession.strip("[]").split("/")[-1].strip()
                logger.debug(f"Looking for accession in cache: {accession}")

                # Instead of opening a Dash modal, trigger the universal modal via JavaScript
                return (
                    False,  # Don't open the Dash modal
                    no_update,
                    no_update
                )

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
            error_title = dmc.Title("Error", order=3)
            return True, error_content, error_title

    return toggle_modal


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

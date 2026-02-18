from dash import html, dcc
from typing import List
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import traceback

from src.config.logging import get_logger
from src.components.data import safe_get_value, safe_get_position
from src.database.sql_manager import get_database_version

logger = get_logger(__name__)


def curated_switch(text="Only search curated Starships", size="sm"):
    """Create a switch component for toggling curated-only searches."""
    return dmc.Switch(id="curated-input", label=text, size=size, checked=True)


def dereplicated_switch(text="Only search dereplicated Starships", size="sm"):
    """Create a switch component for toggling dereplicated-only searches."""
    return dmc.Switch(id="dereplicated-input", label=text, size=size, checked=True)


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
        "default": "gray",
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
                            dmc.Group(
                                [
                                    dmc.Text("Genome Source:", fw=700),
                                    dmc.Text(genome_source),
                                ]
                            ),
                            dmc.Group(
                                [dmc.Text("ContigID:", fw=700), dmc.Text(contig_id)]
                            ),
                            dmc.Group(
                                [
                                    dmc.Text("Element Position:", fw=700),
                                    dmc.Text(position),
                                ]
                            ),
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


def create_database_version_indicator():
    """Create a database version indicator for the bottom-left corner"""
    try:
        db_version = get_database_version()
        version_text = f"v{db_version}" if db_version != "unknown" else "Unknown"
    except Exception as e:
        logger.error(f"Error fetching database version: {str(e)}")
        version_text = "Error"

    return dmc.Notification(
        title="Database Version",
        id="db-version-notify",
        action="show",
        autoClose=30000,  # Auto-close after 30 seconds (longer than feedback button)
        color="blue",
        radius="md",
        message=[
            dmc.Stack(
                [
                    dmc.Group(
                        [
                            DashIconify(icon="mdi:database", width=16, color="white"),
                            dmc.Text(
                                version_text,
                                size="sm",
                                c="white",
                                fw="bold",
                            ),
                        ],
                        gap="xs",
                        align="center",
                    ),
                ],
                gap="xs",
                align="flex-start",
            )
        ],
        className="d-none d-md-block",
        style={
            "width": "200px",
            "backgroundColor": "var(--mantine-color-blue-6)",
            "color": "white",
        },
    )


def create_modal(modal_data: dict, mode: str = "ship_accession") -> tuple:
    """
    Create a flexible modal component from structured data.

    Args:
        modal_data: Dictionary containing modal data fields
        mode: Modal type ("ship_accession", "accession", or custom)

    Returns:
        Tuple of (modal_content: dmc.Stack, modal_title: dmc.Title)
    """
    try:
        # Handle error cases
        if modal_data.get("error"):
            return (
                dmc.Alert(
                    modal_data["error"],
                    title="Error",
                    color="red",
                    variant="light",
                ),
                dmc.Title(modal_data.get("title", "Error"), order=2),
            )

        # Generate title based on mode and data
        modal_title = _generate_modal_title(modal_data, mode)

        # Build sections based on available data
        sections = []

        # Ship information section (if applicable)
        if _has_ship_data(modal_data):
            sections.append(_create_ship_info_section(modal_data))

        # Quality tags and curation section
        if modal_data.get("quality_tags") or modal_data.get("curated_status"):
            sections.append(_create_curation_section(modal_data))

        # Taxonomy section
        if _has_taxonomy_data(modal_data):
            sections.append(_create_taxonomy_section(modal_data))

        # Genome details section
        if _has_genome_data(modal_data):
            sections.append(_create_genome_section(modal_data))

        # Additional sections can be added here for future modal types

        if not sections:
            sections.append(
                dmc.Alert(
                    "No additional information available.",
                    title="Information",
                    color="blue",
                    variant="light",
                )
            )

        return dmc.Stack(sections, gap="md"), modal_title

    except Exception as e:
        logger.error(f"Error in create_modal: {str(e)}")
        logger.error(traceback.format_exc())
        raise


def _generate_modal_title(modal_data: dict, mode: str) -> dmc.Title:
    """Generate modal title based on mode and data."""
    accession_tag = modal_data.get("accession_tag") or modal_data.get(
        "title", "Unknown"
    )

    title_map = {
        "ship_accession": f"Starship Accession: {accession_tag}",
        "accession": f"Accession: {accession_tag}",
    }

    title_text = title_map.get(mode, f"Accession: {accession_tag}")
    return dmc.Title(title_text, order=2)


def _has_ship_data(modal_data: dict) -> bool:
    """Check if modal data contains ship-specific information."""
    return any(
        modal_data.get(key)
        for key in ["familyName", "navis_name", "haplotype_name", "genomes_present"]
    )


def _create_ship_info_section(modal_data: dict) -> dmc.Paper:
    """Create ship information section."""
    fields = []

    # Family name
    if modal_data.get("familyName"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Starship Family:", fw=700),
                    dmc.Text(modal_data["familyName"]),
                ]
            )
        )

    # Genomes present (calculate if not provided)
    genomes_present = modal_data.get("genomes_present")
    if genomes_present is not None:
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Genomes Present:", fw=700),
                    dmc.Badge(str(genomes_present), color="blue"),
                ]
            )
        )

    # Navis name
    if modal_data.get("navis_name"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Starship Navis:", fw=700),
                    dmc.Text(modal_data["navis_name"]),
                ]
            )
        )

    # Haplotype name
    if modal_data.get("haplotype_name"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Starship Haplotype:", fw=700),
                    dmc.Text(modal_data["haplotype_name"]),
                ]
            )
        )

    if not fields:
        return None

    return dmc.Paper(
        p="md",
        withBorder=True,
        children=[
            dmc.SimpleGrid(
                cols={"base": 1, "sm": 2},
                spacing="lg",
                children=fields,
            )
        ],
    )


def _create_curation_section(modal_data: dict) -> dmc.Paper:
    """Create curation status and quality tags section."""
    curation_items = []

    # Curation status
    curated_status = modal_data.get("curated_status", "unknown")
    if curated_status:
        badge_color = _get_curation_badge_color(curated_status)
        curation_items.append(
            dmc.Group(
                [
                    dmc.Text("Curation Status:", fw=700),
                    dmc.Badge(curated_status, color=badge_color),
                ]
            )
        )

    # Quality tags
    quality_tags = modal_data.get("quality_tags", [])
    if quality_tags:
        quality_tag_badges = create_quality_tag_badges(quality_tags)
        curation_items.append(
            dmc.Stack(
                [
                    dmc.Text("Quality Tags:", fw=700),
                    dmc.Flex(quality_tag_badges, wrap="wrap", gap="xs"),
                ]
            )
        )

    if not curation_items:
        return None

    return dmc.Paper(
        p="md",
        withBorder=True,
        children=[dmc.Stack(curation_items)],
    )


def _get_curation_badge_color(curated_status: str) -> str:
    """Get badge color based on curation status."""
    color_map = {
        "curated": "green",
        "uncurated": "gray",
        "needs_review": "yellow",
        "pending": "blue",
        "unknown": "gray",
    }
    return color_map.get(curated_status.lower(), "gray")


def _has_taxonomy_data(modal_data: dict) -> bool:
    """Check if modal data contains taxonomy information."""
    return any(
        modal_data.get(key)
        for key in ["order", "taxonomic_family", "species_name", "tax_id"]
    )


def _create_taxonomy_section(modal_data: dict) -> dmc.Paper:
    """Create taxonomy information section."""
    fields = []

    # Order
    if modal_data.get("order"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Order:", fw=700),
                    dmc.Text(modal_data["order"]),
                ]
            )
        )

    # Family
    if modal_data.get("taxonomic_family"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Family:", fw=700),
                    dmc.Text(modal_data["taxonomic_family"]),
                ]
            )
        )

    # Species
    if modal_data.get("species_name"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Species:", fw=700),
                    dmc.Text(modal_data["species_name"]),
                ]
            )
        )

    # NCBI Taxonomy ID
    tax_id = modal_data.get("tax_id")
    if tax_id and tax_id != "N/A":
        tax_id_component = dmc.Anchor(
            str(tax_id),
            href=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={tax_id}",
            target="_blank",
        )
        fields.append(
            dmc.Group(
                [
                    dmc.Text("NCBI Taxonomy ID:", fw=700),
                    tax_id_component,
                ]
            )
        )
    elif tax_id:
        fields.append(
            dmc.Group(
                [
                    dmc.Text("NCBI Taxonomy ID:", fw=700),
                    dmc.Text(str(tax_id)),
                ]
            )
        )

    if not fields:
        return None

    return dmc.Paper(
        p="md",
        withBorder=True,
        children=[
            dmc.SimpleGrid(
                cols={"base": 1, "sm": 2},
                spacing="lg",
                children=fields,
            )
        ],
    )


def _has_genome_data(modal_data: dict) -> bool:
    """Check if modal data contains genome information."""
    return any(
        modal_data.get(key)
        for key in [
            "assembly_accession",
            "genome_source",
            "contig_id",
            "element_position",
            "element_length",
        ]
    )


def _create_genome_section(modal_data: dict) -> dmc.Paper:
    """Create genome details section."""
    # Check if we have multiple genomes (for genome cards)
    # This would need to be determined from the data structure
    # For now, assume single genome if not specified

    fields = []

    # Assembly Accession
    if modal_data.get("assembly_accession"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Assembly Accession:", fw=700),
                    dmc.Text(modal_data["assembly_accession"]),
                ]
            )
        )

    # Genome Source
    if modal_data.get("genome_source"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Genome Source:", fw=700),
                    dmc.Text(modal_data["genome_source"]),
                ]
            )
        )

    # Contig ID
    if modal_data.get("contig_id"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("ContigID:", fw=700),
                    dmc.Text(modal_data["contig_id"]),
                ]
            )
        )

    # Element Position
    if modal_data.get("element_position"):
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Element Position:", fw=700),
                    dmc.Text(modal_data["element_position"]),
                ]
            )
        )

    # Element Length
    element_length = modal_data.get("element_length")
    if element_length and element_length != "N/A":
        fields.append(
            dmc.Group(
                [
                    dmc.Text("Size:", fw=700),
                    dmc.Text(f"{element_length} bp"),
                ]
            )
        )

    if not fields:
        return None

    # Check if we should use genome cards (if multiple genomes present)
    # This would require checking if modal_data contains a list of genomes
    # For now, use simple grid
    use_genome_cards = False  # TODO: Determine from data structure

    if use_genome_cards and hasattr(modal_data, "genomes"):
        # Use genome cards for multiple genomes
        genome_cards = create_genome_cards(modal_data.get("genomes", []))
        return dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Genome Details", order=4, mb=10),
                dmc.SimpleGrid(
                    cols={"base": 1, "sm": 2, "md": 3},
                    spacing="lg",
                    children=genome_cards,
                ),
            ],
        )
    else:
        # Use simple grid for single genome
        return dmc.Paper(
            p="md",
            withBorder=True,
            children=[
                dmc.Title("Genome Details", order=4, mb=10),
                dmc.SimpleGrid(
                    cols={"base": 1, "sm": 2},
                    spacing="lg",
                    children=fields,
                ),
            ],
        )

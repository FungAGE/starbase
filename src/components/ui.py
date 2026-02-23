from dash import html, dcc
from typing import List
import dash_mantine_components as dmc
from dash_iconify import DashIconify

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
                            *(
                                [
                                    dmc.Group(
                                        [
                                            dmc.Text("Element Position:", fw=700),
                                            dmc.Text(position),
                                        ]
                                    )
                                ]
                                if position != "N/A"
                                else []
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

from dash import html, dcc
from typing import List
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from src.config.logging import get_logger
from src.database.sql_manager import get_database_version

logger = get_logger(__name__)


def curated_switch(text="Only search curated Starships", size="sm"):
    """Create a switch component for toggling curated-only searches."""
    return dmc.Switch(id="curated-input", label=text, size=size, checked=True, color="indigo")


def dereplicated_switch(text="Only search dereplicated Starships", size="sm"):
    """Create a switch component for toggling dereplicated-only searches."""
    return dmc.Switch(id="dereplicated-input", label=text, size=size, checked=True, color="indigo")


download_ships_button = dmc.Anchor(
    dmc.Button(
        [
            dmc.Group(
                [
                    DashIconify(icon="mdi:download", color="indigo"),
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
        variant="filled",
        color="indigo",
        size="lg",
        radius="md",
        fullWidth=False,
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
        "width": "auto",
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
                " data on our ",
                html.A(
                    "GitHub repo",
                    href="https://github.com/FungAGE/starbase",
                    target="_blank",
                ),
                ". This is a work in progress, and we are currently in the process of migrating to a new back-end, which will provide more options for data export.",
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
    return dcc.Upload(
        id=upload_id,
        children=dmc.Stack(
            [
                dmc.Center(
                    DashIconify(icon=icon, width=40, height=40, color="var(--mantine-color-indigo-7)")
                ),
                dmc.Stack(
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
                dmc.Progress(
                    id=f"{upload_id}-progress",
                    value=0,
                    animated=True,
                    style={"display": "none"},
                    c="var(--mantine-primary-color-6)",
                ),
            ],
            className="upload-box",
            style={"boxShadow": "none"},
            **kwargs,
            gap="md",
        ),
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
                    leftSection=DashIconify(icon="octicon:mark-github-16", width=20, color="indigo"),
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
    """Create a database version indicator (legacy notification - use create_footer instead)."""
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
        autoClose=30000,
        color="blue",
        radius="md",
        message=[
            dmc.Stack(
                [
                    dmc.Group(
                        [
                            DashIconify(icon="mdi:database", width=16, color="var(--mantine-color-white)"),
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


def create_footer():
    """Create app footer with database version and feedback link."""
    try:
        db_version = get_database_version()
        version_text = f"v{db_version}" if db_version != "unknown" else "Unknown"
    except Exception as e:
        logger.error(f"Error fetching database version: {str(e)}")
        version_text = "Error"

    return html.Footer(
        dmc.Group(
            [
                dmc.Group(
                    [
                        DashIconify(icon="mdi:database", width=16, color="var(--mantine-color-gray-6)"),
                        dmc.Text(version_text, size="sm", c="dimmed"),
                    ],
                    gap="xs",
                    align="center",
                ),
                dmc.Anchor(
                    "Report an issue",
                    href="https://github.com/FungAGE/starbase/issues",
                    target="_blank",
                    size="sm",
                    c="dimmed",
                    style={"textDecoration": "none"},
                ),
            ],
            justify="space-between",
            align="center",
            p="md",
            style={
                "borderTop": "1px solid var(--mantine-color-indigo-2)",
                "backgroundColor": "var(--mantine-color-indigo-0)",
            },
        ),
        style={
            "position": "fixed",
            "bottom": 0,
            "left": 0,
            "right": 0,
            "zIndex": 100,
        },
    )


source_code_card = dmc.Paper(
    children=[
        dmc.Stack(
            [
                dmc.Text(
                    [
                        "The source code for ",
                        html.Span("starbase", className="logo-text"),
                        " webserver will soon be available on GitHub",
                    ],
                    size="lg",
                ),
                dmc.Image(
                    src="assets/images/starbase-map.png",
                    fit="contain",
                    className="auto-resize-750",
                ),
            ],
            gap="xl",
        ),
    ],
    p="xl",
    radius="md",
    withBorder=True,
    h="100%",
    shadow="sm",
)

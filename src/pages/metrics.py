import dash
from dash import dcc, html, callback, Input, Output, no_update
import dash_mantine_components as dmc
from src.telemetry.utils import get_telemetry_visualizations
from src.config.logging import get_logger
import plotly.graph_objects as go
from datetime import datetime, timedelta
import traceback

dash.register_page(__name__)


logger = get_logger(__name__)


def get_metrics_layout():
    try:
        telemetry_visualizations = get_telemetry_visualizations()
        if telemetry_visualizations is None:
            telemetry_visualizations = {"unique_users": 0, "total_requests": 0}

        time_series_fig = telemetry_visualizations.get("time_series", go.Figure())
        endpoints_fig = telemetry_visualizations.get("endpoints", go.Figure())
        map_fig = telemetry_visualizations.get("map", go.Figure())

        layout = dmc.Container(
            [
                dcc.Interval(
                    id="telemetry-refresh-interval",
                    interval=60 * 1000,  # refresh every minute
                    n_intervals=0,
                ),
                dmc.Stack(
                    style={"paddingTop": "20px"},
                    children=[
                        # Header with controls
                        dmc.Group(
                            [
                                dmc.Title("Analytics Dashboard", order=1),
                                dmc.Badge(
                                    "Live Data",
                                    color="green",
                                    variant="light",
                                    size="lg",
                                ),
                            ],
                            justify="space-between",
                            align="center",
                        ),
                        dmc.Text(
                            "Real-time insights into visitor behavior and platform usage",
                            size="lg",
                            c="dimmed",
                        ),
                        # Key metrics cards
                        dmc.Grid(
                            style={"paddingTop": "20px"},
                            children=[
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Group(
                                                    [
                                                        dmc.Text("👥", size="xl"),
                                                        dmc.Stack(
                                                            [
                                                                dmc.Text(
                                                                    f"{telemetry_visualizations['unique_users']:,}",
                                                                    size="xl",
                                                                    fw="bold",
                                                                ),
                                                                dmc.Text(
                                                                    "Total Unique Users",
                                                                    size="sm",
                                                                    c="dimmed",
                                                                ),
                                                            ],
                                                            gap=0,
                                                        ),
                                                    ],
                                                    align="center",
                                                )
                                            ],
                                            p="lg",
                                            radius="md",
                                            withBorder=True,
                                            style={"borderLeft": "4px solid #6366f1"},
                                        )
                                    ],
                                    span={"base": 12, "sm": 6, "md": 3},
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Group(
                                                    [
                                                        dmc.Text("📊", size="xl"),
                                                        dmc.Stack(
                                                            [
                                                                dmc.Text(
                                                                    f"{telemetry_visualizations['total_requests']:,}",
                                                                    size="xl",
                                                                    fw="bold",
                                                                ),
                                                                dmc.Text(
                                                                    "Total Page Views",
                                                                    size="sm",
                                                                    c="dimmed",
                                                                ),
                                                            ],
                                                            gap=0,
                                                        ),
                                                    ],
                                                    align="center",
                                                )
                                            ],
                                            p="lg",
                                            radius="md",
                                            withBorder=True,
                                            style={"borderLeft": "4px solid #10b981"},
                                        )
                                    ],
                                    span={"base": 12, "sm": 6, "md": 3},
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Group(
                                                    [
                                                        dmc.Text("🌍", size="xl"),
                                                        dmc.Stack(
                                                            [
                                                                dmc.Text(
                                                                    f"{len(set(loc.get('country', 'Unknown') for loc in telemetry_visualizations.get('locations', []) if loc.get('country') and loc.get('country') != 'Unknown'))}",
                                                                    size="xl",
                                                                    fw="bold",
                                                                ),
                                                                dmc.Text(
                                                                    "Countries",
                                                                    size="sm",
                                                                    c="dimmed",
                                                                ),
                                                            ],
                                                            gap=0,
                                                        ),
                                                    ],
                                                    align="center",
                                                )
                                            ],
                                            p="lg",
                                            radius="md",
                                            withBorder=True,
                                            style={"borderLeft": "4px solid #f59e0b"},
                                        )
                                    ],
                                    span={"base": 12, "sm": 6, "md": 3},
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Group(
                                                    [
                                                        dmc.Text("⚡", size="xl"),
                                                        dmc.Stack(
                                                            [
                                                                dmc.Text(
                                                                    "Live",
                                                                    size="xl",
                                                                    fw="bold",
                                                                ),
                                                                dmc.Text(
                                                                    "Auto-refresh",
                                                                    size="sm",
                                                                    c="dimmed",
                                                                ),
                                                            ],
                                                            gap=0,
                                                        ),
                                                    ],
                                                    align="center",
                                                )
                                            ],
                                            p="lg",
                                            radius="md",
                                            withBorder=True,
                                            style={"borderLeft": "4px solid #ef4444"},
                                        )
                                    ],
                                    span={"base": 12, "sm": 6, "md": 3},
                                ),
                            ],
                        ),
                        # Date range picker
                        dmc.Paper(
                            [
                                dmc.Group(
                                    [
                                        dmc.DatePicker(
                                            id="date-range-picker",
                                            label="Date Range",
                                            placeholder="Select date range",
                                            type="range",
                                            value=[
                                                (
                                                    datetime.now() - timedelta(days=14)
                                                ).strftime("%Y-%m-%d"),
                                                datetime.now().strftime("%Y-%m-%d"),
                                            ],
                                            style={"flex": 1},
                                        ),
                                        dmc.Button(
                                            "Refresh",
                                            id="refresh-btn",
                                            variant="light",
                                            leftSection="🔄",
                                        ),
                                    ],
                                    align="end",
                                ),
                                dmc.Text(
                                    id="date-range-display",
                                    size="sm",
                                    c="dimmed",
                                    style={"marginTop": "10px"},
                                ),
                            ],
                            p="md",
                            radius="md",
                            withBorder=True,
                            style={"marginTop": "20px"},
                        ),
                        # Charts
                        dmc.Grid(
                            [
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dcc.Graph(
                                                    id="time-series",
                                                    figure=time_series_fig,
                                                    config={
                                                        "displayModeBar": False,
                                                        "responsive": True,
                                                    },
                                                )
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        ),
                                    ],
                                    span={"base": 12, "lg": 8},
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dcc.Graph(
                                                    id="endpoints",
                                                    figure=endpoints_fig,
                                                    config={
                                                        "displayModeBar": False,
                                                        "responsive": True,
                                                    },
                                                )
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        )
                                    ],
                                    span={"base": 12, "lg": 4},
                                ),
                            ]
                        ),
                        dmc.Grid(
                            [
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dcc.Graph(
                                                    id="map",
                                                    figure=map_fig,
                                                    config={
                                                        "displayModeBar": False,
                                                        "responsive": True,
                                                    },
                                                )
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        ),
                                    ],
                                    span=12,
                                ),
                            ]
                        ),
                    ],
                    gap="lg",
                ),
            ],
            size="xl",
        )
        return layout
    except Exception as e:
        logger.error(f"Error creating metrics layout: {str(e)}")
        # Create a more informative error layout
        error_layout = dmc.Container(
            [
                dmc.Alert(
                    title="Error Loading Metrics",
                    children=[
                        dmc.Text(
                            f"An error occurred while loading the metrics dashboard: {str(e)}"
                        ),
                        dmc.Text(
                            "Please check the logs for more details.",
                            size="sm",
                            c="dimmed",
                        ),
                        dmc.Code(
                            traceback.format_exc(),
                            block=True,
                            style={"marginTop": "10px"},
                        ),
                    ],
                    color="red",
                    style={"marginTop": "20px"},
                )
            ],
            size="xl",
        )
        return error_layout


# Make the layout dynamic by using a callback
@callback(Output("metrics-content", "children"), Input("url", "pathname"))
def display_metrics_page(pathname):
    if pathname == "/metrics":
        return get_metrics_layout()
    return no_update


# Callback to display selected date range
@callback(Output("date-range-display", "children"), Input("date-range-picker", "value"))
def update_date_range_display(date_range):
    if date_range:
        if isinstance(date_range, list) and len(date_range) == 2:
            start_date, end_date = date_range
            return f"Selected range: {start_date} to {end_date}"
        else:
            return f"Selected date: {date_range}"
    return "No date range selected"


# Set up the page layout structure
layout = html.Div(id="metrics-content")


@callback(
    [
        Output("time-series", "figure"),
        Output("endpoints", "figure"),
        Output("map", "figure"),
    ],
    [
        Input("telemetry-refresh-interval", "n_intervals"),
        Input("refresh-btn", "n_clicks"),
        Input("date-range-picker", "value"),
    ],
    prevent_initial_call=True,
)
def refresh_telemetry(_, __, date_range):
    try:
        telemetry_visualizations = get_telemetry_visualizations()

        time_series_fig = telemetry_visualizations.get("time_series", go.Figure())
        endpoints_fig = telemetry_visualizations.get("endpoints", go.Figure())
        map_fig = telemetry_visualizations.get("map", go.Figure())

        return time_series_fig, endpoints_fig, map_fig
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return no_update, no_update, no_update
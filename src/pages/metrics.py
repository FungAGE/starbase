import dash
from dash import dcc, html, callback, Input, Output, no_update
import dash_mantine_components as dmc
from src.telemetry.utils import analyze_telemetry
from src.config.logging import get_logger

dash.register_page(__name__)


logger = get_logger(__name__)


def get_metrics_layout():
    try:
        telemetry_data = analyze_telemetry()
        if telemetry_data is None:
            telemetry_data = {"unique_users": 0, "total_requests": 0}

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
                        dmc.Title("Telemetry Dashboard", order=1),
                        dmc.Text(
                            "This dashboard shows visitor statistics and usage patterns.",
                            size="lg",
                        ),
                        dmc.Grid(
                            style={"paddingTop": "20px"},
                            children=[
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Text(
                                                    f"Total Unique Users: {telemetry_data['unique_users']}",
                                                    size="md",
                                                ),
                                                dmc.Text(
                                                    "Number of unique IP addresses across all pages",
                                                    size="sm",
                                                    c="dimmed",
                                                ),
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        )
                                    ],
                                    span={
                                        "base": 12,
                                        "sm": 6,
                                    },
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dmc.Text(
                                                    f"Sum of Unique Visitors per Page: {telemetry_data['total_requests']}",
                                                    size="md",
                                                ),
                                                dmc.Text(
                                                    "Sum of unique visitors counted separately for each page",
                                                    size="sm",
                                                    c="dimmed",
                                                ),
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        )
                                    ],
                                    span={
                                        "base": 12,
                                        "sm": 6,
                                    },
                                ),
                            ],
                        ),
                        dmc.Grid(
                            [
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dcc.Graph(
                                                    id="time-series",
                                                    figure=telemetry_data[
                                                        "time_series"
                                                    ],
                                                )
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        ),
                                    ],
                                    span={
                                        "base": 12,
                                        "lg": 6,
                                    },
                                ),
                                dmc.GridCol(
                                    [
                                        dmc.Paper(
                                            [
                                                dcc.Graph(
                                                    id="endpoints",
                                                    figure=telemetry_data["endpoints"],
                                                )
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        )
                                    ],
                                    span={
                                        "base": 12,
                                        "lg": 6,
                                    },
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
                                                    figure=telemetry_data["map"],
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
                    gap=3,
                ),
            ],
            size="xl",
        )
        return layout
    except Exception as e:
        logger.error(f"Error creating metrics layout: {str(e)}")
        return html.Div("Error loading metrics")


layout = get_metrics_layout()


@callback(
    [
        Output("time-series", "figure"),
        Output("endpoints", "figure"),
        Output("map", "figure"),
    ],
    Input("telemetry-refresh-interval", "n_intervals"),
)
def refresh_telemetry(_):
    try:
        telemetry_data = analyze_telemetry()
        return (
            telemetry_data["time_series"],
            telemetry_data["endpoints"],
            telemetry_data["map"],
        )
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return no_update, no_update, no_update

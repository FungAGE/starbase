import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import logging
import plotly.express as px
import pandas as pd
from src.utils.telemetry import analyze_telemetry

dash.register_page(__name__)

logger = logging.getLogger(__name__)

# create dynamic layout with dash mantine components
def get_metrics_layout():
  try:
    telemetry_data = analyze_telemetry()
    if telemetry_data is None:
      telemetry_data = {
                "unique_users": 0,
                "endpoints": [],
                "total_requests": 0
            }

    layout = dmc.Container([
        dmc.Stack(style={"paddingTop": "20px"}, children=[
            dmc.Title('Telemetry Dashboard', order=1),
            dmc.Text('This dashboard shows the number of unique users and the most accessed endpoints.', size='lg'),
            dmc.Grid(style={"paddingTop": "20px"}, children=[
                dmc.GridCol([
                    dmc.Paper([
                        dmc.Text(f"Unique Users: {telemetry_data['unique_users']}", size="md"),
                    ], p='md', radius='md', withBorder=True)
                ], span={
                    "base": 12,
                    "sm": 6,
                }),
                dmc.GridCol([
                    dmc.Paper([
                        dmc.Text(f"Total Requests: {telemetry_data['total_requests']}", size="md"),
                    ], p='md', radius='md', withBorder=True)
                ], span={
                    "base": 12,
                    "sm": 6,
                }),
            ]),
            dmc.Grid([
                dmc.GridCol([
                    dmc.Paper([
                        dcc.Graph(id='time-series', figure=telemetry_data['time_series'])
                    ], p='md', radius='md', withBorder=True),
                ], span={
                    "base": 12,
                    "lg": 6,
                }),
                dmc.GridCol([
                    dmc.Paper([
                        dcc.Graph(id='endpoints', figure=telemetry_data['endpoints'])
                    ], p='md', radius='md', withBorder=True)
                ], span={
                    "base": 12,
                    "lg": 6,
                }),
            ]),
            dmc.Grid([
                dmc.GridCol([
                    dmc.Paper([
                        dcc.Graph(id='map', figure=telemetry_data['map'])
                    ], p='md', radius='md', withBorder=True),
                ], span=12),
            ]),
        ], gap=3)
    ], size="xl")
    return layout
  except Exception as e:
      logger.error(f"Error creating metrics layout: {str(e)}")
      return html.Div("Error loading metrics")

layout = get_metrics_layout()
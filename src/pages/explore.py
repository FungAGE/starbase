import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, html, dash_table, callback
from dash.dependencies import Input, Output

import plotly.express as px
import pandas as pd

# Sample data for sunburst chart and table
df = pd.DataFrame(
    {
        "category": ["A", "A", "A", "B", "B", "C", "C", "C", "C"],
        "subcategory": ["A1", "A2", "A3", "B1", "B2", "C1", "C2", "C3", "C4"],
        "subsubcategory": [
            "A11",
            "A21",
            "A31",
            "B11",
            "B21",
            "C11",
            "C21",
            "C31",
            "C41",
        ],
        "value": [10, 20, 30, 40, 50, 15, 25, 35, 45],
    }
)

dash.register_page(__name__)

# Layout for Dash app with two sunburst charts and a table
layout = html.Div(
    [
        html.Div(
            [
                dcc.Graph(id="sunburst1"),
                dcc.Graph(id="sunburst2"),
            ],
            style={"display": "flex", "justify-content": "space-around"},
        ),
        # Data Table for showing filtered results
        html.Div(
            [
                dash_table.DataTable(
                    id="table",
                    columns=[
                        {"name": "Category", "id": "category"},
                        {"name": "Subcategory", "id": "subcategory"},
                        {"name": "Subsubcategory", "id": "subsubcategory"},
                        {"name": "Value", "id": "value"},
                    ],
                    data=[],  # The data will be updated dynamically
                    page_size=10,  # Number of rows per page
                )
            ],
            style={
                "margin-top": "20px",
                "width": "80%",
                "margin-left": "auto",
                "margin-right": "auto",
            },
        ),
    ]
)


# Callback to update both charts and the table based on user interaction
@callback(
    [
        Output("sunburst1", "figure"),
        Output("sunburst2", "figure"),
        Output("table", "data"),
    ],
    [Input("sunburst1", "clickData"), Input("sunburst2", "clickData")],
)
def update_sunburst_and_table(sunburst1_click, sunburst2_click):
    # Get the clicked data for crossfiltering
    filter_value = None
    if sunburst1_click:
        filter_value = sunburst1_click["points"][0]["label"]
    elif sunburst2_click:
        filter_value = sunburst2_click["points"][0]["label"]

    # Filter data based on the clicked segment, if available
    if filter_value:
        filtered_df = df[df["category"] == filter_value]
    else:
        filtered_df = df

    # Sunburst chart 1
    fig1 = px.sunburst(
        df,
        path=["category", "subcategory", "subsubcategory"],
        values="value",
        title="Sunburst Chart 1",
    )

    # Sunburst chart 2 (filtered by click)
    fig2 = px.sunburst(
        filtered_df,
        path=["category", "subcategory", "subsubcategory"],
        values="value",
        title="Sunburst Chart 2 (Filtered)",
    )

    # Update the table data (filtering applied)
    table_data = filtered_df.to_dict("records")

    return fig1, fig2, table_data

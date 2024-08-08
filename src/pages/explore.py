import dash
from dash import dash_table, dcc, html, callback, clientside_callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px

from Bio import Phylo
import plotly.graph_objects as go

from src.components.tree import plot_tree

dash.register_page(__name__)

df = pd.read_csv("src/data/joined_ships.csv")
df_sub = df[
    [
        "starshipID",
        "starship_family",
        "starship_navis",
        "starship_haplotype",
        "genus",
        "species",
    ]
]

ship_count = df[["starshipID"]].nunique()
species = df["genus"] + "-" + df["species"]
species_count = species.nunique()


def load_svg(svg_file):
    with open(svg_file, "r") as file:
        svg_content = file.read()
    return svg_content


# Load the SVG content
svg_content = load_svg(
    "/home/adrian/Systematics/Starship_Database/starbase/tmp/test.svg"
)


layout = html.Div(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.Stack(
                    [
                        dbc.Stack(
                            [
                                dbc.Row(
                                    justify="center",
                                    align="start",
                                    children=[
                                        dbc.Col(
                                            lg=3,
                                            sm=6,
                                            children=[
                                                dbc.Card(
                                                    [
                                                        dbc.CardHeader(
                                                            [
                                                                html.H4(
                                                                    html.P(
                                                                        [
                                                                            "Total number of Starships in ",
                                                                            html.Span(
                                                                                "starbase",
                                                                                className="logo-text",
                                                                            ),
                                                                            ":",
                                                                        ]
                                                                    ),
                                                                    className="card-title",
                                                                ),
                                                            ]
                                                        ),
                                                        dbc.CardBody(
                                                            [
                                                                html.H1(ship_count),
                                                            ],
                                                        ),
                                                    ],
                                                    style={
                                                        "display": "flex",
                                                        "textAlign": "center",
                                                        "fontSize": "1vw",
                                                    },
                                                    color="primary",
                                                    inverse=True,
                                                ),
                                            ],
                                        ),
                                        dbc.Col(
                                            lg=3,
                                            sm=8,
                                            children=[
                                                dbc.Card(
                                                    [
                                                        dbc.CardHeader(
                                                            [
                                                                html.H4(
                                                                    "Fungal species with Starships",
                                                                    className="card-title",
                                                                ),
                                                            ]
                                                        ),
                                                        dbc.CardBody(
                                                            [
                                                                html.H1(
                                                                    species_count,
                                                                ),
                                                            ],
                                                        ),
                                                    ],
                                                    style={
                                                        "display": "flex",
                                                        "textAlign": "center",
                                                        "fontSize": "1vw",
                                                    },
                                                    color="secondary",
                                                    inverse=True,
                                                ),
                                            ],
                                        ),
                                    ],
                                    style={"paddingTop": "20px"},
                                ),
                                dcc.Location(id="url", refresh=False),
                                dbc.Row(
                                    justify="center",
                                    align="center",
                                    children=[
                                        dbc.Col(
                                            lg=3,
                                            sm=8,
                                            children=[
                                                dcc.Loading(
                                                    id="loading-2",
                                                    type="default",
                                                    children=[
                                                        dcc.Graph(
                                                            id="pie-chart2",
                                                            config={
                                                                "displayModeBar": False
                                                            },
                                                        ),
                                                    ],
                                                )
                                            ],
                                        ),
                                        dbc.Col(
                                            lg=3,
                                            sm=8,
                                            children=[
                                                dcc.Loading(
                                                    id="loading-1",
                                                    type="default",
                                                    children=[
                                                        dcc.Graph(
                                                            id="pie-chart1",
                                                            config={
                                                                "displayModeBar": False,
                                                            },
                                                        ),
                                                    ],
                                                )
                                            ],
                                        ),
                                    ],
                                ),
                                dbc.Row(
                                    justify="center",
                                    align="center",
                                    children=[
                                        dbc.Col(
                                            lg=6,
                                            sm=10,
                                            children=[
                                                html.H2(
                                                    [
                                                        "All Starships in ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                    ]
                                                ),
                                            ],
                                        ),
                                    ],
                                ),
                                dbc.Row(
                                    justify="center",
                                    align="center",
                                    children=[
                                        dbc.Col(
                                            lg=6,
                                            sm=10,
                                            children=[
                                                dash_table.DataTable(
                                                    id="table",
                                                    columns=[
                                                        {
                                                            "name": i,
                                                            "id": i,
                                                            "deletable": False,
                                                            "selectable": True,
                                                        }
                                                        for i in df_sub.columns
                                                    ],
                                                    data=df_sub.to_dict("records"),
                                                    editable=False,
                                                    filter_action="native",
                                                    sort_action="native",
                                                    sort_mode="multi",
                                                    # column_selectable="single",
                                                    row_selectable="multi",
                                                    row_deletable=False,
                                                    selected_columns=[],
                                                    selected_rows=[],
                                                    page_action="native",
                                                    page_current=0,
                                                    page_size=10,
                                                ),
                                                html.Div(id="table-container"),
                                            ],
                                        ),
                                    ],
                                ),
                            ],
                            gap=4,
                            direction="vertical",
                        ),
                        # html.Div(
                        #     [
                        #         html.H1("Interactable SVG in Dash"),
                        #         html.Div(id="svg-container", children=svg_content),
                        #         html.Div(id="output-container"),
                        #     ]
                        # ),
                        plot_tree(),
                    ],
                    direction="horizontal",
                    gap=1,
                )
            ],
        )
    ]
)


def agg_df(df, groups):
    agg = df.groupby(groups).starshipID.agg(
        ["count", "nunique", lambda x: x.size - x.nunique()]
    )
    agg = agg.reset_index()

    return agg


# Callback to update the sunburst figure based on selected row in the DataTable
@callback(
    [Output("pie-chart1", "figure"), Output("pie-chart2", "figure")],
    [Input("url", "pathname")],
)
def update_sunburst(selected_rows):

    ship_selection = agg_df(df, ["starship_family", "starship_navis"])
    tax_selection = agg_df(df, ["order", "family"])

    ship_pie = px.sunburst(
        ship_selection,
        path=["starship_family", "starship_navis"],
        values="count",
    )
    ship_pie.update_layout(
        autosize=True,
        title_font=dict(size=24),
        title={
            "text": "Starships by captain family",
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(t=50, l=0, r=0, b=0),
    )

    tax_pie = px.sunburst(
        tax_selection,
        path=["order", "family"],
        values="count",
    )

    tax_pie.update_layout(
        autosize=True,
        title_font=dict(size=24),
        title={
            "text": "Starships by Order",
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(t=50, l=0, r=0, b=0),
    )

    return ship_pie, tax_pie


# # Callback to handle node clicks (you need to implement custom JS in the HTML to send node data)
# @callback(Output("node-info", "children"), [Input("tree", "n_clicks")])
# def display_node_info(n_clicks):
#     # Placeholder: Return information about the clicked node
#     # This part requires custom JS in the HTML to capture the node click and send data to Dash
#     if n_clicks:
#         return f"Node clicked {n_clicks} times"
#     return "Click on a node to see details."


# JavaScript to listen for custom events and update Dash Store
# clientside_callback(
#     """
#     function(n_clicks) {
#         document.addEventListener('elementClicked', function(event) {
#             var clickedElement = event.detail.id;
#             var store = document.querySelector('[data-dash-is-loading="clicked-element"]');
#             store.__value = clickedElement;
#             store.dispatchEvent(new Event('change', { bubbles: true }));
#         });
#         return window.dash_clientside.no_update;
#     }
#     """,
#     Output("clicked-element", "data"),
#     Input("svg-container", "n_clicks"),
# )


# # Callback to handle interactions
# @callback(Output("output-container", "children"), Input("clicked-element", "data"))
# def update_output(clicked_element):
#     if clicked_element is None:
#         return "No SVG element clicked yet"
#     return f"{clicked_element} was clicked"

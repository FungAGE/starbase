import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, html, callback, no_update
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import pandas as pd
from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot
from src.utils.tree import plot_tree, hex_to_rgba, default_highlight_colors
from src.components.callbacks import curated_switch

dash.register_page(__name__)

species_card = dbc.Card(
    [
        dbc.CardHeader([html.Div(id="species-card-header")]),
        dbc.CardBody(
            [
                html.H1(
                    html.Div(id="species-count"),
                ),
            ],
        ),
    ],
    style={
        "textAlign": "center",
        "fontSize": "1vw",
        "height": 150,
        "maxWidth": "500px",
    },
    color="secondary",
    inverse=True,
)

tax_card = dbc.Card(
    [
        dbc.CardHeader([html.Div(id="ship-card-header")]),
        dbc.CardBody(
            [
                html.H1(html.Div(id="ship-count")),
            ],
        ),
    ],
    style={
        "textAlign": "center",
        "fontSize": "1vw",
        "height": 150,
        "maxWidth": "500px",
    },
    color="primary",
    inverse=True,
)

layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dmc.Container(
            fluid=True,
            children=[
                dmc.Grid(
                    justify="center",
                    align="center",
                    children=[
                        dmc.GridCol(
                            span={"lg": 6, "sm": 12},
                            children=[
                                dmc.Grid(
                                    justify="center",
                                    align="start",
                                    style={"padding": "10px"},
                                    children=[
                                        dmc.GridCol(
                                            span={"lg": 6, "sm": 10},
                                            children=curated_switch(text="Subset to curated Starships",size="large"),
                                        ),
                                        dmc.GridCol(
                                            span={"lg": 6, "sm": 10},
                                            children=dbc.Button(
                                                "Reset",
                                                id="reset-button",
                                                n_clicks=0,
                                                color="success",
                                                class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
                                            ),
                                        ),
                                    ],
                                ),
                                dmc.Grid(
                                    justify="center",
                                    align="start",
                                    style={"padding": "10px"},
                                    children=[
                                        dmc.GridCol(
                                            span={"lg": 6, "sm": 10},
                                            children=[
                                                dmc.Stack(
                                                    [
                                                        species_card,
                                                        dcc.Loading(
                                                            id="loading-2",
                                                            type="default",
                                                            children=[
                                                                dcc.Graph(
                                                                    id="pie-chart2",
                                                                    config={
                                                                        "displayModeBar": False,
                                                                    },
                                                                    style={
                                                                        "padding": "10px"
                                                                    },
                                                                ),
                                                            ],
                                                        ),
                                                    ]
                                                )
                                            ],
                                        ),
                                        dmc.GridCol(
                                            span={"lg": 6, "sm": 10},
                                            children=[
                                                dmc.Stack(
                                                    [
                                                        tax_card,
                                                        dcc.Graph(
                                                            id="pie-chart1",
                                                            config={
                                                                "displayModeBar": False,
                                                            },
                                                            style={"padding": "10px"},
                                                        ),
                                                    ]
                                                )
                                            ],
                                        ),
                                    ],
                                ),
                                dmc.GridCol(
                                    span=12,
                                    children=[
                                        dcc.Loading(
                                            id="loading-3",
                                            type="default",
                                            children=[
                                                html.Div(
                                                    id="explore-table",
                                                    className="center-content",
                                                ),
                                            ],
                                        ),
                                        # html.Div(
                                        #     id="explore-table-cache-error"
                                        # ),
                                    ],
                                ),
                            ],
                        ),
                        dmc.GridCol(
                            span={"lg": 4, "sm": 8},
                            children=[
                                dcc.Loading(
                                    id="loading-4",
                                    type="default",
                                    children=[
                                        dcc.Graph(
                                            id="explore-phylogeny",
                                            className="div-card",
                                            # figure=plot_tree(highlight_families=None),
                                        ),
                                        # html.Div(id="phylogeny-cache-error"),
                                    ],
                                ),
                            ],
                            className="justify-content-center",
                        ),
                    ],
                ),
            ],
        ),
    ]
)

# @callback(
#     [Output('pie-chart1', 'clickData'),
#     Output('pie-chart2', 'clickData'),],
#     Input('explore-table', 'selected_rows'),
#     State('explore-table', 'data')
# )
# def create_sunburst_plot_click_data(selected_rows, table_data):
#     if not selected_rows:
#         return None

#     clicked_points = []
#     for row_idx in selected_rows:
#         row_data = table_data[row_idx]
#         clicked_points.append({
#             "label": row_data["label"],
#             "parent": row_data["parent"],
#             "value": row_data["value"],
#             "id": row_data["id"]
#         })

#     # Return as Sunburst clickData-like structure
#     return {
#         "points": [{"label": point["label"], "value": point["value"], "id": point["id"]} for point in clicked_points]
#     }


@callback(
    [
        Output("ship-card-header", "children"),
        Output("ship-count", "children"),
        Output("species-card-header", "children"),
        Output("species-count", "children"),
    ],
    [
        Input("curated-status", "data"),
        Input("curated-dataset", "data"),
    ],
)
def make_ship_cards(curated_status, cached_data):
    df = pd.DataFrame(cached_data)
    ship_count = df["starshipID"].nunique()
    species = df["genus"] + "-" + df["species"]
    species_count = species.nunique()

    ship_card_header = html.H4(
        html.P(
            [
                f"Total number {curated_status}of Starships in ",
                html.Span(
                    "starbase",
                    className="logo-text",
                ),
                ":",
            ]
        ),
        className="card-title",
    )
    species_card_header = html.H4(
        f"Fungal species with {curated_status}Starships",
        className="card-title",
    )
    return ship_card_header, ship_count, species_card_header, species_count


def update_phylogeny(phylo, selected_clades):
    print(selected_clades)
    if "shapes" in phylo["layout"]:
        for shape in phylo["layout"]["shapes"]:
            if "fillcolor" in shape:
                for key, val in default_highlight_colors.items():
                    if key in selected_clades:
                        print(key)
                        print(hex_to_rgba(default_highlight_colors[key]))
                        shape["fillcolor"] = hex_to_rgba(default_highlight_colors[key])
                    else:
                        shape["fillcolor"] = "rgba(255, 0, 0, 0)"
    # if "markers" in phylo["layout"]:
    #     for scatter in phylo["layout"]["markers"]:
    #         if (
    #             "name" in scatter
    #             and scatter["name"] in default_highlight_colors.keys()
    #         ):
    #             scatter["marker"]["color"] = default_highlight_colors[
    #                 selected_clade
    #             ]
    if "annotations" in phylo["layout"]:
        for text in phylo["layout"]["annotations"]:
            if text["text"] in selected_clades:
                for selected_clade in selected_clades:
                    text["font_size"] = 24
                else:
                    text["font_size"] = 0

    return phylo


@callback(
    [
        Output("pie-chart1", "figure"),
        Output("pie-chart2", "figure"),
        Output("explore-phylogeny", "figure"),
        Output("explore-table", "children"),
    ],
    [
        Input("reset-button", "n_clicks"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("explore-phylogeny", "clickData"),
        Input("explore-table", "derived_virtual_data"),
        Input("explore-table", "derived_virtual_selected_rows"),
        Input("curated-dataset", "data"),
        Input("phylogeny-cache", "data"),
        Input("pie1-cache", "data"),
        Input("pie2-cache", "data"),
        Input("explore-table-cache", "data"),
    ],
)
def update_ui(
    n_clicks,
    clickData1,
    clickData2,
    clickData_phylo,
    table_data,
    table_rows,
    cached_data,
    phylogeny_cache,
    pie1_cache,
    pie2_cache,
    table_cache,
):
    initial_df = pd.DataFrame(cached_data)
    pie1_output = pie1_cache
    pie2_output = pie2_cache
    tree_output = phylogeny_cache
    table_output = table_cache

    if n_clicks:
        return pie1_output, pie2_output, tree_output, table_output

    # If pie-chart1 was clicked, update pie-chart2, tree, and table
    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        plot_df = initial_df[initial_df["familyName"].isin(path1)]
        pie1_output = no_update
        pie2_output = create_sunburst_plot(plot_df, "tax")
        # tree_output = update_phylogeny(phylogeny_cache, path1)
        tree_output = plot_tree(path1)
        table_output = make_ship_table(df=plot_df, id="explore-table")

    # If pie-chart2 was clicked, update pie-chart1 and tree
    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        plot_df = initial_df[initial_df["order"].isin(path2)]
        pie1_output = create_sunburst_plot(plot_df, "ship")
        pie2_output = no_update
        # tree_output = update_phylogeny(phylogeny_cache, path2)
        tree_output = plot_tree(path2)
        table_output = make_ship_table(df=plot_df, id="explore-table")

    # If phylogeny was clicked, update pie-chart1, pie-chart2, and table
    if clickData_phylo:
        selected_clades = [clickData_phylo["points"][0]["text"]]
        plot_df = initial_df[initial_df["familyName"].isin(selected_clades)]
        pie1_output = create_sunburst_plot(plot_df, "ship")
        pie2_output = create_sunburst_plot(plot_df, "tax")
        tree_output = no_update
        table_output = make_ship_table(df=plot_df, id="explore-table")

    # If table rows were selected, update pie-chart1, pie-chart2, and tree
    if table_rows:
        table_df = pd.DataFrame(table_data).iloc[table_rows]
        selected_clades = table_df["familyName"].tolist()
        pie1_output = create_sunburst_plot(table_df, "ship")
        pie2_output = create_sunburst_plot(table_df, "tax")
        # tree_output = update_phylogeny(phylogeny_cache, selected_clades)
        tree_output = plot_tree(selected_clades)
        table_output = no_update

    # Return the updated or cached figures and table
    return pie1_output, pie2_output, tree_output, table_output

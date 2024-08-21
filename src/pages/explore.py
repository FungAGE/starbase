import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, html, callback
from dash.dependencies import Output, Input
import dash_bootstrap_components as dbc

import pandas as pd
from src.utils.tree import plot_tree, tree_file, metadata, default_highlight_clades
from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot

dash.register_page(__name__)

specified_columns = [
    "starshipID",
    "starship_family",
    "starship_navis",
    "starship_haplotype",
    "genus",
    "species",
]

ship_groups = ["starship_family", "starship_navis"]
ship_title = "Starships by Family/Navis"
tax_groups = ["order", "family"]
tax_title = "Starships by Order/Family"

curated_switch = dbc.Row(
    justify="center",
    align="start",
    style={"paddingTop": "20px"},
    children=[
        dbc.Col(
            lg=6,
            sm=8,
            children=[
                dbc.RadioItems(
                    options=[
                        {"label": "Use just curated Starships", "value": 1},
                    ],
                    value=[],
                    id="curated-input",
                    switch=True,
                    inline=True,
                    label_style={
                        "fontSize": "24px",
                    },  # Controls label size
                    style={
                        "fontSize": "20px",
                        "transform": "scale(1.5)",
                    },  # Controls radio button size
                ),
            ],
        )
    ],
)

species_card = (
    dbc.Card(
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
        },
        color="secondary",
        inverse=True,
    ),
)

tax_card = (
    dbc.Card(
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
        },
        color="primary",
        inverse=True,
    ),
)

layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dcc.Store("curated-dataset"),
        # dcc.Store(id="phylogeny-cache"),
        # dcc.Store(id="pie-chart-cache"),
        dbc.Container(
            fluid=True,
            children=[
                dbc.Row(
                    justify="center",
                    align="start",
                    style={"paddingTop": "20px"},
                    children=[
                        dbc.Col(
                            lg=6,
                            sm=8,
                            children=[
                                dbc.Stack(
                                    [
                                        curated_switch,
                                        dbc.Row(
                                            justify="center",
                                            align="start",
                                            children=[
                                                dbc.Col(
                                                    lg=6,
                                                    sm=10,
                                                    children=species_card,
                                                ),
                                                dbc.Col(
                                                    lg=6,
                                                    sm=10,
                                                    children=tax_card,
                                                ),
                                            ],
                                        ),
                                        dbc.Row(
                                            justify="center",
                                            align="center",
                                            children=[
                                                dbc.Col(
                                                    lg=6,
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
                                                                # html.Div(
                                                                #     id="pie-cache-error"
                                                                # ),
                                                            ],
                                                        )
                                                    ],
                                                ),
                                                dbc.Col(
                                                    lg=6,
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
                                                                # html.Div(
                                                                #     id="pie-cache-error"
                                                                # ),
                                                            ],
                                                        )
                                                    ],
                                                ),
                                            ],
                                        ),
                                        dbc.Button(
                                            "Reset",
                                            id="reset-button",
                                            n_clicks=0,
                                            className="d-grid gap-2 col-3 mx-auto",
                                        ),
                                        dbc.Row(
                                            justify="center",
                                            align="center",
                                            children=[
                                                dbc.Col(
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
                                                        html.Div(id="explore-table"),
                                                        # html.Div(
                                                        #     id="table-cache-error"
                                                        # ),
                                                    ],
                                                ),
                                            ],
                                        ),
                                    ],
                                    gap=4,
                                    direction="vertical",
                                ),
                            ],
                        ),
                        dbc.Col(
                            lg=4,
                            sm=8,
                            children=[
                                dcc.Loading(
                                    id="loading-3",
                                    type="default",
                                    children=[
                                        dcc.Graph(
                                            id="explore-phylogeny",
                                            className="div-card",
                                            # figure=plot_tree(tree_file, metadata, None),
                                        ),
                                        # html.Div(id="phylogeny-cache-error"),
                                    ],
                                ),
                            ],
                            className="justify-content-center",
                        ),
                    ],
                )
            ],
        ),
    ]
)


@callback(
    [
        Output("pie-chart1", "figure"),
        Output("pie-chart2", "figure"),
    ],
    [
        Input("curated-dataset", "data"),
        Input("reset-button", "n_clicks"),
        Input("explore-phylogeny", "clickData"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("explore-table", "derived_virtual_data"),
        Input("explore-table", "derived_virtual_selected_rows"),
    ],
)
def update_sunburst(
    cached_data,
    n_clicks,
    phylo_clickData,
    selected1,
    selected2,
    table_data,
    selected_rows,
):
    initial_df = pd.DataFrame(cached_data)
    plot_df = initial_df

    if phylo_clickData:
        selected_clades = [phylo_clickData["points"][0]["text"]]
        plot_df = initial_df[initial_df["starship_family"].isin(selected_clades)]

    if selected1:
        path1 = [point["label"] for point in selected1["points"]]
        plot_df = plot_df[plot_df[ship_groups[0]].isin(path1)]
    if selected2:
        path2 = [point["label"] for point in selected2["points"]]
        plot_df = plot_df[plot_df[tax_groups[0]].isin(path2)]

    if selected_rows:
        plot_df = pd.DataFrame(table_data).iloc[selected_rows]

    if n_clicks:
        plot_df = initial_df

    ship_pie = create_sunburst_plot(plot_df, ship_groups, ship_title)
    tax_pie = create_sunburst_plot(plot_df, tax_groups, tax_title)

    return ship_pie, tax_pie


@callback(
    Output("explore-phylogeny", "figure"),
    [
        Input("curated-dataset", "data"),
        Input("reset-button", "n_clicks"),
        Input("explore-phylogeny", "clickData"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("explore-table", "derived_virtual_data"),
        Input("explore-table", "derived_virtual_selected_rows"),
    ],
)
def update_phylogeny_tree(
    cached_data,
    n_clicks,
    phylo_clickData,
    clickData1,
    clickData2,
    table_data,
    selected_rows,
):
    initial_df = pd.DataFrame(cached_data)
    selected_clades = []

    if n_clicks:
        return plot_tree(tree_file, metadata, None)

    if phylo_clickData:
        selected_clades = [phylo_clickData["points"][0]["text"]]

    if selected_rows:
        table_df = pd.DataFrame(table_data)
        rows = table_df.iloc[selected_rows]
        for row in rows["starship_family"]:
            selected_clades.append(row)

    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        selection_df = initial_df[initial_df[ship_groups[0]].isin(path1)]
        selected_clades = selection_df["starship_family"].unique()

    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        selection_df = initial_df[initial_df[tax_groups[0]].isin(path2)]
        selected_clades = selection_df["starship_family"].unique()
    print(selected_clades)
    fig = plot_tree(tree_file, metadata, selected_clades)
    return fig


@callback(
    Output("explore-table", "children"),
    [
        Input("curated-dataset", "data"),
        Input("reset-button", "n_clicks"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
    ],
)
def update_table(cached_data, n_clicks, clickData1, clickData2):
    initial_df = pd.DataFrame(cached_data)
    filtered_df = initial_df

    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        filtered_df = initial_df[initial_df[ship_groups[0]].isin(path1)]
    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        filtered_df = initial_df[initial_df[tax_groups[0]].isin(path2)]
    if n_clicks:
        filtered_df = initial_df

    return make_ship_table(filtered_df, "table", specified_columns)

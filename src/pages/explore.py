import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px

from src.components.tree import plot_tree, tree_file, metadata, default_highlight_clades
from src.components.shipTable import make_table
from src.data.joined_ships import df

dash.register_page(__name__)

df_sub = df[
    [
        "starshipID",
        "captain_superfamily",
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

layout = html.Div(
    [
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
                                        dbc.Row(
                                            justify="center",
                                            align="start",
                                            children=[
                                                dbc.Col(
                                                    lg=6,
                                                    sm=10,
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
                                                                        html.H1(
                                                                            ship_count
                                                                        ),
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
                                                    ],
                                                ),
                                                dbc.Col(
                                                    lg=6,
                                                    sm=10,
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
                                                                "textAlign": "center",
                                                                "fontSize": "1vw",
                                                                "height": 150,
                                                            },
                                                            color="secondary",
                                                            inverse=True,
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
                                                        make_table(df_sub),
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
                                            id="phylogeny",
                                            className="div-card",
                                            # figure=plot_tree(tree_file, metadata, None),
                                        ),
                                    ],
                                ),
                            ],
                            className="justify-content-center",
                        ),
                    ],
                )
            ],
        ),
        dcc.Store(
            id="initial-data",
            data=df.to_dict("records"),
        ),
        dcc.Store(id="phylogeny-cache"),
        dcc.Store(id="pie-chart-cache"),
    ]
)


def agg_df(df, groups):
    agg = df.groupby(groups).starshipID.agg(
        count="count", nunique="nunique", duplicates=lambda x: x.size - x.nunique()
    )

    agg = agg.reset_index()

    return agg


def pie_plot(df, path, title):
    selection = agg_df(df, path)

    pie = px.sunburst(
        selection,
        path=path,
        values="count",
    )

    pie.update_layout(
        autosize=True,
        title_font=dict(size=24),
        title={
            "text": title,
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(t=50, l=0, r=0, b=0),
    )
    return pie


ship_groups = ["starship_family", "starship_navis"]
ship_title = "Starships by Superfamily/Navis"
tax_groups = ["order", "family"]
tax_title = "Starships by Order/Family"


@callback(
    [
        Output("initial-phylogeny", "data"),
    ],
    Input("reset-button", "n_clicks"),
)
def cache_phylogeny(n_clicks):
    tree = plot_tree(tree_file, metadata)
    return tree


@callback(
    Output("pie-chart-cache", "data"),  # Cache pie charts data
    [
        Input("reset-button", "n_clicks"),
        Input("phylogeny", "clickData"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("table", "derived_virtual_data"),
        Input("table", "derived_virtual_selected_rows"),
    ],
    [State("initial-data", "data")],
)
def update_sunburst(
    n_clicks,
    phylo_clickData,
    selected1,
    selected2,
    table_data,
    selected_rows,
    initial_data,
):
    df_initial = pd.DataFrame(initial_data)
    plot_df = df_initial

    if phylo_clickData:
        selected_clades = [phylo_clickData["points"][0]["text"]]
        plot_df = df[df["captain_superfamily"].isin(selected_clades)]

    if selected1:
        path1 = [point["label"] for point in selected1["points"]]
        plot_df = plot_df[plot_df[ship_groups[0]].isin(path1)]
    if selected2:
        path2 = [point["label"] for point in selected2["points"]]
        plot_df = plot_df[plot_df[tax_groups[0]].isin(path2)]

    if selected_rows:
        plot_df = pd.DataFrame(table_data).iloc[selected_rows]

    if n_clicks:
        plot_df = df

    ship_pie = pie_plot(plot_df, ship_groups, ship_title)
    tax_pie = pie_plot(plot_df, tax_groups, tax_title)

    # Cache both pie charts
    return {"ship_pie": ship_pie, "tax_pie": tax_pie}


@callback(
    [Output("pie-chart1", "figure"), Output("pie-chart2", "figure")],
    [Input("pie-chart-cache", "data")],
)
def display_cached_pie_charts(cached_data):
    if cached_data is not None:
        return cached_data["ship_pie"], cached_data["tax_pie"]
    return {}, {}  # Placeholder if no data yet


@callback(
    Output("phylogeny-cache", "data"),
    [
        Input("reset-button", "n_clicks"),
        Input("phylogeny", "clickData"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("table", "derived_virtual_data"),
        Input("table", "derived_virtual_selected_rows"),
    ],
    [State("initial-data", "data")],
)
def update_phylogeny_tree(
    n_clicks,
    phylo_clickData,
    clickData1,
    clickData2,
    table_data,
    selected_rows,
    initial_data,
):
    df_initial = pd.DataFrame(initial_data)
    selected_clades = default_highlight_clades

    if n_clicks:
        return plot_tree(tree_file, metadata)

    if phylo_clickData:
        selected_clades.append(phylo_clickData["points"][0]["text"])

    if selected_rows:
        df = pd.DataFrame(table_data)
        rows = df.iloc[selected_rows]
        for row in rows["captain_superfamily"]:
            selected_clades.append(row)

    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        selection_df = df_initial[df_initial[ship_groups[0]].isin(path1)]
        selected_clades.extend(selection_df["captain_superfamily"].unique())

    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        selection_df = df_initial[df_initial[tax_groups[0]].isin(path2)]
        selected_clades.extend(selection_df["captain_superfamily"].unique())

    fig = plot_tree(tree_file, metadata, selected_clades)
    return fig


@callback(Output("phylogeny", "figure"), [Input("phylogeny-cache", "data")])
def display_cached_phylogeny(cached_fig):
    if cached_fig is not None:
        return cached_fig
    return {}


@callback(
    Output("table", "data"),
    [
        Input("reset-button", "n_clicks"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
    ],
    [State("initial-data", "data")],
)
def update_table(n_clicks, clickData1, clickData2, initial_data):
    df_initial = pd.DataFrame(initial_data)
    filtered_df = df_initial

    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        filtered_df = filtered_df[filtered_df[ship_groups[0]].isin(path1)]
    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        filtered_df = filtered_df[filtered_df[tax_groups[0]].isin(path2)]
    if n_clicks:
        filtered_df = df

    return filtered_df.to_dict("records")

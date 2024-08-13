import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px

from src.components.tree import plot_tree, tree_file, metadata

dash.register_page(__name__)

df = pd.read_csv("src/data/joined_ships.csv")
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
                                                    lg=4,
                                                    sm=8,
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
                                                    lg=4,
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
                                        dbc.Row(
                                            justify="center",
                                            align="center",
                                            children=[
                                                dbc.Col(
                                                    lg=4,
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
                                                    lg=4,
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
                                                            data=df_sub.to_dict(
                                                                "records"
                                                            ),
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
                                            id="phylogeny-graph",
                                            className="div-card",
                                            # figure=plot_tree(tree_file, metadata, None),
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                )
            ],
        )
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
    [Output("pie-chart1", "figure"), Output("pie-chart2", "figure")],
    [
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("table", "derived_virtual_data"),
        Input("table", "derived_virtual_selected_rows"),
    ],
)
def update_sunburst(selected1, selected2, table_data, selected_rows):
    plot_df = df

    if selected1 or selected2:
        if selected1:
            path1 = [point["label"] for point in selected1["points"]]
            plot_df = df[df[ship_groups[0]].isin(path1)]
        if selected2:
            path2 = [point["label"] for point in selected2["points"]]
            plot_df = df[df[tax_groups[0]].isin(path2)]
    elif selected_rows:
        plot_df = pd.DataFrame(table_data).iloc[selected_rows]

    ship_pie = pie_plot(plot_df, ship_groups, ship_title)
    tax_pie = pie_plot(plot_df, tax_groups, tax_title)

    return ship_pie, tax_pie


@callback(
    Output("phylogeny-graph", "figure"),
    [
        Input("phylogeny-graph", "clickData"),
        Input("phylogeny-graph", "figure"),
        Input("pie-chart1", "clickData"),
        Input("pie-chart2", "clickData"),
        Input("table", "derived_virtual_data"),
        Input("table", "derived_virtual_selected_rows"),
    ],
)
def update_phylogeny_tree(
    clickData, fig, selected1, selected2, table_data, selected_rows
):
    original_fig = plot_tree(tree_file, metadata, [])

    if not selected_rows and clickData is None:
        return original_fig

    selected_clades = []
    if clickData is not None:
        selected_clades.append(clickData["points"][0]["text"])

    if selected_rows is not None:
        df = pd.DataFrame(table_data)
        rows = df.iloc[selected_rows]
        for row in rows["captain_superfamily"]:
            selected_clades.append(row)

    if selected1 or selected2:
        if selected1:
            path1 = [point["label"] for point in selected1["points"]]
            selection_df = df[df[ship_groups[0]].isin(path1)]
        if selected2:
            path2 = [point["label"] for point in selected2["points"]]
            selection_df = df[df[tax_groups[0]].isin(path2)]
        for row in selection_df["captain_superfamily"]:
            selected_clades.append(row)

    fig = plot_tree(tree_file, metadata, selected_clades)

    return fig

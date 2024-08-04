import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px

from Bio import Phylo
import plotly.graph_objects as go

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

layout = html.Div(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.Stack(
                    [
                        dbc.Row(
                            justify="center",
                            align="start",
                            children=[
                                dbc.Col(
                                    lg=3,
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
                                                    config={"displayModeBar": False},
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
            ],
        ),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(
                    width=4,
                    children=[
                        dcc.Loading(
                            id="loading-2",
                            type="default",
                            children=dcc.Graph(
                                id="pie-chart2",
                                config={"displayModeBar": False},
                            ),
                        ),
                    ],
                ),
                dbc.Col(
                    width=4,
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="default",
                            children=dcc.Graph(
                                id="pie-chart1",
                                config={"displayModeBar": False},
                            ),
                        ),
                    ],
                ),
            ],
        ),
        html.Br(),
        dbc.Row(
            children=[
                dbc.Col(
                    width=8,
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
                            page_size=25,
                        ),
                        html.Div(id="table-container"),
                    ],
                )
            ],
        ),
        dbc.Row(
            dbc.Col(
                [
                    html.Div(
                        [
                            html.Img(
                                src="assets/images/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.new_colored.treefile.png",
                                width="50%",
                            )
                            # dcc.Graph(id="tree-fig"),
                        ]
                    )
                ]
            )
        ),
    ],
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

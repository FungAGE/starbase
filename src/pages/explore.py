import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input
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

ship_agg = df.groupby(["starship_family", "starship_navis"]).starshipID.agg(
    ["count", "nunique", lambda x: x.size - x.nunique()]
)
ship_agg = ship_agg.reset_index()

tax_agg = df.groupby(["order", "family"]).starshipID.agg(
    ["count", "nunique", lambda x: x.size - x.nunique()]
)
tax_agg = tax_agg.reset_index()

styles = {"pre": {"border": "thin lightgrey solid", "overflowX": "scroll"}}

# Read the tree from file
tree = Phylo.read("tmp/microsynt.treefile", "newick")


# Convert the tree into a Plotly-compatible format
def plotly_tree(tree):
    layout = go.Layout(
        title="Phylogenetic Tree of Captain Superfamilies",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    fig = go.Figure(layout=layout)

    def add_node(node):
        if isinstance(node, Phylo.BaseTree.Clade):
            children = node.clades
            if children:
                for child in children:
                    add_node(child)
        else:
            fig.add_trace(
                go.Scatter(
                    x=[node.branch_length, node.branch_length],
                    y=[node.clades[0].y, node.clades[1].y],
                    mode="lines",
                    line=dict(width=1),
                )
            )

    add_node(tree.root)
    fig.update_layout(showlegend=False)
    return fig


layout = dbc.Container(
            fluid=True,
            className="justify-content-start",
            children=[
                dbc.Row(
                    children=[
                        dbc.Col(
                            width=8,
                            children=[
                                dbc.Stack(
                                    [
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
                                                        )
                                                    ]
                                                ),
                                                dbc.CardBody(
                                                    [
                                                        html.H1(
                                                            ship_count,
                                                            className="card-text",
                                                            style={
                                                                "textAlign": "center"
                                                            },
                                                        ),
                                                    ]
                                                ),
                                            ],
                                            style={"font-size": "1vw", "width": "50%"},
                                            color="primary",
                                            inverse=True,
                                        ),
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H4(
                                                        "Fungal species with Starships:",
                                                        className="card-title",
                                                    ),
                                                ),
                                                dbc.CardBody(
                                                    html.H1(
                                                        species_count,
                                                        className="card-text",
                                                        style={"textAlign": "center"},
                                                    ),
                                                ),
                                            ],
                                            style={"font-size": "1vw", "width": "50%"},
                                            color="secondary",
                                            inverse=True,
                                        ),
                                    ],
                                    gap=3,
                                    direction="horizontal",
                                )
                            ],
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

@callback(
    [
        Output("pie-chart1", "figure"),
        Output("pie-chart2", "figure"),
    ],
    [Input("table", "derived_virtual_selected_rows")],
)
def update_sunburst(selected_rows):
    if selected_rows is None or len(selected_rows) == 0:
        ship_agg_filt = ship_agg
        tax_agg_filt = tax_agg
    else:
        # TODO: filter by specific column, not index
        ship_agg_category = ship_agg.iloc[selected_rows[0]]["starship_family"]
        tax_agg_category = tax_agg.iloc[selected_rows[0]]["order"]

        ship_agg_filt = ship_agg[ship_agg["starship_family"] == ship_agg_category]
        tax_agg_filt = tax_agg[tax_agg["order"] == tax_agg_category]

    ship_pie = px.sunburst(
        ship_agg_filt,
        path=["starship_family", "starship_navis"],
        values="count",
    )
    ship_pie.update_layout(
        width=800,  # Increase the width of the plot
        height=600,  # Increase the height of the plot
        title_font=dict(size=24),  # Increase the text size of the title
        title={
            "text": "Starships by captain family",
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(l=40, r=40, t=40, b=40),  # Add margins to the plot
    )

    tax_pie = px.sunburst(
        tax_agg_filt,
        path=["order", "family"],
        values="count",
    )

    tax_pie.update_layout(
        width=800,  # Increase the width of the plot
        height=600,  # Increase the height of the plot
        title_font=dict(size=24),  # Increase the text size of the title
        title={
            "text": "Starships by Order",
            "y": 1,
            "x": 0.5,
            "xanchor": "center",
            "yanchor": "top",
        },
        margin=dict(l=40, r=40, t=40, b=40),  # Add margins to the plot
    )

    return ship_pie, tax_pie


# Callback to highlight selected part of the sunburst figure based on selected row in the DataTable
@callback(
    Output("pie-chart1", "selectedData"),
    [Input("table", "derived_virtual_selected_rows")],
)
def update_sunburst_selection(selected_rows):
    if selected_rows is None or len(selected_rows) == 0:
        return None

    ship_agg_category = ship_agg.iloc[selected_rows[0]]["starship_family"]
    tax_agg_category = tax_agg.iloc[selected_rows[0]]["order"]

    ship_agg_filt = ship_agg[ship_agg["starship_family"] == ship_agg_category]
    tax_agg_filt = tax_agg[tax_agg["order"] == tax_agg_category]

    ship_agg_selected_subcategories = ship_agg_filt["starship_family"].tolist()
    # tax_agg_selected_subcategories = tax_agg_filt["Subcategory"].tolist()
    return [
        {
            "points": [
                {"label": subcategory}
                for subcategory in ship_agg_selected_subcategories
            ]
        }
    ]

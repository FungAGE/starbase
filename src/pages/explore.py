import dash
from dash import dash_table, dcc, html, callback
from dash.dependencies import Output, Input
import dash_bootstrap_components as dbc

import pandas as pd
import plotly.express as px

dash.register_page(__name__)

df = pd.read_csv("src/assets/joined_ships.csv")
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

layout = html.Div(
    [
        html.Div(
            className="title-bar",
            children=[
                html.H1(
                    "Explore Starship Diversity",
                    className="title-text",
                    style={"text-align": "center"},
                ),
            ],
        ),
        html.Div(
            style={
                "height": "150vh",
                "display": "flex",
                "flex-direction": "column",
                "justify-content": "center",
                "align-items": "center",
            },
            children=[
                dbc.Card(
                    [
                        dbc.CardBody(
                            [
                                html.Table(
                                    style={"width": "50%"},
                                    children=[
                                        html.Tr(
                                            [
                                                html.Td(
                                                    style={"width": "50%"},
                                                    children=[
                                                        dbc.Card(
                                                            dbc.CardBody(
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
                                                                    html.P(
                                                                        ship_count,
                                                                        className="card-text",
                                                                    ),
                                                                ]
                                                            ),
                                                            style={
                                                                "font-size": "26px",
                                                                "width": "75%",
                                                                "justify-content": "center",
                                                                "align-items": "center",
                                                            },
                                                            color="primary",
                                                            inverse=True,
                                                        )
                                                    ],
                                                ),
                                                html.Td(
                                                    style={"width": "50%"},
                                                    children=[
                                                        dbc.Card(
                                                            dbc.CardBody(
                                                                [
                                                                    html.H4(
                                                                        "Fungal species with Starships:",
                                                                        className="card-title",
                                                                    ),
                                                                    html.P(
                                                                        species_count,
                                                                        className="card-text",
                                                                    ),
                                                                ]
                                                            ),
                                                            style={
                                                                "font-size": "26px",
                                                                "width": "75%",
                                                                "justify-content": "center",
                                                                "align-items": "center",
                                                            },
                                                            color="secondary",
                                                            inverse=True,
                                                        )
                                                    ],
                                                ),
                                            ]
                                        ),
                                        html.Br(),
                                        html.Tr(
                                            [
                                                html.Td(
                                                    style={
                                                        "width": "50%",
                                                        "justify-content": "center",
                                                        "align-items": "center",
                                                    },
                                                    children=[
                                                        dcc.Graph(
                                                            id="pie-chart2",
                                                            config={
                                                                "displayModeBar": False
                                                            },
                                                        )
                                                    ],
                                                ),
                                                html.Td(
                                                    style={
                                                        "width": "50%",
                                                        "justify-content": "center",
                                                        "align-items": "center",
                                                    },
                                                    children=[
                                                        dcc.Graph(
                                                            id="pie-chart1",
                                                            config={
                                                                "displayModeBar": False
                                                            },
                                                        )
                                                    ],
                                                ),
                                            ]
                                        ),
                                    ],
                                ),
                                html.Br(),
                                html.Div(
                                    [
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
                                    style={
                                        "textAlign": "center",
                                        "width": "100%",
                                        "display": "inline-block",
                                    },
                                ),
                            ]
                        )
                    ]
                ),
            ],
        ),
    ]
)


# Callback to update the sunburst figure based on selected row in the DataTable
@callback(
    [Output("pie-chart1", "figure"), Output("pie-chart2", "figure")],
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
        title="Starships in starbase by captain family",
    )
    ship_pie.update_layout(
        width=800,  # Increase the width of the plot
        height=600,  # Increase the height of the plot
        title_font=dict(size=24),  # Increase the text size of the title
        margin=dict(l=40, r=40, t=40, b=40),  # Add margins to the plot
    )

    tax_pie = px.sunburst(
        tax_agg_filt,
        path=["order", "family"],
        values="count",
        title="Starships in starbase by Order",
    )

    tax_pie.update_layout(
        width=800,  # Increase the width of the plot
        height=600,  # Increase the height of the plot
        title_font=dict(size=24),  # Increase the text size of the title
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

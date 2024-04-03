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
        html.H1("Explore Starship Diversity"),
        html.Div(
            style={
                "display": "flex",
                "justify-content": "center",
                "align-items": "center",
            },
            children=[
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
                                            style={"font-size": "26px"},
                                            color="primary",
                                            inverse=True,
                                        )
                                    ],
                                ),
                                html.Td(style={"width": "25%"}),
                                html.Td(
                                    style={"width": "50%"},
                                    children=[
                                        dbc.Card(
                                            dbc.CardBody(
                                                [
                                                    html.H4(
                                                        "Fungal species with Starships",
                                                        className="card-title",
                                                    ),
                                                    html.P(
                                                        species_count,
                                                        className="card-text",
                                                    ),
                                                ]
                                            ),
                                            style={"font-size": "26px"},
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
                                    style={"width": "50%"},
                                    children=[
                                        dcc.Graph(
                                            id="pie-chart2",
                                            config={"displayModeBar": False},
                                        )
                                    ],
                                ),
                                html.Td(
                                    style={"width": "25%"},
                                ),
                                html.Td(
                                    style={"width": "50%"},
                                    children=[
                                        dcc.Graph(
                                            id="pie-chart1",
                                            config={"displayModeBar": False},
                                        )
                                    ],
                                ),
                            ]
                        ),
                    ],
                )
            ],
        ),
        html.Br(),
        html.Div(
            style={
                "textAlign": "center",
                "width": "75%",
            },
            children=[
                html.H3("Table for metadata of all Starships in starbase"),
                html.Div(
                    style={
                        "textAlign": "center",
                        "width": "75%",
                        "display": "inline-block",
                    },
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
                            column_selectable="single",
                            row_selectable="multi",
                            row_deletable=False,
                            selected_columns=[],
                            selected_rows=[],
                            page_action="native",
                            page_current=0,
                            page_size=25,
                        ),
                    ],
                ),
                html.Div(id="table-container"),
            ],
        ),
    ]
)


@callback(
    [
        Output("pie-chart1", "figure"),
        Output("pie-chart2", "figure"),
        Output("table", "selected_rows"),
    ],
    [Input("table", "selected_rows")],
)
def update_plots(selected_rows):
    if selected_rows is None or len(selected_rows) == 0:
        ship_agg_filt = ship_agg
        tax_agg_filt = tax_agg
    else:
        # TODO: filter by specific column, not index
        ship_agg_filt = ship_agg.iloc[selected_rows]
        tax_agg_filt = tax_agg.iloc[selected_rows]

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

    return ship_pie, tax_pie, selected_rows


# this callback defines 3 figures
# as a function of the intersection of their 3 selections
# @callback(
#     Output("g1", "figure"),
#     Output("g2", "figure"),
#     Input("g1", "selectedData"),
#     Input("g2", "selectedData"),
# )
# def callback(selection1, selection2):
#     selectedpoints = df.index
#     for selected_data in [selection1, selection2]:
#         if selected_data and selected_data["points"]:
#             selectedpoints = np.intersect1d(
#                 selectedpoints, [p["customdata"] for p in selected_data["points"]]
#             )

#     return [
#         get_figure(df, "Col 1", "Col 2", selectedpoints, selection1),
#         get_figure(df, "Col 3", "Col 4", selectedpoints, selection2),
#     ]


# layout = html.Div([
#     dbc.Container([
#         dbc.Row([
#             dbc.Col([dcc.Graph(figure=ship_pie)]),
#             dbc.Col([dcc.Graph(figure=tax_pie)])
#         ])
#     ]),
# ])

# @callback(
#     Output('table', 'style_data_conditional'),
#     Input('table', 'selected_columns')
# )
# def update_styles(selected_columns):
#     return [{
#         'if': { 'column_id': i },
#         'background_color': '#D2F3FF'
#     } for i in selected_columns]

# @callback(
#     Output('table-container', "children"),
#     Input('table', "derived_virtual_data"),
#     Input('table', "derived_virtual_selected_rows"))

# # TODO: make graph dynamic to selection
# def update_graphs(rows, derived_virtual_selected_rows):
#     # When the table is first rendered, `derived_virtual_data` and
#     # `derived_virtual_selected_rows` will be `None`. This is due to an
#     # idiosyncrasy in Dash (unsupplied properties are always None and Dash
#     # calls the dependent callbacks when the component is first rendered).
#     # So, if `rows` is `None`, then the component was just rendered
#     # and its value will be the same as the component's dataframe.
#     # Instead of setting `None` in here, you could also set
#     # `derived_virtual_data=df.to_rows('dict')` when you initialize
#     # the component.

#     if derived_virtual_selected_rows is None:
#         derived_virtual_selected_rows = []

#     dff = df if rows is None else pd.DataFrame(rows)

#     if derived_virtual_selected_rows:
#         dff_selected = dff.iloc[derived_virtual_selected_rows]
#     else:
#         dff_selected = dff

#     graphs = []
#     for column in ["starship_family"]:
#         if column in dff_selected:
#             dff_agg = dff_selected.groupby([column]).starshipID.agg(['count', 'nunique', lambda x: x.size - x.nunique()])
#             dff_agg = dff_agg.reset_index()

#             graph = dcc.Graph(
#                 id=column,
#                 figure={
#                     "data": [
#                         {
#                             "x": dff_agg[column],
#                             "y": dff_agg["count"],
#                             "type": "bar"
#                         }
#                     ],
#                     "layout": {
#                         "xaxis": {"title": column},
#                         "yaxis": {"title": "count"},
#                         "height": 250,
#                         "margin": {"t": 10, "l": 10, "r": 10},
#                     },
#                 },
#             )
#             graphs.append(graph)

#     return graphs

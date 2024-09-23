import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

import dash
from dash import dcc, callback, html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

import base64
import pandas as pd

from src.components.caching import cache
from src.components.sqlite import engine

from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot
from src.utils.tree import plot_tree


download_ships_button = dbc.Button(
    html.Div(
        [
            "Download Starships from the latest version of ",
            html.Span(
                "starbase",
                className="logo-text",
            ),
            ".",
        ],
    ),
    href="/dl",
    color="primary",
    class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
)

download_ships_card = dbc.Card(
    [
        dbc.CardHeader(
            [
                html.Div(
                    "Data Availability",
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                )
            ],
            className="card-header-custom",
        ),
        dbc.CardBody(
            [
                dmc.Text(
                    [
                        "We have been maintaining ",
                        html.Span(
                            "starbase",
                            className="logo-text",
                        ),
                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export",
                    ],
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                    style={"paddingBottom": "20px"},
                ),
                dmc.Center(download_ships_button),
            ],
        ),
    ],
    className="auto-resize-750",
    # style={"height": "350px"},
    #
)

download_starbase_button = html.Div(
    [
        dbc.Button(
            html.P(
                [
                    "Download Starships from the latest version of ",
                    html.Span(
                        "starbase",
                        className="logo-text",
                    ),
                    ".",
                ]
            ),
            id="dl-button",
            color="primary",
            className="mt-2",
        ),
    ],
    className="text-center",
    style={
        "font-size": "0.875rem",
    },
)


def curated_switch(text, size="normal"):
    if size == "normal":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "justify-content": "center",
        }
    if size == "large":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "transform": "scale(1.5)",
            "justify-content": "center",
        }
    switch = dbc.Row(
        justify="center",
        align="start",
        style={"paddingTop": "20px"},
        children=[
            dbc.Col(
                lg=6,
                sm=8,
                children=[
                    dbc.Switch(
                        id="curated-input",
                        label=text,
                        value=False,
                        style=style,
                    ),
                ],
            )
        ],
    )
    return switch


def caching(app):
    @cache.cached(timeout=300, key_prefix="joined-ships")
    def global_store(value):
        query = """
        SELECT j.*, t."order", t.family, f.longFamilyID, f.familyName, a.accession_tag
        FROM joined_ships j
        JOIN taxonomy t ON j.taxid = t.id
        JOIN family_names f ON j.ship_family_id = f.id
        JOIN accessions a ON j.ship_id = a.id
        """
        df = pd.read_sql_query(query, engine)
        if value == "curated":
            df_filtered = df[df["curated_status"] == "curated"]
        else:
            df_filtered = df
        unique_df = df_filtered.drop_duplicates(subset=["accession_tag"])
        return unique_df

    @app.callback(
        [
            Output("curated-status", "data"),
            Output("curated-dataset", "data"),
        ],
        [Input("curated-input", "value")],
    )
    def global_store(switches_value):
        try:
            initial_df = global_store(switches_value)
            df_filtered = initial_df
            curated_status = ""

            if switches_value:
                df_filtered = initial_df[initial_df["curated_status"] == "curated"]
                curated_status = "curated"
            return curated_status, df_filtered.to_dict(orient="records")
        except Exception as e:
            print(f"Error processing curated data: {e}")
            return "", []

    def update_graph_1(data):
        return create_sunburst_plot(df=data, type="ship")

    def update_graph_2(data):
        return create_sunburst_plot(df=data, type="tax")

    def update_ship_table_data(data):
        return data.to_dict("records")

    def update_phylogeny(highlight_families="all"):
        tree = plot_tree(highlight_families)
        return tree

    @app.callback(
        [
            Output("pie1-cache", "data"),
            Output("pie2-cache", "data"),
            Output("phylogeny-cache", "data"),
            Output("explore-table-cache", "data"),
        ],
        Input("curated-dataset", "data"),
    )
    def cache_components(cached_data):
        try:
            initial_df = pd.DataFrame(cached_data)
            unique_df = initial_df.drop_duplicates(subset=["accession_tag"])

            # Define caching functions
            def cache_ship_table():
                return make_ship_table(
                    df=unique_df,
                    id="explore-table",
                    columns=[
                        "accession_tag",
                        "familyName",
                        "order",
                        "family",
                        "genus",
                        "species",
                    ],
                )

            def cache_ship_sunburst():
                return create_sunburst_plot(df=unique_df, type="ship")

            def cache_tax_sunburst():
                return create_sunburst_plot(df=unique_df, type="tax")

            def cache_tree():
                return plot_tree(highlight_families="all")

            pie1_data = cache_ship_sunburst()
            pie2_data = cache_tax_sunburst()
            phylogeny_data = cache_tree()
            explore_table_data = cache_ship_table()

            return pie1_data, pie2_data, phylogeny_data, explore_table_data

        except Exception as e:
            print(f"Error caching components: {e}")
            return [], [], [], []

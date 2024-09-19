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

from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot
from src.utils.tree import plot_tree


download_ships_button = dbc.Button(
    html.Div(
        [
            "Download the latest version of ",
            html.Span(
                "starbase",
                className="logo-text",
            ),
            ".",
        ],
    ),
    id="open-modal",
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
                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the meantime, you can retrieve Starship sequences below:",
                    ],
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
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
        dcc.Download(id="dl-package"),
    ],
    className="text-center",
    style={
        "font-size": "0.875rem",
    },
)

def curated_switch(text,size="normal"):
    if size == "normal":
        style={
                            "display": "flex",
                            "alignItems": "baseline",
                            "justify-content": "center",
                        }
    if size == "large":
        style={
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

modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Choose Starships to Download")),
        dbc.ModalBody(
            html.Div(
                [
                    dmc.Grid(
                        justify="space-within",
                        children=[
                            dmc.GridCol(
                                span="content",
                                children=[
                                    dbc.Button(
                                        "Download All Starships", id="download-all-btn"
                                    )
                                ],
                            ),
                            dmc.GridCol(
                                span="content",
                                children=[
                                    dbc.Button(
                                        "Download Selected Starships",
                                        id="download-selected-btn",
                                    )
                                ],
                            ),
                            html.Div(id="download-table"),
                            dcc.Download(id="download-fasta"),
                        ],
                    )
                ]
            )
        ),
        dbc.ModalFooter(
            dbc.Button(
                "Close",
                id="close",
                className="ms-auto",
                n_clicks=0,
            )
        ),
    ],
    id="download-modal",
    is_open=False,
    size="lg",
)


def caching(app, engine, cache):
    @app.callback(
        Output("joined-ships", "data"),
    )
    @cache.memoize()
    def ship_metadata():
        query = """
        SELECT j.*, t."order", t.family, f.longFamilyID, f.familyName, a.accession_tag
        FROM joined_ships j
        JOIN taxonomy t ON j.taxid = t.id
        JOIN family_names f ON j.ship_family_id = f.id
        JOIN accessions a ON j.ship_id = a.id
        """
        try:
            df = pd.read_sql_query(query, engine)
            return df.to_dict(orient="records")
        except Exception as e:
            print(f"Error fetching ship metadata: {e}")
            return []

    @app.callback(
        [
            Output("curated-status", "data"),
            Output("curated-dataset", "data"),
        ],
        [
            Input("curated-input", "value"),
            Input("joined-ships", "data"),
        ]
    )
    def make_curated(switches_value, cached_data):
        try:
            initial_df = pd.DataFrame(cached_data)
            df_filtered = initial_df
            curated_status = ""

            if switches_value:
                df_filtered = initial_df[initial_df["curated_status"] == "curated"]
                curated_status = "curated"
            return curated_status, df_filtered.to_dict(orient="records")
        except Exception as e:
            print(f"Error processing curated data: {e}")
            return "", []

    @app.callback(
        Output("paper-cache", "data"),
        Input("url","href")
    )
    @cache.memoize()
    def ship_papers(url):
        if url:
            query = """
            SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, p.shortCitation, f.familyName, f.type_element_reference
            FROM papers p
            JOIN family_names f ON p.shortCitation = f.type_element_reference
            """
            try:
                df = pd.read_sql_query(query, engine)
                data = df.to_dict(orient="records")
                return data
            except Exception as e:
                print(f"Error fetching ship papers: {e}")
                return []

    @app.callback(
        [
            Output("pie1-cache", "data"),
            Output("pie2-cache", "data"),
            Output("phylogeny-cache", "data"),
            Output("explore-table-cache", "data"),
        ],
        Input("curated-dataset", "data")
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


# TODO: update to generate based on sql query
def dl_package(app):
    @app.callback(Output("dl-package", "data"), [Input("dl-button", "n_clicks")])
    def generate_download(n_clicks):
        if n_clicks and n_clicks > 0:
            return dcc.send_file(
                "database_folder/Starships/ships/fna/blastdb/concatenated.fa"
            )
        else:
            return dash.no_update

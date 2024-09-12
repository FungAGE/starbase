import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

import dash
from dash import dcc, callback, html
import dash_bootstrap_components as dbc

from dash.dependencies import Output, Input, State
import base64
import pandas as pd
import sqlite3

from src.utils.parsing import parse_fasta, parse_gff


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

download_starbase_button = (
    html.Div(
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
    ),
)


def dl_package(app):
    @app.callback(Output("dl-package", "data"), [Input("dl-button", "n_clicks")])
    def generate_download(n_clicks):
        if n_clicks and n_clicks > 0:
            return dcc.send_file(
                "database_folder/Starships/ships/fna/blastdb/concatenated.fa"
            )
        else:
            return dash.no_update


def update_fasta_upload(app):
    @app.callback(
        Output("fasta-sequence-upload", "children"),
        [
            Input("fasta-upload", "contents"),
            Input("fasta-upload", "filename"),
        ],
    )
    def update_fasta_details(seq_content, seq_filename):
        if seq_content is None:
            return [
                html.Div(
                    html.P(
                        ["Select a FASTA file to upload"],
                    )
                )
            ]
        else:
            try:
                # "," is the delimeter for splitting content_type from content_string
                content_type, content_string = seq_content.split(",")
                query_string = base64.b64decode(content_string).decode("utf-8")
                children = parse_fasta(query_string, seq_filename)
                return children

            except Exception as e:
                logging.error(e)
                return html.Div(["There was an error processing this file."])


def update_gff_upload(app):
    @app.callback(
        Output("output-gff-upload", "children"),
        [
            Input("upload-gff", "contents"),
            Input("upload-gff", "filename"),
        ],
    )
    def update_gff_details(anno_content, anno_filename):
        if anno_content is None:
            return [
                html.Div(["Select a GFF file to upload"]),
            ]
        else:
            try:
                children = parse_gff(anno_content, anno_filename)
                return children

            except Exception as e:
                logging.error(e)
                return html.Div(["There was an error processing this file."])


def load_ship_metadata(app):
    @app.callback(Output("joined-ships", "data"), Input("url", "href"))
    def ship_metadata(url):
        if url:
            try:
                conn = sqlite3.connect("database_folder/starbase.sqlite")
                query = """
                SELECT j.*, o."order", o.family, f.longFamilyID, f.familyName
                FROM joined_ships j
                JOIN genome_taxonomy o ON j.taxid = o.taxID
                JOIN family_names f ON j.ship_family_id = f.id
                """
                df = pd.read_sql_query(query, conn)
                data = df.to_dict(orient="records")

                return data

            except sqlite3.Error as error:
                logging.error("Failed to retrieve data from SQLite table:", error)
                return None

            finally:
                if conn:
                    conn.close()

        return None


def update_dataset(app):

    @app.callback(
        [
            Output("curated-status", "data"),
            Output("curated-dataset", "data"),
        ],
        [
            Input("curated-input", "value"),
            Input("joined-ships", "data"),
        ],
    )
    def curated_switch(switches_value, cached_data):
        initial_df = pd.DataFrame(cached_data)
        df_filtered = initial_df
        curated_status = ""

        if switches_value:
            df_filtered = initial_df[initial_df["curated_status"] == "curated"]
            curated_status = "curated "
        data = df_filtered.to_dict(orient="records")
        return curated_status, data


def load_ship_papers(app):
    @app.callback(
        Output("paper-cache", "data"),
        Input("url", "href"),
    )
    def ship_papers(url):
        if url:
            try:
                conn = sqlite3.connect("database_folder/starbase.sqlite")
                query = """
                SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, p.shortCitation, f.familyName, f.type_element_reference
                FROM papers p
                JOIN family_names f ON p.shortCitation = f.type_element_reference
                """
                df = pd.read_sql_query(query, conn)
                data = df.to_dict(orient="records")

                return data

            except sqlite3.Error as error:
                logging.error("Failed to retrieve data from SQLite table:", error)
                return None

            finally:
                if conn:
                    conn.close()

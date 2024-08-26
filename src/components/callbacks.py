import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, callback, html
from dash.dependencies import Output, Input, State
import base64
import pandas as pd
import sqlite3

from src.utils.parsing import parse_fasta, parse_gff


def dl_package(app):
    @app.callback(Output("dl-package", "data"), [Input("dl-button", "n_clicks")])
    def generate_download(n_clicks):
        if n_clicks is None:
            return dash.no_update
        return dcc.send_file(
            "database_folder/Starships/ships/fna/blastdb/concatenated.fa"
        )


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
                        style={
                            "justify-content": "center",
                            "textAlign": "center",
                        },
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
                print(e)
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
                print(e)
                return html.Div(["There was an error processing this file."])


def load_ship_metadata(app):
    @app.callback(
        Output("joined-ships", "data"),  # Output to dcc.Store
        Input("url", "href"),  # Trigger the callback when the URL changes
    )
    def starship_metadata_table(url):
        if url:
            try:
                conn = sqlite3.connect("database_folder/starbase.sqlite")
                query = """
                SELECT j.*, o."order", o.family, f.longFamilyID, f.familyName, f.type_element_reference
                FROM joined_ships j
                JOIN genome_taxonomy o ON j.taxid = o.taxID
                JOIN family_names f ON j.ship_family_id = f.id
                """
                df = pd.read_sql_query(query, conn)
                data = df.to_dict(orient="records")

                return data

            except sqlite3.Error as error:
                print("Failed to retrieve data from SQLite table:", error)
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
        Input("curated-input", "value"),
        State("joined-ships", "data"),
    )
    def curated_switch(switches_value, cached_data):
        initial_df = pd.DataFrame(cached_data)

        if switches_value:
            df_filtered = initial_df[initial_df["curated_status"] == "curated"]
            curated_status = "curated "
        else:
            df_filtered = initial_df
            curated_status = ""
        data = df_filtered.to_dict(orient="records")

        return (
            curated_status,
            data,
        )

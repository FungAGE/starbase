import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

import dash
from dash import dcc, callback, html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.dependencies import Output, Input, State
import base64
import pandas as pd
import sqlite3
from dash.exceptions import PreventUpdate

from src.utils.parsing import parse_fasta, parse_gff
from src.components.tables import make_ship_table


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


# Callback to handle modal opening and table creation
def modal_download(app):
    @app.callback(
        [
            Output("download-modal", "is_open"),
            Output("download-table", "children"),
        ],
        [
            Input("open-modal", "n_clicks"),
            Input("joined-ships", "data"),
            Input("close", "n_clicks"),
        ],
        State("download-modal", "is_open"),
    )
    def create_modal_table(dl_click, cached_data, close_modal, is_open):
        if dl_click and not close_modal:
            modal = not is_open
        elif close_modal:
            modal = False
        else:
            raise PreventUpdate

        if modal:
            initial_df = pd.DataFrame(cached_data)
            table = make_ship_table(
                df=initial_df,
                id="download-table",
                pg_sz=25,
                columns=[
                    "starshipID",
                    "familyName",
                    "genus",
                    "species",
                ],
            )
            return modal, table

        return modal, None


# Callback to handle FASTA download
def dl_fa(app):
    @app.callback(
        Output("download-fasta", "data"),
        Input("download-all-btn", "n_clicks"),
        Input("download-selected-btn", "n_clicks"),
        [
            State("download-table", "derived_virtual_data"),
            State("download-table", "derived_virtual_selected_rows"),
        ],
        prevent_initial_call=True,
    )
    def download_fasta(click_all, click_selected, rows, selected_rows):
        if click_all:
            selected_rows = list(range(len(rows)))

        if click_selected:
            if selected_rows is None:
                raise ValueError("Select rows in the table first")

        if click_all or click_selected:
            ship_names = [
                rows[selected_row]["starshipID"] for selected_row in selected_rows
            ]
            try:
                conn = sqlite3.connect("database_folder/starbase.sqlite")
                query = f"SELECT genome_name, genome_sequence FROM genome_genome WHERE genome_name "
                if len(ship_names) > 1:
                    placeholders = ",".join(["?"] * len(ship_names))
                    query += f"IN ({placeholders})"
                else:
                    query += "= ?"

                df = pd.read_sql_query(query, conn, params=ship_names)

                # Create FASTA content
                fasta_content = [
                    f">{row['genome_name']}\n{row['genome_sequence']}"
                    for _, row in df.iterrows()
                ]
                fasta_str = "\n".join(fasta_content)

                # Send the FASTA file for download
                return dcc.send_string(fasta_str, filename="starships.fasta")

            except sqlite3.Error as error:
                print("Failed to retrieve data from SQLite table:", error)
                return None

            finally:
                if conn:
                    conn.close()
        else:
            return None

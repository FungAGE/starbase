import warnings

warnings.filterwarnings("ignore")

import dash
from dash import dcc, callback, html
from dash.dependencies import Output, Input
import base64
import pandas as pd
from Bio import SeqIO
import io
import dash_bootstrap_components as dbc
import sqlite3


def parse_gff(contents, filename):
    content_type, content_string = contents.split(",")
    decoded = base64.b64decode(content_string)
    gff = pd.read_csv(io.StringIO(decoded.decode("utf-8")), sep="\t")
    # cols = [
    #     "seqid",
    #     "source",
    #     "type",
    #     "start",
    #     "end",
    #     "score",
    #     "strand",
    #     "phase",
    #     "attributes",
    # ]

    # annotations = pd.DataFrame(gff)
    nanno = len(gff)

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of annotations: {nanno}"),
            ]
        )
    ]


def parse_fasta(contents, filename):
    sequences = SeqIO.parse(io.StringIO(contents), "fasta")

    records = []
    nseq = 0
    for sequence in sequences:
        records.append({"ID": sequence.id, "Sequence": str(sequence.seq)})
        nseq += 1

    return [
        html.Div(
            [
                html.H6(f"File name: {filename}"),
                html.H6(f"Number of sequences: {nseq}"),
            ],
        )
    ]


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
            Output("curated-dataset", "data"),
            Output("ship-card-header", "children"),
            Output("ship-count", "children"),
            Output("species-card-header", "children"),
            Output("species-count", "children"),
            Output("curated-status", "data"),
        ],
        [Input("joined-ships", "data"), Input("curated-input", "value")],
    )
    def curated_switch(cached_data, switches_value):
        initial_df = pd.DataFrame(cached_data)

        if switches_value:
            df_filtered = initial_df[initial_df["curated_status"] == "curated"]
            curated_status = "curated "
        else:
            df_filtered = initial_df
            curated_status = ""

        ship_count = df_filtered["starshipID"].nunique()
        species = df_filtered["genus"] + "-" + df_filtered["species"]
        species_count = species.nunique()

        data = df_filtered.to_dict(orient="records")

        ship_card_header = html.H4(
            html.P(
                [
                    f"Total number {curated_status}of Starships in ",
                    html.Span(
                        "starbase",
                        className="logo-text",
                    ),
                    ":",
                ]
            ),
            className="card-title",
        )
        species_card_header = html.H4(
            f"Fungal species with {curated_status}Starships",
            className="card-title",
        )

        # Return the filtered data, card headers, and counts
        return (
            data,
            ship_card_header,
            ship_count,
            species_card_header,
            species_count,
            curated_status,
        )

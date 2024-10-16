import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)

import dash
from dash import html, dcc, dash_table, callback
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash.dependencies import Output, Input
import pandas as pd

from src.components.sqlite import engine

dash.register_page(__name__)

table_columns = [
    {
        "name": "Accession",
        "id": "accession_tag",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Starship Family",
        "id": "familyName",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Order",
        "id": "order",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Family",
        "id": "family",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Species",
        "id": "species",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
]


layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dmc.Grid(
            justify="center",
            align="center",
            style={"paddingTop": "20px"},
            gutter="xl",
            children=[
                dmc.GridCol(
                    span=12,
                    children=[dmc.Center([dmc.Title("Choose Starships to Download:")])],
                ),
                dmc.GridCol(
                    span=12,
                    children=[
                        dmc.Center(
                            [
                                dbc.Button(
                                    "Download All Starships",
                                    id="download-all-btn",
                                    class_name="text-custom text-custom-sm text-custom-md",
                                ),
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    span=12,
                    children=[
                        dmc.Center(
                            [
                                dbc.Button(
                                    "Download Selected Starships",
                                    id="download-selected-btn",
                                    class_name="text-custom text-custom-sm text-custom-md",
                                ),
                            ]
                        )
                    ],
                ),
                dcc.Download(id="dl-package"),
                html.Div(id="dl-warning"),
                dmc.GridCol(
                    span=10,
                    children=[
                        dmc.Center(
                            [
                                html.Div(
                                    [
                                        dash_table.DataTable(
                                            id="dl-table",
                                            data=[],
                                            columns=table_columns,
                                            selected_columns=[],
                                            selected_rows=[],
                                            editable=False,
                                            filter_action="native",
                                            sort_action="native",
                                            sort_mode="multi",
                                            row_selectable="multi",
                                            row_deletable=False,
                                            page_action="native",
                                            page_current=0,
                                            page_size=24,
                                            markdown_options={"html": True},
                                            style_table={
                                                "overflowX": "auto",
                                            },
                                            style_data={
                                                "whiteSpace": "minimal",
                                            },
                                            style_cell={
                                                "minWidth": "100px",
                                                "maxWidth": "100%",
                                                "textAlign": "left",
                                            },
                                        ),
                                    ],
                                    className="auto-resize-900",
                                )
                            ]
                        )
                    ],
                ),
            ],
        ),
    ],
)


@callback(Output("dl-table", "data"), Input("url", "href"))
def make_dl_table(url):
    try:
        query = """
        SELECT a.accession_tag, f.familyName, t."order", t.family, t.species 
        FROM accessions a
        LEFT JOIN joined_ships j ON a.id = j.ship_id
        LEFT JOIN taxonomy t ON j.taxid = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        WHERE j.orphan IS NULL
        """
        df = pd.read_sql_query(query, engine)
        if df.empty:
            logging.error("Query returned an empty DataFrame.")
        else:
            logging.info(f"Retrieved {len(df)} records from the database.")

        df.fillna("", inplace=True)  # Replace NaN with an empty string
        return df.to_dict("records")

    except Exception as e:
        logging.error(f"Failed to execute query in make_dl_table. Details: {e}")
        return []


@callback(
    [Output("dl-package", "data"), Output("dl-warning", "children")],
    [
        Input("download-all-btn", "n_clicks"),
        Input("download-selected-btn", "n_clicks"),
        Input("dl-table", "derived_virtual_selected_rows"),
        Input("dl-table", "derived_virtual_selected_data"),
    ],
)
def generate_download(dl_all, dl_select, rows, data):
    try:
        if dl_all or dl_select:
            logging.info(f"dl_all={dl_all}, dl_select={dl_select}, rows={rows}")

            if dl_select and (rows is None or rows == []):
                logging.warning(
                    "Download selected was triggered but no rows are selected."
                )
                return None, dbc.Alert(
                    "Make a selection in the table first", color="warning"
                )

            # Validate if data is available
            if not data:
                logging.error("No data available from the table for processing.")
                return None, dbc.Alert("No data available for download", color="danger")

            # Extract accessions safely
            try:
                accessions = pd.DataFrame(data).iloc[rows]["accession_tag"]
                logging.info(f"Selected accessions: {accessions.tolist()}")
            except Exception as e:
                logging.error(f"Failed to extract accessions. Details: {e}")
                return None, dbc.Alert(
                    "Error extracting accessions from the selection", color="danger"
                )

            # Build and execute the SQL query
            try:
                query = """
                SELECT s.*
                FROM ships s
                LEFT JOIN accessions a ON s.accession = a.id
                """
                if accessions is not None and len(accessions) > 0:
                    placeholders = ",".join(["?"] * len(accessions))
                    query += f" WHERE accession_tag IN ({placeholders})"
                    df = pd.read_sql_query(query, engine, params=accessions.tolist())
                    logging.info(f"Query returned {len(df)} records for accessions.")
                else:
                    df = pd.read_sql_query(query, engine)
                    logging.info("Retrieved all records from the database.")
            except Exception as e:
                logging.error(f"Failed to execute database query. Details: {e}")
                return None, dbc.Alert(
                    "Error retrieving data from the database", color="danger"
                )

            if df.empty:
                logging.warning("Query returned no matching records.")
                return None, dbc.Alert(
                    "No matching records found for the selected accessions",
                    color="warning",
                )

            # Create FASTA content
            try:
                fasta_content = [
                    f">{row['accession_tag']}\n{row['ship_sequence']}"
                    for _, row in df.iterrows()
                ]
                fasta_str = "\n".join(fasta_content)
                logging.info("FASTA content created successfully.")
                return dcc.send_string(fasta_str, filename="starships.fasta"), None
            except Exception as e:
                logging.error(f"Failed to create FASTA content. Details: {e}")
                return None, dbc.Alert(
                    "Error creating FASTA file for download", color="danger"
                )
        else:
            logging.info("No download action triggered.")
            return dash.no_update, None

    except Exception as e:
        logging.error(
            f"An unexpected error occurred in generate_download. Details: {e}"
        )
        return None, dbc.Alert(
            "An unexpected error occurred during the download process", color="danger"
        )

import dash
from dash import html, dcc, dash_table, callback
import dash_mantine_components as dmc
from dash.dependencies import Output, Input
import pandas as pd

from src.components.cache import cache
from src.components.sql_engine import starbase_engine
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import fetch_download_data, fetch_all_ships

import logging

logger = logging.getLogger(__name__)

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
        dmc.Center(
            children=[
                dmc.Title(
                    ["Choose Starships to Download:"], style={"paddingTop": "20px"}
                ),
            ],
        ),
        dmc.Grid(
            justify="center",
            align="center",
            style={"paddingTop": "20px"},
            gutter="xl",
            children=[
                dmc.GridCol(
                    span=12,
                    children=[
                        dmc.Center(
                            children=[
                                dmc.Stack(
                                    children=[
                                        html.Div(
                                            [
                                                dmc.NotificationProvider(),
                                                html.Div(id="dl-notify"),
                                                dmc.SimpleGrid(
                                                    cols={"md": 2, "sm": 1},
                                                    spacing={"base": 10, "sm": "xl"},
                                                    verticalSpacing={
                                                        "base": "md",
                                                        "sm": "xl",
                                                    },
                                                    children=[
                                                        dmc.Button(
                                                            "Download All Starships",
                                                            id="download-all-btn",
                                                            className="text-custom text-custom-sm text-custom-md",
                                                            style={"width": "100%"},
                                                        ),
                                                        dmc.Button(
                                                            "Download Selected Starships",
                                                            id="download-selected-btn",
                                                            className="text-custom text-custom-sm text-custom-md",
                                                            style={"width": "100%"},
                                                        ),
                                                    ],
                                                ),
                                            ]
                                        ),
                                        dcc.Download(id="dl-package"),
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
                                                        "overflowY": "auto",
                                                        "maxHeight": "400px",
                                                    },
                                                    style_data={
                                                        "height": "auto",
                                                        "whiteSpace": "minimal",
                                                        "overflow": "hidden",
                                                    },
                                                    style_cell={
                                                        "minWidth": "0px",
                                                        "maxWidth": "300px",
                                                        "whiteSpace": "minimal",
                                                        "textAlign": "left",
                                                    },
                                                    style_header={
                                                        "backgroundColor": "lightgrey",
                                                        "fontWeight": "bold",
                                                        "textAlign": "left",
                                                    },
                                                ),
                                            ],
                                        ),
                                    ]
                                )
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@cache.memoize()
@callback(Output("dl-table", "data"), Input("url", "href"))
def make_dl_table(url):
    try:
        df = load_from_cache("download_data")
        if df is None:
            df = fetch_download_data()
        logger.info(f"Retrieved {len(df)} records from the database.")

        df.fillna("", inplace=True)  # Replace NaN with an empty string
        return df.to_dict("records")

    except Exception as e:
        logger.error(f"Failed to execute query in make_dl_table. Details: {e}")
        return []


def notification_base(title, message):
    dmc.Notification(
        title=title,
        id="simple-notify",
        action="show",
        message=message,
    )


@callback(
    [
        Output("dl-package", "data"),
        Output("dl-notify", "children"),
        Output("download-all-btn", "disabled"),
        Output("download-selected-btn", "disabled"),
    ],
    [
        Input("download-all-btn", "n_clicks"),
        Input("download-selected-btn", "n_clicks"),
        Input("dl-table", "data"),
        Input("dl-table", "derived_virtual_selected_rows"),
    ],
)
def generate_download(dl_all, dl_select, table_data, selected_rows):
    logger.info(
        f"dl_all={dl_all}, dl_select={dl_select}, table_data_length={len(table_data)}, selected_rows={selected_rows}"
    )

    if not table_data or len(table_data) == 0:
        logger.error("No data available from the table for processing.")
        return (
            dash.no_update,
            notification_base(title="error", message="No data available for download"),
            False,
            False,
        )

    if not dl_all and not dl_select:
        logger.debug("No download action triggered. dl_all or dl_select not active.")
        return dash.no_update, None, False, False

    table_df = pd.DataFrame(table_data)
    try:

        if dl_all:
            accessions = table_df["accession_tag"].to_list()
            logger.info("Using all table data for download.")

            df = load_from_cache("all_ships")
            if df is None:
                df = fetch_all_ships()

        elif dl_select:
            if not selected_rows or len(selected_rows) == 0:
                logger.warning(
                    "Download selected was triggered but no rows are selected."
                )
                return (
                    dash.no_update,
                    notification_base(
                        title="Warning:", message="Make a selection in the table first"
                    ),
                    False,
                    False,
                )

            else:
                if isinstance(selected_rows, list) and all(
                    isinstance(i, int) for i in selected_rows
                ):
                    selected_df = table_df.iloc[selected_rows]
                else:
                    logger.error(
                        "Invalid index type for iloc. Must be a list of integers."
                    )

                accessions = selected_df["accession_tag"].to_list()
                logger.info(f"Using selected table data: {accessions}")

                df = df[df["accession_tag"].isin(accessions)]

            if df.empty:
                logger.error("No data available for selected rows.")
                return (
                    dash.no_update,
                    notification_base(
                        title="Error:", message="No data available for selected rows."
                    ),
                    False,
                    False,
                )

        if df.empty:
            logger.warning("Query returned no matching records.")
            return (
                dash.no_update,
                notification_base(
                    title="Error:",
                    message="No matching records found for the selected accessions",
                ),
                False,
                False,
            )
        else:
            logger.info(f"Retrieved {len(df)} records from the database.")
            try:
                fasta_content = [
                    f">{row['accession_tag']}\n{row['sequence']}"
                    for _, row in df.iterrows()
                ]
                fasta_str = "\n".join(fasta_content)
                logger.info("FASTA content created successfully.")
                return (
                    dcc.send_string(fasta_str, filename="starships.fasta"),
                    None,
                    True,
                    False,
                )
            except Exception as e:
                logger.error(f"Failed to create FASTA content. Details: {e}")
                return (
                    dash.no_update,
                    notification_base(
                        title="Error:",
                        message="Error when creating FASTA file for download",
                    ),
                    False,
                    False,
                )
    except Exception as e:
        logger.error(f"Failed to execute database query. Details: {e}")
        return (
            dash.no_update,
            notification_base(
                title="Error:",
                message=f"Failed to execute database query. Details: {e}",
            ),
            False,
            False,
        )

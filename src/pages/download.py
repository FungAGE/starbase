import dash
from dash import html, dcc, dash_table, callback
import dash_mantine_components as dmc
from dash.dependencies import Output, Input
from dash_iconify import DashIconify
import pandas as pd

from src.components.callbacks import curated_switch
from src.components.cache import cache
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
    size="xl", 
    children=[
        dcc.Location(id="url", refresh=False),
        
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title(
                    "Download Starship Sequences",
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    "Select individual Starships or download the complete dataset",
                    size="lg",
                    c="dimmed",
                ),
                # Add curated switch here
                dmc.Stack([
                    dmc.Group(
                        children=[
                            curated_switch(text="Only show curated Starships", size="md"),
                        ],
                        mt="md",
                    ),
                    # Download Options
                    dmc.Center(
                        dmc.Group(
                            gap="xl",
                            children=[
                                dmc.Button(
                                    "Download All Starships",
                                    id="download-all-btn",
                                    variant="gradient",
                                    gradient={"from": "indigo", "to": "cyan"},
                                    leftSection=html.I(className="bi bi-cloud-download"),
                                    size="md",
                                ),
                                dmc.Button(
                                    "Download Selected Starships",
                                    id="download-selected-btn",
                                    variant="gradient",
                                    gradient={"from": "teal", "to": "lime"},
                                    leftSection=html.I(className="bi bi-download"),
                                    size="md",
                                ),
                            ],
                            ),
                        ),
                ], gap="xl"),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
        ),
        
        # Main Content
        # Notification Area
        html.Div(id="dl-notify"),
        dcc.Download(id="dl-package"),
        
        # Table Section
        html.Div(
            children=[
                # Table Header with Stats
                dmc.Group(
                    gap="apart",
                    mt="md",
                    children=[
                        dmc.Text(
                            id="table-stats",
                            size="sm",
                            c="dimmed",
                        ),
                        dmc.Text(
                            "Click rows to select Starships",
                            size="sm",
                            c="dimmed",
                            style={"fontStyle": "italic"},
                        ),
                    ],
                ),

                dash_table.DataTable(
                    id="dl-table",
                    columns=table_columns,
                    data=[],
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    row_selectable="multi",
                    page_action="native",
                    page_current=0,
                    page_size=20,
                    markdown_options={"html": True},
                    style_table={
                        "overflowX": "auto",
                        "overflowY": "auto",
                        "maxHeight": "60vh",  # Responsive height
                    },
                    style_data={
                        "height": "auto",
                        "lineHeight": "20px",
                        "padding": "10px",
                    },
                    style_cell={
                        "fontFamily": "Arial, sans-serif",
                        "textAlign": "left",
                        "minWidth": "100px",
                        "maxWidth": "300px",
                        "overflow": "hidden",
                        "textOverflow": "ellipsis",
                    },
                    style_header={
                        "backgroundColor": "#f8f9fa",
                        "fontWeight": "bold",
                        "borderBottom": "2px solid #dee2e6",
                        "textAlign": "left",
                        "padding": "12px",
                    },
                    style_filter={
                        "backgroundColor": "#f8f9fa",
                        "padding": "8px",
                    },
                    style_data_conditional=[
                        {
                            "if": {"row_index": "odd"},
                            "backgroundColor": "#f8f9fa",
                        },
                        {
                            "if": {"state": "selected"},
                            "backgroundColor": "#e3f2fd",
                            "border": "1px solid #2196f3",
                        },
                    ],
                ),
            ],
        ),
    ],
    py="xl",
)


@cache.memoize()
@callback(
    Output("dl-table", "data"), 
    [Input("url", "href"),
     Input("curated-input", "checked")]
)
def make_dl_table(url, curated=True):
    try:
        df = fetch_download_data(curated=curated)
        if df is None:
            return []
            
        logger.info(f"Retrieved {len(df)} records from the database (curated={curated}).")
        df.fillna("", inplace=True)
        return df.to_dict("records")

    except Exception as e:
        logger.error(f"Failed to execute query in make_dl_table. Details: {e}")
        return []


@callback(
    Output("download-selected-btn", "disabled"),
    [Input("dl-table", "derived_virtual_selected_rows")],
)
def update_download_selected_button(selected_rows):
    # Button is disabled if no rows are selected
    return not selected_rows or len(selected_rows) == 0


@callback(
    [
        Output("dl-package", "data"),
        Output("notifications-container", "children"),
        Output("download-all-btn", "disabled"),
    ],  # Removed download-selected-btn from outputs
    [
        Input("download-all-btn", "n_clicks"),
        Input("download-selected-btn", "n_clicks"),
        Input("dl-table", "data"),
        Input("dl-table", "derived_virtual_selected_rows"),
    ],
    prevent_initial_call=True,
)
def generate_download(dl_all_clicks, dl_select_clicks, table_data, selected_rows):
    # Determine which button was clicked using callback context
    ctx = dash.callback_context
    if not ctx.triggered or not any([dl_all_clicks, dl_select_clicks]):
        raise dash.exceptions.PreventUpdate
    
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    if not table_data or len(table_data) == 0:
        return (
            dash.no_update,
            dmc.Notification(
                title="Error",
                message="No data available for download",
                color="red",
                icon=DashIconify(icon="ic:round-error"),
                action="show",
            ),
            False,  # Only return one boolean for download-all-btn
        )

    try:
        table_df = pd.DataFrame(table_data)
        
        # Get all ships data first
        df = load_from_cache("all_ships")
        if df is None:
            df = fetch_all_ships()
            
        if df is None or df.empty:
            raise ValueError("Failed to fetch ship data from database")

        # Handle download based on which button was clicked
        if button_id == "download-all-btn":
            accessions = table_df["accession_tag"].to_list()
            logger.info("Using all table data for download.")
            
            # Filter all ships by the accessions in the table
            df = df[df["accession_tag"].isin(accessions)]

        elif button_id == "download-selected-btn":
            if not selected_rows or len(selected_rows) == 0:
                logger.warning("Download selected was triggered but no rows are selected.")
                return (
                    dash.no_update,
                    dmc.Notification(
                        title="Warning",
                        message="Make a selection in the table first",
                        color="yellow",
                        icon=DashIconify(icon="ic:round-warning"),
                        action="show",
                    ),
                    False,
                    False,
                )

            # Get selected accessions
            selected_df = table_df.iloc[selected_rows]
            accessions = selected_df["accession_tag"].to_list()
            logger.info(f"Using selected table data: {accessions}")
            
            # Filter all ships by selected accessions
            df = df[df["accession_tag"].isin(accessions)]
        else:
            return dash.no_update, None, False, False

        if df.empty:
            logger.warning("No matching records found.")
            return (
                dash.no_update,
                dmc.Notification(
                    title="Error",
                    message="No matching records found for the selected accessions",
                    color="red",
                    icon=DashIconify(icon="ic:round-error"),
                    action="show",
                ),
                False,
                False,
            )

        # Create FASTA file
        try:
            fasta_content = [
                f">{row['accession_tag']}\n{row['sequence']}"
                for _, row in df.iterrows()
            ]
            fasta_str = "\n".join(fasta_content)
            logger.info(f"FASTA content created successfully for {len(df)} sequences.")
            
            return (
                dcc.send_string(fasta_str, filename="starships.fasta"),
                dmc.Notification(
                    title="Success",
                    message=f"Downloaded {len(df)} Starship sequences",
                    color="green",
                    icon=DashIconify(icon="ic:round-check-circle"),
                    action="show",
                ),
                True,
                False,
            )
            
        except Exception as e:
            logger.error(f"Failed to create FASTA content. Details: {e}")
            return (
                dash.no_update,
                dmc.Notification(
                    title="Error",
                    message="Error when creating FASTA file for download",
                    color="red",
                    icon=DashIconify(icon="ic:round-error"),
                    action="show",
                ),
                False,
                False,
            )
            
    except Exception as e:
        logger.error(f"Failed to execute database query. Details: {e}")
        return (
            dash.no_update,
            dmc.Notification(
                title="Error",
                message="Failed to download sequences",
                color="red",
                icon=DashIconify(icon="ic:round-error"),
                action="show",
            ),
            False,
            False,
        )

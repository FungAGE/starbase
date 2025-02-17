import dash
from dash import html, dcc, dash_table, callback
import dash_mantine_components as dmc
from dash.dependencies import Output, Input, State
from dash_iconify import DashIconify
import pandas as pd

from src.components.callbacks import curated_switch, create_accession_modal, create_modal_callback, dereplicated_switch
from src.config.cache import cache

from src.database.sql_manager import fetch_download_data, fetch_all_ships
from src.components.tables import make_dl_table
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
        "cellStyle": {"cursor": "pointer", "color": "#1976d2"}
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
                            dereplicated_switch(text="Only show dereplicated Starships", size="md"),
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
                dmc.Modal(
                    id="accession-modal",
                    opened=False,
                    centered=True,
                    overlayProps={"blur": 3},
                    size="lg",
                    children=[
                        dmc.Title(id="modal-title", order=3),
                        dmc.Space(h="md"),
                        html.Div(id="modal-content"),
                    ],
                ),
                html.Div(
                    id="dl-table-container",
                    children=[make_dl_table(
                        df=pd.DataFrame(),  # Empty initial DataFrame
                        id="dl-table",
                        table_columns=table_columns
                    )]
                ),
            ],
        ),
    ],
    py="xl",
)


@callback(
    [Output("dl-table", "rowData"),
     Output("table-stats", "children")],
    [Input("url", "href"),
     Input("curated-input", "checked"),
     Input("dereplicated-input", "checked")]
)
def update_dl_table(url, curated=True, dereplicate=False):
    logger.debug(f"update_dl_table called with curated={curated}, dereplicate={dereplicate}")
    try:
        df = fetch_download_data(curated=curated, dereplicate=dereplicate)
        if df is None or df.empty:
            logger.warning("fetch_download_data returned None or empty DataFrame")
            return [], "No records found"
            
        logger.info(f"Retrieved {len(df)} records (curated={curated}, dereplicated={dereplicate}).")
        df = df.fillna("")  # Explicitly fill NA values
        
        # Convert to records and ensure all values are strings
        records = [{k: str(v) if pd.notnull(v) else "" for k, v in record.items()} 
                  for record in df.to_dict("records")]
        
        return records, f"Showing {len(records)} records"
        
    except Exception as e:
        logger.error(f"Failed to execute query in make_dl_table. Details: {e}")
        return [], "Error loading data"

@callback(
    [Output("dl-package", "data"),
     Output("notifications-container", "children")],
    [Input("download-all-btn", "n_clicks"),
     Input("download-selected-btn", "n_clicks"),
     Input("dl-table", "rowData"),
     Input("dl-table", "selectedRows")],
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
                style={
                    "width": "100%",
                    "maxWidth": "500px",  # Limit width on larger screens
                    "margin": "0 auto",   # Center the notification
                    "@media (max-width: 600px)": {
                        "maxWidth": "100%"  # Full width on mobile
                    }
                }
            )
        )

    try:
        table_df = pd.DataFrame(table_data)
        
        # Get all ships data first
        df = cache.get("all_ships")
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
                    )
                )

            # Get selected accessions directly from selected_rows
            accessions = [row["accession_tag"] for row in selected_rows]
            logger.info(f"Using selected table data: {accessions}")
            
            # Filter all ships by selected accessions
            df = df[df["accession_tag"].isin(accessions)]
        else:
            return dash.no_update, None

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
                )
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
                )
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
                )
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
            )
        )

@callback(
    Output("download-selected-btn", "disabled"),
    [Input("dl-table", "selectedRows")]
)
def update_download_selected_button(selected_rows):
    # Disable the button if no rows are selected
    return not selected_rows or len(selected_rows) == 0

@callback(
    [Output("accession-modal", "opened"),
     Output("modal-content", "children"),
     Output("modal-title", "children")],
    [Input("dl-table", "cellClicked")],
    [State("dl-table", "rowData")],
    prevent_initial_call=True
)
def toggle_modal(cell_clicked, row_data):
    if not cell_clicked or cell_clicked["colId"] != "accession_tag":
        raise dash.exceptions.PreventUpdate
        
    accession = cell_clicked["value"]
    return create_accession_modal(accession)
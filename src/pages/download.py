import dash
from dash import html, dcc, dash_table, callback
import dash_mantine_components as dmc
from dash.dependencies import Output, Input, State
from dash_iconify import DashIconify
import pandas as pd
from datetime import datetime

from src.components.callbacks import curated_switch, create_accession_modal, create_modal_callback, dereplicated_switch
from src.config.cache import cache

from src.database.sql_manager import fetch_download_data, fetch_ships
from src.components.tables import make_dl_table, table_loading_alert, table_no_results_alert, table_error
import logging
from src.utils.seq_utils import create_ncbi_style_header

logger = logging.getLogger(__name__)

dash.register_page(__name__)

table_columns = [
    {
        "name": "Accession",
        "id": "accession_tag",
        "deletable": False,
        "selectable": True,
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
        "name": "Reference",
        "id": "shortCitation",
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
        "id": "name",
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
                                    loaderProps={"variant": "dots", "color": "white"},
                                ),
                                dmc.Button(
                                    "Download Selected Starships",
                                    id="download-selected-btn",
                                    variant="gradient",
                                    gradient={"from": "teal", "to": "lime"},
                                    leftSection=html.I(className="bi bi-download"),
                                    size="md",
                                    loaderProps={"variant": "dots", "color": "white"},
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
                    children=[make_dl_table(pd.DataFrame(), "dl-table", table_columns)]
                ),
            ],
        ),
    ],
    py="xl",
)


@callback(
    [Output("dl-table-container", "children"),
     Output("table-stats", "children")],
    [Input("curated-input", "checked"),
     Input("dereplicated-input", "checked")]
)
def update_dl_table(curated, dereplicate):
    logger.debug(f"update_dl_table called with curated={curated}, dereplicate={dereplicate}")
    
    try:
        df = fetch_download_data(curated=curated, dereplicate=dereplicate)
        if df is None or df.empty:
            logger.warning("fetch_download_data returned None or empty DataFrame")
            return table_no_results_alert(), "No records found"
            
        logger.info(f"Retrieved {len(df)} records (curated={curated}, dereplicated={dereplicate}).")
        df = df.fillna("")  # Explicitly fill NA values
        
        # Convert to records and ensure all values are strings
        records = [{k: str(v) if pd.notnull(v) else "" for k, v in record.items()} 
                  for record in df.to_dict("records")]
        
        if not records:  # Add explicit check for empty records
            return table_no_results_alert(), "No records found"
            
        # Create the table with the data
        table = make_dl_table(
            df=records,  # Pass the records directly
            id="dl-table",
            table_columns=table_columns
        )
        
        return table, f"Showing {len(records)} records"
        
    except Exception as e:
        logger.error(f"Failed to execute query to generate table. Details: {e}")
        return table_error(e), f"Error loading data: {str(e)}"

def generate_download_helper(rows, curated, dereplicate):
    """Helper function containing the common download logic"""
    try:
        if not rows:
            raise ValueError("No rows selected for download")
            
        accessions = [row["accession_tag"] for row in rows]
        dl_df = fetch_ships(accession_tags=accessions, curated=curated, dereplicate=dereplicate)
            
        if dl_df is None or dl_df.empty:
            raise ValueError("No sequences found for download")
            
        # Count occurrences of each accession tag
        accession_counts = dl_df['accession_tag'].value_counts()
        
        fasta_content = []
        for _, row in dl_df.drop_duplicates(subset=['accession_tag', 'sequence']).iterrows():
            count = accession_counts[row['accession_tag']]
            
            if count > 1:
                # Simplified header for multiple representatives
                header = (
                    f">{row['accession_tag']} "
                    f"[family={row['familyName']}] "
                    f"[representatives={count}]"
                )
            else:
                # Full header for single entries
                header = create_ncbi_style_header(row)
            if header is None or header == "None":
                header = (
                    f">{row['accession_tag']}"
                    + (f" [family={row['familyName']}]" if row.get('familyName') else "")
                )
            fasta_content.append(f"{header}\n{row['sequence']}")
            
        fasta_str = "\n".join(fasta_content)
        logger.info(f"FASTA content created successfully for {len(fasta_content)} sequences.")
        
        # Return both the FASTA content and download info
        return dict(
            content=fasta_str,
            filename=f"starships_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
            type="text/plain"
        ), len(fasta_content)
            
    except ValueError as e:
        logger.warning(f"Download helper warning: {str(e)}")
        return None, 0
    except Exception as e:
        logger.error(f"Failed to execute database query. Details: {e}")
        return None, 0

@callback(
    [Output("dl-package", "data", allow_duplicate=True),
     Output("dl-notify", "children", allow_duplicate=True)],
    [Input("download-all-btn", "n_clicks")],
    [State("dl-table", "rowData"),
     State("curated-input", "checked"),
     State("dereplicated-input", "checked")],
    prevent_initial_call=True,
    running=[
        (Output("download-all-btn", "loading"), True, False),
    ],
)
def generate_download_all(dl_all_clicks, table_data, curated, dereplicate):
    if not dl_all_clicks:
        raise dash.exceptions.PreventUpdate
    
    download_data, num_sequences = generate_download_helper(table_data, curated, dereplicate)
    
    if not download_data:
        return None, dmc.Alert(
            children="No sequences found for download",
            color="red",
            title="Error"
        )
        
    return download_data, dmc.Alert(
        children=f"Downloading {num_sequences} sequences",
        color="green",
        title="Success"
    )

@callback(
    [Output("dl-package", "data"),
     Output("dl-notify", "children")],
    [Input("download-selected-btn", "n_clicks")],
    [State("dl-table", "rowData"),
     State("dl-table", "selectedRows"),
     State("curated-input", "checked"),
     State("dereplicated-input", "checked")],
    prevent_initial_call=True,
    running=[
        (Output("download-selected-btn", "loading"), True, False),
    ],
)
def generate_download_selected(dl_select_clicks, table_data, selected_rows, curated, dereplicate):
    if not dl_select_clicks or not selected_rows:
        raise dash.exceptions.PreventUpdate
    
    # Filter table_data to get only the selected rows
    selected_data = [row for row in table_data if row.get("accession_tag") in 
                    [selected.get("accession_tag") for selected in selected_rows]]
    
    download_data, num_sequences = generate_download_helper(selected_data, curated, dereplicate)
    
    if not download_data:
        return None, dmc.Alert(
            children="No sequences found for selected rows",
            color="red",
            title="Error"
        )
        
    return download_data, dmc.Alert(
        children=f"Downloading {num_sequences} sequences",
        color="green",
        title="Success"
    )

@callback(
    Output("download-selected-btn", "disabled"),
    [Input("dl-table", "selectedRows")]  # Only selectedRows, no rowData
)
def update_download_selected_button(selected_rows):
    # Disable the button if no rows are selected
    return not selected_rows or len(selected_rows) == 0

toggle_modal = create_modal_callback(
    "dl-table",
    "accession-modal",
    "modal-content",
    "modal-title"
)
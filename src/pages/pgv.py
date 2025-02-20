import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update
from dash.dependencies import Output, Input, State

import os
import tempfile

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Gff
from matplotlib.lines import Line2D
from pygenomeviz.utils import ColorCycler
from pygenomeviz.align import Blast, AlignCoord

from src.components.callbacks import create_modal_callback
from src.components.error_boundary import handle_callback_error

from pathlib import Path

from src.config.logging import get_logger

logger = get_logger(__name__)

dash.register_page(__name__)

table_columns = [
    {
        "id": "accession_tag",
        "name": "Accession",
        "selectable": True,
    },
    {
        "id": "familyName",
        "name": "Starship Family",
        "selectable": True,
    },
    {
        "id": "name",
        "name": "Species",
        "selectable": True,
    },
]

# Update the layout to match clinker's requirements
layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="cluster-data"),
        
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title("Starship Genome Viewer", order=1, mb="md"),
                dmc.Text(
                    "Compare and visualize up to 4 Starship sequences",
                    size="lg",
                    c="dimmed",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        # Main content
        dmc.Grid(
            children=[
                dmc.GridCol(
                    span={"base": 12, "md": 6},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                dmc.Group(
                                    pos="apart",
                                    children=[
                                        dmc.Title("Select Starships", order=2),
                                        dmc.Button(
                                            dmc.Text("Show Selected Starship(s)", size="lg"),
                                            id="update-button",
                                            variant="gradient",
                                            gradient={"from": "indigo", "to": "cyan"},
                                            leftSection=html.I(className="bi bi-eye"),
                                        ),
                                    ],
                                ),
                                # Clinker Settings
                                dmc.Group(
                                    children=[
                                        dmc.NumberInput(
                                            label="Identity Threshold",
                                            id="identity-threshold",
                                            value=0.3,
                                            min=0,
                                            max=1,
                                            step=0.1,
                                        ),
                                        dmc.Switch(
                                            id="use-file-order",
                                            label="Use File Order",
                                            checked=False,
                                        ),
                                    ],
                                    gap="md",
                                ),
                                # Table
                                dcc.Loading(
                                    id="loading",
                                    type="circle",
                                    children=html.Div(id="pgv-table"),
                                ),
                            ], gap="md"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            h="100%",
                        ),
                    ],
                ),
                # Right Column - Visualization Section
                dmc.GridCol(
                    span={"base": 12, "md": 6},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                html.Div(id="pgv-message", style={"textAlign": "center"}),
                                dcc.Loading(
                                    id="loading-1",
                                    type="circle",
                                    children=html.Div(
                                        id="pgv-figure",
                                        style={
                                            "height": "800px",
                                            "width": "100%",
                                            "overflow": "auto",
                                            "backgroundColor": "#f8f9fa",
                                            "border": "1px solid #dee2e6",
                                            "borderRadius": "4px",
                                        },
                                    ),
                                ),
                                ],
                                gap="md",
                            ),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            h="100%",
                        ),
                    ],
                ),
            ],
            gutter="xl",
        ),
    ],
    py="xl",
)

@callback(
    [Output("pgv-table", "children"), Output("pgv-table-loading", "visible")],
    Input("url", "href"),
    prevent_initial_call=False,
)
@handle_callback_error
def load_ship_table(href):
    """Load and display the ship selection table, filtered to only show entries with GenBank files"""
    from src.database.sql_manager import fetch_ship_table
    from src.components.tables import make_pgv_table
    
    # Get the path to the GenBank files directory
    gbk_dir = Path("/home/adrian/Systematics/Starship_Database/starbase/src/database/db/ships/gbks")
    
    # Get list of available GenBank files and extract starshipIDs
    available_starships = {path.stem.split('_', 2)[2] for path in gbk_dir.glob("*.gbk")}
    
    # Fetch all ships and filter for those with GenBank files
    table_df = fetch_ship_table(curated=True)
    if table_df is not None and not table_df.empty:
        # Filter the DataFrame to only include rows where starshipID is in available_starships
        table_df = table_df[table_df['starshipID'].isin(available_starships)]
        
        if table_df.empty:
            logger.error("No matching entries found with GenBank files")
            return html.Div("No Starships found with corresponding GenBank files")
        
        table = make_pgv_table(
            df=table_df,
            columns=table_columns,
            id="pgv-table",
            select_rows=True,
            pg_sz=10,
        )
        logger.info(f"Table created successfully with {len(table_df)} entries")
        return table
    
    logger.error("No data available or empty DataFrame")
    return html.Div("No data available")

@callback(
    [Output("pgv-figure", "children"), Output("pgv-message", "children")],
    [Input("update-button", "n_clicks")],
    [
        State("pgv-table", "selected_rows"),
        State("pgv-table", "data"),
        State("identity-threshold", "value"),
        State("use-file-order", "value"),
    ],
)
def update_visualization(n_clicks, selected_rows, table_data, identity_threshold, use_file_order):
    from src.utils.clinker import process_gbk_files, create_clustermap_data, custom_save_html
    import tempfile
    
    if not n_clicks:
        return (
            no_update,
            "Select Starships from the table and click 'Show Selected Starships'",
        )

    if not selected_rows or len(selected_rows) == 0:
        return html.H4("Please select at least one Starship."), None

    if len(selected_rows) > 4:
        return html.H4("Please select no more than 4 Starships."), None

    try:
        # Get the selected starshipIDs
        selected_ships = [table_data[idx]["starshipID"] for idx in selected_rows]
        
        # Create temporary directory for visualization
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as tmp_file:
            # Process the selected GenBank files
            gbk_dir = Path("/home/adrian/Systematics/Starship_Database/starbase/src/database/db/ships/gbks")
            selected_gbks = [str(next(gbk_dir.glob(f"*_{ship}.gbk"))) for ship in selected_ships]
            
            clusters = process_gbk_files(selected_gbks)
            if not clusters:
                return html.H4("Error processing GenBank files."), None

            # Create visualization
            viz_html = custom_save_html(
                clusters,
                identity_threshold or 0.3,
                use_file_order or False,
                tmp_file.name
            )

            return html.Iframe(
                srcDoc=viz_html,
                style={
                    "width": "100%",
                    "height": "800px",
                    "border": "none",
                    "overflow": "hidden"
                }
            ), None

    except Exception as e:
        logger.error(f"Error in visualization: {str(e)}")
        return html.H4(f"Error: {str(e)}"), None

import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update, clientside_callback, ClientsideFunction
from dash.dependencies import Output, Input, State
import json
import tempfile
from pathlib import Path

from src.components.callbacks import create_modal_callback
from src.components.error_boundary import handle_callback_error
from src.config.logging import get_logger
from src.config.settings import GBK_PATH

logger = get_logger(__name__)

dash.register_page(__name__, path="/synteny")

# TODO: handle multiple ships grouped under a single accession_tag in gff files
# TODO: table should include information about multiple ships under a single accession_tag

table_columns = [
    {
        "id": "accession_display",
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

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="synteny-url", refresh=False),
        dcc.Store(id="synteny-data-store"),
        
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title("Starship Synteny Viewer", order=1, mb="md"),
                dmc.Text(
                    "Interactive synteny visualization using ClusterMap.js - Compare and visualize up to 4 Starship sequences",
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
                # Left Column - Controls
                dmc.GridCol(
                    span={"base": 12, "md": 4},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                dmc.Group(
                                    pos="apart",
                                    children=[
                                        dmc.Title("Select Starships", order=2),
                                        dmc.Button(
                                            dmc.Text("Generate Visualization", size="lg"),
                                            id="synteny-update-button",
                                            variant="gradient",
                                            gradient={"from": "indigo", "to": "cyan"},
                                            leftSection=html.I(className="bi bi-diagram-3"),
                                        ),
                                    ],
                                ),
                                
                                # Configuration Controls
                                dmc.Accordion(
                                    variant="filled",
                                    chevronPosition="right",
                                    chevronSize=16,
                                    children=[                                        
                                    dmc.AccordionItem(
                                        value="visualization-configuration",
                                        children=[
                                                    dmc.Group(
                                                        pos="apart",
                                                        children=[
                                                            dmc.Stack([
                                                                dmc.NumberInput(
                                                                    label="Identity Threshold",
                                                                    id="synteny-identity-threshold",
                                                                    value=0.3,
                                                                    min=0.1,
                                                                    max=1,
                                                                    step=0.1,
                                                                    description="Minimum identity for showing links",
                                                                ),
                                                                dmc.NumberInput(
                                                                    label="Scale Factor",
                                                                    id="synteny-scale-factor",
                                                                    value=10,
                                                                    min=1,
                                                                    max=50,
                                                                    description="Scaling factor for visualization",
                                                                ),
                                                                dmc.NumberInput(
                                                                    label="Cluster Spacing",
                                                                    id="synteny-cluster-spacing",
                                                                    value=50,
                                                                    min=10,
                                                                    max=200,
                                                                    description="Vertical spacing between clusters",
                                                                ),
                                                                dmc.Switch(
                                                                    id="synteny-use-file-order",
                                                                    label="Use File Order",
                                                                    checked=False,
                                                                    description="Maintain original file order",
                                                                ),
                                                                dmc.Switch(
                                                                    id="synteny-show-links",
                                                                    label="Show Links",
                                                                    checked=True,
                                                                    description="Display synteny links",
                                                                ),
                                                                dmc.Switch(
                                                                    id="synteny-show-gene-labels",
                                                                    label="Show Gene Labels",
                                                                    checked=False,
                                                                    description="Display gene labels",
                                                                ),
                                                            ], gap="sm")
                                                        ]
                                                    )   
                                                ],
                                            ),
                                        ],
                                    ),
                                # Table
                                dcc.Loading(
                                    id="synteny-loading",
                                    type="circle",
                                    children=html.Div(id="synteny-table"),
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
                    span={"base": 12, "md": 8},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                html.Div(id="synteny-message", style={"textAlign": "center"}),
                                dcc.Loading(
                                    id="synteny-loading-viz",
                                    type="circle",
                                    children=[
                                        html.Div(
                                            id="synteny-visualization",
                                            style={
                                                "height": "800px",
                                                "width": "100%",
                                                "overflow": "auto",
                                                "backgroundColor": "#f8f9fa",
                                                "border": "1px solid #dee2e6",
                                                "borderRadius": "4px",
                                            },
                                        ),
                                        # Controls for the visualization
                                        dmc.Group([
                                            dmc.Button(
                                                "Save as SVG",
                                                id="synteny-save-svg",
                                                variant="light",
                                                leftSection=html.I(className="bi bi-download"),
                                                style={"display": "none"}
                                            ),
                                            dmc.Button(
                                                "Reset View",
                                                id="synteny-reset-view", 
                                                variant="light",
                                                leftSection=html.I(className="bi bi-arrow-clockwise"),
                                                style={"display": "none"}
                                            ),
                                        ], justify="center", mt="md", id="synteny-controls")
                                    ]
                                ),
                            ], gap="md"),
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
    Output("synteny-table", "children"),
    Input("synteny-url", "href"),
    prevent_initial_call=False,
)
@handle_callback_error
def load_synteny_table(href):
    """Load and display the ship selection table, filtered to only show entries with GenBank files"""
    from src.database.sql_manager import fetch_ship_table
    from src.components.tables import make_pgv_table  
    
    # Get list of available GenBank files and extract accession_displays (with version)
    available_starships = {path.stem for path in Path(GBK_PATH).glob("*.gbk")}
    logger.info(f"Found {len(available_starships)} GenBank files: {list(available_starships)[:5]}...")
    
    # Fetch all ships and filter for those with GenBank files
    table_df = fetch_ship_table(curated=True)
    if table_df is not None and not table_df.empty:
        logger.info(f"Fetched {len(table_df)} ships from database")
        logger.info(f"Sample accession_tags from database: {table_df['accession_tag'].head().tolist()}")
        logger.info(f"Sample accession_displays from database: {table_df['accession_display'].head().tolist()}")
        
        # Filter the DataFrame to only include rows where accession_display is in available_starships
        table_df = table_df[table_df['accession_display'].isin(available_starships)]
        
        if table_df.empty:
            logger.error("No matching entries found with GenBank files")
            return html.Div("No Starships found with corresponding GenBank files")
        
        table = make_pgv_table(
            df=table_df,
            columns=table_columns,
            id="synteny-table",
            select_rows=True,
            pg_sz=10,
        )
        logger.info(f"Table created successfully with {len(table_df)} entries")
        return table
    
    logger.error("No data available or empty DataFrame")
    return html.Div("No data available")

@callback(
    Output("synteny-data-store", "data"),
    [Input("synteny-update-button", "n_clicks")],
    [
        State("synteny-table", "selectedRows"),
        State("synteny-table", "rowData"),
        State("synteny-identity-threshold", "value"),
        State("synteny-use-file-order", "value"),
    ],
)
@handle_callback_error
def prepare_synteny_data(n_clicks, selected_rows, table_data, identity_threshold, use_file_order):
    """Prepare the data for visualization"""
    from src.utils.clinker import process_gbk_files, create_clustermap_data
    
    if not n_clicks:
        return no_update

    logger.info(f"Button clicked: {n_clicks}, Selected rows: {selected_rows}")
    logger.info(f"Selected rows type: {type(selected_rows)}, Length: {len(selected_rows) if selected_rows else 0}")

    if not selected_rows or len(selected_rows) == 0:
        return {"error": "Please select at least one Starship."}

    if len(selected_rows) > 4:
        return {"error": "Please select no more than 4 Starships."}

    try:
        # Get the selected accession_displays (AG Grid returns the actual row data, not indices)
        selected_ships = [row["accession_display"] for row in selected_rows]
        
        # Process the selected GenBank files
        selected_gbks = []
        
        for ship in selected_ships:
            gbk_file = Path(GBK_PATH) / f"{ship}.gbk"
            if gbk_file.exists():
                selected_gbks.append(str(gbk_file))
            else:
                logger.warning(f"No GenBank file found for ship: {ship}")
        
        if not selected_gbks:
            return {"error": "No GenBank files found for selected ships."}
            
        globaligner = process_gbk_files(selected_gbks)
        if not globaligner:
            return {"error": "Error processing GenBank files."}

        # Create the data for clustermap.js
        clustermap_data = create_clustermap_data(globaligner, use_file_order or False)
        
        # Add configuration
        config = {
            "identity_threshold": identity_threshold or 0.3,
            "use_file_order": use_file_order or False,
        }
        
        return {
            "data": clustermap_data,
            "config": config,
            "success": True
        }

    except Exception as e:
        logger.error(f"Error in preparing synteny data: {str(e)}")
        return {"error": f"Error: {str(e)}"}

# Client-side callback to render the visualization
clientside_callback(
    """
    function(store_data, scale_factor, cluster_spacing, show_links, show_gene_labels, identity_threshold) {
        if (!store_data || store_data.error) {
            return [
                store_data && store_data.error ? store_data.error : "Select Starships and click 'Generate Visualization'",
                {"display": "none"}
            ];
        }
        
        if (!store_data.success || !store_data.data) {
            return ["No data available", {"display": "none"}];
        }
        
        // Check if synteny helper functions are available
        if (typeof window.syntenyViz === 'undefined') {
            console.error("Synteny visualization helpers not loaded");
            return ["Error: Visualization library not loaded", {"display": "none"}];
        }
        
        try {
            // Configuration for the visualization
            const config = {
                plot: {
                    scaleFactor: scale_factor || 30,
                    scaleGenes: true,
                },
                cluster: {
                    spacing: cluster_spacing || 50,
                    alignLabels: true,
                },
                gene: {
                    label: {
                        show: show_gene_labels || false,
                    }
                },
                legend: {
                    show: true,
                },
                link: {
                    show: show_links !== false,
                    threshold: identity_threshold || 0.3,
                }
            };
            
            // Initialize the ClusterMap visualization
            const chart = window.syntenyViz.initialize('synteny-visualization', store_data.data, config);
            
            if (chart) {
                return ["", {"display": "block"}];
            } else {
                return ["Error initializing visualization", {"display": "none"}];
            }
        } catch (error) {
            console.error("Error in synteny visualization:", error);
            return [`Error: ${error.message}`, {"display": "none"}];
        }
    }
    """,
    [Output("synteny-message", "children"), Output("synteny-controls", "style")],
    [
        Input("synteny-data-store", "data"),
        Input("synteny-scale-factor", "value"),
        Input("synteny-cluster-spacing", "value"),
        Input("synteny-show-links", "value"),
        Input("synteny-show-gene-labels", "value"),
        Input("synteny-identity-threshold", "value"),
    ]
)

# Client-side callback for saving SVG
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks && typeof window.syntenyViz !== 'undefined') {
            try {
                window.syntenyViz.exportSVG('synteny-visualization', 'starship-synteny.svg');
            } catch (error) {
                console.error("Error exporting SVG:", error);
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("synteny-save-svg", "n_clicks"),
    Input("synteny-save-svg", "n_clicks"),
    prevent_initial_call=True
)

# Client-side callback for resetting view
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks && typeof window.syntenyViz !== 'undefined') {
            try {
                window.syntenyViz.resetView('synteny-visualization');
            } catch (error) {
                console.error("Error resetting view:", error);
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("synteny-reset-view", "n_clicks"),
    Input("synteny-reset-view", "n_clicks"),
    prevent_initial_call=True
)

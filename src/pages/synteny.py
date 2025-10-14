import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update, clientside_callback, ClientsideFunction
from dash.dependencies import Output, Input, State
import json
import tempfile
from pathlib import Path

from src.components.callbacks import create_modal_callback, handle_callback_error
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
                    "Interactive synteny visualization using ClusterMap.js - Compare and visualize up to 4 Starship sequences with GFF annotation data",
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
                                dmc.Stack([
                                    dmc.Title("Select Starships", order=2),
                                    dcc.Loading(
                                        id="synteny-loading",
                                        type="circle",
                                        children=html.Div(id="synteny-table"),
                                    ),
                                ], gap="sm"),                                
                                
                                # Configuration Controls
                                dmc.Accordion(
                                    variant="filled",
                                    chevronPosition="right",
                                    chevronSize=16,
                                    value="visualization-settings",
                                    children=[
                                        dmc.AccordionItem(
                                            value="visualization-settings",
                                            children=[
                                                dmc.AccordionControl(dmc.Title("Visualization Settings", order=3)),
                                                dmc.AccordionPanel([
                                                    dmc.Stack([
                                                        dmc.NumberInput(
                                                            label="Identity Threshold",
                                                            id="synteny-identity-threshold",
                                                            value=0.3,
                                                            min=0.1,
                                                            max=1,
                                                            step=0.1,
                                                        ),
                                                        dmc.NumberInput(
                                                            label="Scale Factor",
                                                            id="synteny-scale-factor",
                                                            value=10,
                                                            min=1,
                                                            max=50,
                                                        ),
                                                        dmc.NumberInput(
                                                            label="Cluster Spacing",
                                                            id="synteny-cluster-spacing",
                                                            value=50,
                                                            min=10,
                                                            max=200,
                                                        ),
                                                        dmc.Switch(
                                                            id="synteny-use-file-order",
                                                            label="Maintain original file order",
                                                            checked=False,
                                                        ),
                                                        dmc.Switch(
                                                            id="synteny-show-links",
                                                            label="Display synteny links",
                                                            checked=True,
                                                        ),
                                                        dmc.Switch(
                                                            id="synteny-show-gene-labels",
                                                            label="Display gene labels",
                                                            checked=False,
                                                        ),
                                                    ], gap="sm")
                                                ])
                                            ],
                                        ),
                                    ],
                                ),
                                dmc.Group([
                                    dmc.Button(
                                        dmc.Text("Generate Visualization", size="lg"),
                                        id="synteny-update-button",
                                        variant="gradient",
                                        gradient={"from": "indigo", "to": "cyan"},
                                        leftSection=html.I(className="bi bi-diagram-3")
                                    ),
                                ], justify="center", gap="md", grow=False),
                            ], gap="md"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            h="100%",
                        ),
                    ],
                ),
                
                # Right Column - Static Visualization Container
                dmc.GridCol(
                    span={"base": 12, "md": 8},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack([
                                html.Div(id="synteny-message"),
                                dcc.Loading(
                                    id="synteny-loading-viz",
                                    type="default",
                                    fullscreen=False,
                                    children=[
                                        html.Div(id="synteny-loading-trigger"),  # Hidden div to trigger loading
                                        # STATIC CONTAINER - never changes ID
                                        html.Div(
                                            "Select Starships with GFF annotation data and click 'Generate Visualization'",
                                            id="synteny-static-viz",  # Fixed ID
                                            style={
                                                "height": "800px",
                                                "width": "100%",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "justifyContent": "center",
                                                "color": "#6c757d",
                                                "fontSize": "16px",
                                                "backgroundColor": "#f8f9fa",
                                                "border": "1px solid #dee2e6",
                                                "borderRadius": "4px",
                                            },
                                        ),
                                        # Controls
                                        dmc.Group([
                                            dmc.Button(
                                                "Save as SVG",
                                                id="synteny-save-svg",
                                                variant="light",
                                                color="green",
                                                leftSection=html.I(className="bi bi-download"),
                                            ),
                                            dmc.Button(
                                                "Reset View",
                                                id="synteny-reset-view", 
                                                variant="light",
                                                color="blue",
                                                leftSection=html.I(className="bi bi-arrow-clockwise"),
                                            ),
                                            dmc.Button(
                                                "Clear Visualization",
                                                id="synteny-clear-viz",
                                                variant="light",
                                                color="red",
                                                leftSection=html.I(className="bi bi-x-circle"),
                                            ),
                                        ], justify="center", gap="md", grow=True, id="synteny-controls")
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
    """Load and display the ship selection table, filtered to only show entries with GFF annotation data"""
    from src.database.sql_manager import fetch_ship_table
    from src.components.tables import make_pgv_table  
    
    # Fetch ships that have GFF annotation data
    table_df = fetch_ship_table(curated=False, with_sequence=True, with_gff_entries=True)
    if table_df is not None and not table_df.empty:
        logger.info(f"Fetched {len(table_df)} ships with GFF annotation data from database")
        logger.info(f"Sample accession_tags from database: {table_df['accession_tag'].head().tolist()}")
        logger.info(f"Sample accession_displays from database: {table_df['accession_display'].head().tolist()}")
        
        table = make_pgv_table(
            df=table_df,
            columns=table_columns,
            id="synteny-table-grid",
            select_rows=True,
            pg_sz=10,
        )
        logger.info(f"Table created successfully with {len(table_df)} entries")
        return table
    
    logger.error("No data available or empty DataFrame")
    return html.Div("No Starships found with GFF annotation data")

@callback(
    Output("synteny-data-store", "data"),
    [Input("synteny-update-button", "n_clicks")],
    [
        State("synteny-table-grid", "selectedRows"),
        State("synteny-table-grid", "rowData"),
        State("synteny-identity-threshold", "value"),
        State("synteny-use-file-order", "value"),
    ],
)
@handle_callback_error
def prepare_synteny_data(n_clicks, selected_rows, table_data, identity_threshold, use_file_order):
    """Prepare the data for visualization with enhanced error handling"""
    from src.utils.clinker import process_gbk_files, create_clustermap_data
    import time
    
    if not n_clicks:
        return no_update
    
    logger.info(f"Button clicked: {n_clicks}, Selected rows: {selected_rows}")
    
    if not selected_rows or len(selected_rows) == 0:
        return {"error": "Please select at least one Starship."}

    if len(selected_rows) > 4:
        return {"error": "Please select no more than 4 Starships."}

    try:
        # Extract base accession_tag (not accession_display) for database lookup
        selected_ships = [row["accession_tag"] for row in selected_rows]
        selected_displays = [row.get("accession_display", row["accession_tag"]) for row in selected_rows]
        
        # Generate GenBank files on-the-fly from GFF data
        logger.info(f"Generating GenBank files for selected ships: {selected_ships} (display: {selected_displays})")
        globaligner = process_gbk_files(None, accession_tags=selected_ships)
        
        if not globaligner:
            return {
                "error": "Failed to process GenBank files from GFF data.",
                "reason": "No valid clusters generated",
                "details": {"accessions": selected_ships}
            }

        clustermap_data = create_clustermap_data(globaligner, use_file_order or False)
        
        # Ensure data is serializable and clean
        clean_data = json.loads(json.dumps(clustermap_data, default=str))
        
        config = {
            "identity_threshold": identity_threshold or 0.3,
            "use_file_order": use_file_order or False,
        }
        
        logger.info(f"Successfully created clustermap data with {len(clean_data.get('clusters', []))} clusters and {len(clean_data.get('links', []))} links")
        
        return {
            "data": clean_data,
            "config": config,
            "success": True,
            "timestamp": time.time()  # Add timestamp for cache busting
        }

    except ValueError as e:
        # Validation errors - these are user-actionable
        logger.warning(f"Validation error in synteny data: {str(e)}")
        return {
            "error": f"Data validation error: {str(e)}",
            "reason": "Invalid or missing data",
            "details": {"accessions": selected_ships}
        }
    except Exception as e:
        # Unexpected errors - these need investigation
        logger.error(f"Unexpected error in preparing synteny data: {str(e)}", exc_info=True)
        return {
            "error": f"Processing error: {str(e)}",
            "reason": "Unexpected failure during visualization generation",
            "details": {"accessions": selected_ships}
        }

@callback(
    Output("synteny-loading-trigger", "children"),
    Input("synteny-data-store", "data"),
)
def update_loading_trigger(data):
    """Update loading trigger to ensure loading indicator shows"""
    if data:
        return ""
    return ""
# Remove this callback as it references non-existent IDs and conflicts with the static container approach


# Client-side callback to render the visualization
clientside_callback(
    """
    function(store_data, scale_factor, cluster_spacing, show_links, show_gene_labels, identity_threshold, clear_clicks) {
        // Always use the same container ID
        const containerId = 'synteny-static-viz';
        
        // If clear button was clicked, just clear and return
        if (clear_clicks) {
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = 'Select Starships with GFF annotation data and click "Generate Visualization"';
                container.style.display = 'flex';
                container.style.alignItems = 'center';
                container.style.justifyContent = 'center';
                container.style.color = '#6c757d';
                container.style.fontSize = '16px';
            }
            return [None, {"display": "none"}];
        }
        
        // Handle null/cleared data store
        if (!store_data) {
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = 'Select Starships with GFF annotation data and click "Generate Visualization"';
                container.style.display = 'flex';
                container.style.alignItems = 'center';
                container.style.justifyContent = 'center';
                container.style.color = '#6c757d';
                container.style.fontSize = '16px';
            }
            return ["", {"display": "none"}];
        }
        
        if (store_data.error) {
            const container = document.getElementById(containerId);
            if (container) {
                let errorHtml = `<div style="color: red; padding: 20px; text-align: center;">
                    <strong>${store_data.error}</strong>`;
                if (store_data.reason) {
                    errorHtml += `<br><small>${store_data.reason}</small>`;
                }
                if (store_data.details) {
                    errorHtml += `<br><small style="font-family: monospace;">Details: ${JSON.stringify(store_data.details)}</small>`;
                }
                errorHtml += `</div>`;
                container.innerHTML = errorHtml;
            }
            return [store_data.error, {"display": "none"}];
        }
        
        if (!store_data.success || !store_data.data) {
            return ["No data available", {"display": "none"}];
        }
        
        try {
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
            
            const container = document.getElementById(containerId);
            if (!container) {
                console.error("Static container not found");
                return ["Error: Container not found", {"display": "none"}];
            }
            
            // COMPLETE RESET - clear everything
            container.innerHTML = '';
            container.style.display = 'block';
            container.style.alignItems = 'unset';
            container.style.justifyContent = 'unset';
            container.style.color = 'unset';
            container.style.fontSize = 'unset';
            
            // Clear any global D3 selections that might reference this container
            if (typeof d3 !== 'undefined') {
                d3.selectAll('#' + containerId + ' *').remove();
            }
            
            // Wait a bit, then create new visualization
            setTimeout(function() {
                try {
                    const chart = ClusterMap.ClusterMap().config(config);
                    d3.select('#' + containerId)
                        .datum(store_data.data)
                        .call(chart);
                    
                    console.log("ClusterMap created successfully");
                } catch (error) {
                    console.error("Error creating ClusterMap:", error);
                    container.innerHTML = '<div style="color: red; padding: 20px; text-align: center;">Error: ' + error.message + '</div>';
                }
            }, 200);
            
            return ["", {"display": "block"}];
            
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
        Input("synteny-clear-viz", "n_clicks"),  # Add clear button input
    ]
)


clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) {
            // Clear the data store to allow new visualizations
            return null;
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("synteny-data-store", "data", allow_duplicate=True),
    Input("synteny-clear-viz", "n_clicks"),
    prevent_initial_call=True
)

# Client-side callback for saving SVG
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) {
            try {
                // Find the SVG in the static visualization container
                const svg = document.querySelector('#synteny-static-viz svg');
                if (svg) {
                    // Create a blob and download
                    const serializer = new XMLSerializer();
                    const svgString = serializer.serializeToString(svg);
                    const blob = new Blob([svgString], {type: 'image/svg+xml'});
                    
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement('a');
                    link.href = url;
                    link.download = 'starship-synteny.svg';
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(url);
                } else {
                    console.warn("No SVG found to export");
                }
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
        if (n_clicks) {
            try {
                // Find the SVG and reset its transform in the static visualization container
                const svg = document.querySelector('#synteny-static-viz svg');
                if (svg) {
                    const g = svg.querySelector('g');
                    if (g) {
                        g.setAttribute('transform', 'translate(0,0) scale(1)');
                    }
                } else {
                    console.warn("No SVG found to reset");
                }
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

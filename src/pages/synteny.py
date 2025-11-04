import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update, clientside_callback, ClientsideFunction
from dash.dependencies import Output, Input, State
from dash_iconify import DashIconify
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
        dcc.Store(id="synteny-loading-state", data=False),  # Track loading state
        
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
                                    dmc.Group([
                                        dmc.Title("Select Starships", order=2),
                                    ], justify="space-between", align="center"),
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
                                                        dmc.Group([
                                                            dmc.NumberInput(
                                                                label="Identity Threshold",
                                                                id="synteny-identity-threshold",
                                                                value=0.3,
                                                                min=0.1,
                                                                max=1,
                                                                step=0.1,
                                                                w=200,
                                                            ),
                                                            dmc.Tooltip(
                                                                label="Minimum sequence identity for showing synteny links (0.1-1.0)",
                                                                children=dmc.ActionIcon(
                                                                    DashIconify(icon="mdi:help-circle"),
                                                                    variant="transparent",
                                                                    size="sm"
                                                                )
                                                            )
                                                        ], align="center"),
                                                        dmc.Group([
                                                            dmc.NumberInput(
                                                                label="Scale Factor",
                                                                id="synteny-scale-factor",
                                                                value=30,
                                                                min=5,
                                                                max=100,
                                                                w=200,
                                                            ),
                                                            dmc.Tooltip(
                                                                label="Zoom level for gene visualization (higher = more detail)",
                                                                children=dmc.ActionIcon(
                                                                    DashIconify(icon="mdi:help-circle"),
                                                                    variant="transparent",
                                                                    size="sm"
                                                                )
                                                            )
                                                        ], align="center"),
                                                        dmc.Group([
                                                            dmc.NumberInput(
                                                                label="Cluster Spacing",
                                                                id="synteny-cluster-spacing",
                                                                value=50,
                                                                min=20,
                                                                max=200,
                                                                w=200,
                                                            ),
                                                            dmc.Tooltip(
                                                                label="Vertical spacing between different Starships",
                                                                children=dmc.ActionIcon(
                                                                    DashIconify(icon="mdi:help-circle"),
                                                                    variant="transparent",
                                                                    size="sm"
                                                                )
                                                            )
                                                        ], align="center"),
                                                        dmc.Tooltip(
                                                            label="Preserve the order of selected Starships instead of clustering by similarity",
                                                            children=dmc.Switch(
                                                                id="synteny-use-file-order",
                                                                label="Maintain original file order",
                                                                checked=False,
                                                            )
                                                        ),
                                                        dmc.Tooltip(
                                                            label="Show/hide the colored lines connecting similar genes between Starships",
                                                            children=dmc.Switch(
                                                                id="synteny-show-links",
                                                                label="Display synteny links",
                                                                checked=True,
                                                            )
                                                        ),
                                                        dmc.Tooltip(
                                                            label="Show/hide gene names on the visualization (can be cluttered)",
                                                            children=dmc.Switch(
                                                                id="synteny-show-gene-labels",
                                                                label="Display gene labels",
                                                                checked=False,
                                                            )
                                                        ),
                                                    ], gap="md")
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
                                dmc.Group([
                                    html.Div(id="synteny-message"),
                                    dmc.Badge(
                                        "Ready",
                                        id="synteny-status-badge",
                                        color="gray",
                                        variant="light",
                                        size="sm"
                                    )
                                ], justify="space-between", align="center"),
                                dcc.Loading(
                                    id="synteny-loading-viz",
                                    type="circle",
                                    children=[
                                        # STATIC CONTAINER - never changes ID
                                        html.Div(id="synteny-viz-message"),
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
                                                "position": "relative",  # Added for centering
                                            },
                                        ),
                                    ]
                                ),
                                # Controls - moved outside loading indicator
                                dmc.Group([
                                    dmc.Tooltip(
                                        label="Download the current visualization as an SVG file",
                                        children=dmc.Button(
                                            "Save as SVG",
                                            id="synteny-save-svg",
                                            variant="light",
                                            color="green",
                                            leftSection=html.I(className="bi bi-download"),
                                            disabled=True,  # Initially disabled
                                        )
                                    ),
                                    dmc.Tooltip(
                                        label="Reset zoom and pan to the default view",
                                        children=dmc.Button(
                                            "Reset View",
                                            id="synteny-reset-view",
                                            variant="light",
                                            color="blue",
                                            leftSection=html.I(className="bi bi-arrow-clockwise"),
                                            disabled=True,  # Initially disabled
                                        )
                                    ),
                                    dmc.Tooltip(
                                        label="Clear the current visualization and start over",
                                        children=dmc.Button(
                                            "Clear Visualization",
                                            id="synteny-clear-viz",
                                            variant="light",
                                            color="red",
                                            leftSection=html.I(className="bi bi-x-circle"),
                                        )
                                    ),
                                ], justify="center", gap="md", grow=True, id="synteny-controls", style={"display": "none"}),  # Initially hidden
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
            id="synteny-table",
            select_rows=True,
            pg_sz=10,
        )
        logger.info(f"Table created successfully with {len(table_df)} entries")
        return table
    
    logger.error("No data available or empty DataFrame")
    return html.Div("No Starships found with GFF annotation data")

@callback(
    [
        Output("synteny-data-store", "data"),
        Output("synteny-loading-state", "data"),
        Output("synteny-message", "children"),
        Output("synteny-update-button", "disabled"),
        Output("synteny-update-button", "children"),
    ],
    [Input("synteny-update-button", "n_clicks")],
    [
        State("synteny-table", "selectedRows"),
        State("synteny-identity-threshold", "value"),
        State("synteny-use-file-order", "value"),
    ],
)
@handle_callback_error
def prepare_synteny_data(n_clicks, selected_rows, identity_threshold, use_file_order):
    """Prepare the data for visualization with enhanced error handling"""
    from src.utils.clinker import process_gbk_files, create_clustermap_data
    import time

    if not n_clicks:
        return no_update, False, "", False, dmc.Text("Generate Visualization", size="lg")

    # Input validation
    if not selected_rows or len(selected_rows) == 0:
        return ({"error": "Please select at least one Starship."}, False,
                "Please select at least one Starship.", False,
                dmc.Text("Generate Visualization", size="lg"))

    if len(selected_rows) > 4:
        return ({"error": "Please select no more than 4 Starships."}, False,
                "Please select no more than 4 Starships.", False,
                dmc.Text("Generate Visualization", size="lg"))

    # Validate identity threshold
    if identity_threshold is not None:
        if not isinstance(identity_threshold, (int, float)) or identity_threshold < 0.1 or identity_threshold > 1.0:
            return ({"error": "Identity threshold must be between 0.1 and 1.0"}, False,
                    "Identity threshold must be between 0.1 and 1.0", False,
                    dmc.Text("Generate Visualization", size="lg"))

    logger.info(f"Button clicked: {n_clicks}, Selected rows: {selected_rows}")

    try:
        selected_ships = [row["accession_display"] for row in selected_rows]

        # Validate selected ships
        if not all(isinstance(ship, str) and ship.strip() for ship in selected_ships):
            return ({"error": "Invalid ship selection data"}, False,
                    "Invalid ship selection data", False,
                    dmc.Text("Generate Visualization", size="lg"))

        # Start loading - disable button and show loading state
        logger.info(f"Generating GenBank files for selected ships: {selected_ships}")

        globaligner = process_gbk_files(None, accession_tags=selected_ships)

        if not globaligner:
            return ({"error": "Error processing GenBank files from GFF data."}, False,
                    "Error processing GenBank files from GFF data.", False,
                    dmc.Text("Generate Visualization", size="lg"))

        clustermap_data = create_clustermap_data(globaligner, use_file_order or False)

        # Ensure data is serializable and clean
        clean_data = json.loads(json.dumps(clustermap_data, default=str))

        config = {
            "identity_threshold": identity_threshold or 0.3,
            "use_file_order": use_file_order or False,
        }

        result = {
            "data": clean_data,
            "config": config,
            "success": True,
            "timestamp": time.time()  # Add timestamp for cache busting
        }

        return (result, False, f"Successfully processed {len(selected_ships)} Starships", False,
                dmc.Text("Generate Visualization", size="lg"))

    except Exception as e:
        logger.error(f"Error in preparing synteny data: {str(e)}")
        error_msg = f"Error processing data: {str(e)}"
        return ({"error": error_msg}, False, error_msg, False,
                dmc.Text("Generate Visualization", size="lg"))
# Remove this callback as it references non-existent IDs and conflicts with the static container approach


# Client-side callback to render the visualization
clientside_callback(
    """
    function(store_data, loading_state, scale_factor, cluster_spacing, show_links, show_gene_labels, identity_threshold, clear_clicks) {
        // Clear visualization container immediately when data changes
        const container = document.getElementById('synteny-static-viz');
        if (container) {
            if (store_data && store_data.success) {
                // Clear for successful data - will be filled by rendering
                container.innerHTML = '';
                container.style.display = 'block';
                container.style.alignItems = 'unset';
                container.style.justifyContent = 'unset';
                container.style.color = 'unset';
                container.style.fontSize = 'unset';
                container.style.backgroundColor = 'transparent';
                container.style.border = 'none';
                container.style.borderRadius = 'unset';
                container.style.position = 'relative';
            }
        }
        const containerId = 'synteny-static-viz';
        const controlsId = 'synteny-controls';
        const saveBtnId = 'synteny-save-svg';
        const resetBtnId = 'synteny-reset-view';

        // If clear button was clicked, clear everything and reset
        if (clear_clicks) {
            // Clear container and reset to default state
            const clearContainer = document.getElementById(containerId);
            if (clearContainer) {
                clearContainer.innerHTML = 'Select Starships with GFF annotation data and click "Generate Visualization"';
                clearContainer.style.display = 'flex';
                clearContainer.style.alignItems = 'center';
                clearContainer.style.justifyContent = 'center';
                clearContainer.style.color = '#6c757d';
                clearContainer.style.fontSize = '16px';
                clearContainer.style.backgroundColor = '#f8f9fa';
                clearContainer.style.border = '1px solid #dee2e6';
                clearContainer.style.borderRadius = '4px';
            }

            // Hide and disable controls
            const controls = document.getElementById(controlsId);
            if (controls) controls.style.display = 'none';

            const saveBtn = document.getElementById(saveBtnId);
            if (saveBtn) saveBtn.disabled = true;

            const resetBtn = document.getElementById(resetBtnId);
            if (resetBtn) resetBtn.disabled = true;

            return ["Visualization cleared", {"display": "none"}, true, true];
        }

        // Handle null/cleared data store
        if (!store_data) {
            // Clear container and reset to default state
            const nullContainer = document.getElementById(containerId);
            if (nullContainer) {
                nullContainer.innerHTML = 'Select Starships with GFF annotation data and click "Generate Visualization"';
                nullContainer.style.display = 'flex';
                nullContainer.style.alignItems = 'center';
                nullContainer.style.justifyContent = 'center';
                nullContainer.style.color = '#6c757d';
                nullContainer.style.fontSize = '16px';
                nullContainer.style.backgroundColor = '#f8f9fa';
                nullContainer.style.border = '1px solid #dee2e6';
                nullContainer.style.borderRadius = '4px';
            }

            // Hide and disable controls
            const controls = document.getElementById(controlsId);
            if (controls) controls.style.display = 'none';

            const saveBtn = document.getElementById(saveBtnId);
            if (saveBtn) saveBtn.disabled = true;

            const resetBtn = document.getElementById(resetBtnId);
            if (resetBtn) resetBtn.disabled = true;

            return ["", {"display": "none"}, true, true];
        }

        // Handle errors
        if (store_data.error) {
            // Clear container and show error
            const errorContainer = document.getElementById(containerId);
            if (errorContainer) {
                errorContainer.innerHTML = '<div style="color: red; padding: 20px; text-align: center; background: #fff5f5; border: 1px solid #feb2b2; border-radius: 4px;"><strong>Error:</strong><br>' + store_data.error + '</div>';
                errorContainer.style.display = 'flex';
                errorContainer.style.alignItems = 'center';
                errorContainer.style.justifyContent = 'center';
                errorContainer.style.backgroundColor = 'transparent';
                errorContainer.style.border = 'none';
                errorContainer.style.borderRadius = 'unset';
            }

            // Hide and disable controls
            const controls = document.getElementById(controlsId);
            if (controls) controls.style.display = 'none';

            const saveBtn = document.getElementById(saveBtnId);
            if (saveBtn) saveBtn.disabled = true;

            const resetBtn = document.getElementById(resetBtnId);
            if (resetBtn) resetBtn.disabled = true;

            return [store_data.error, {"display": "none"}, true, true];
        }

        if (!store_data.success || !store_data.data) {
            return ["No data available", {"display": "none"}, true, true];
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

            const renderContainer = document.getElementById(containerId);
            if (!renderContainer) {
                console.error("Static container not found");
                return ["Error: Container not found", {"display": "none"}, true, true];
            }

            // Show loading state if still processing
            if (loading_state) {
                // Clear container and show loading state
                renderContainer.innerHTML = '<div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%;"><div class="synteny-loading-spinner"></div><div style="color: #6c757d; font-weight: 500; margin-top: 20px;">Rendering visualization...</div></div>';
                renderContainer.style.display = 'flex';
                renderContainer.style.alignItems = 'center';
                renderContainer.style.justifyContent = 'center';
                renderContainer.style.backgroundColor = '#f8f9fa';
                renderContainer.style.border = '1px solid #dee2e6';
                renderContainer.style.borderRadius = '4px';

                // Hide controls during loading
                const controls = document.getElementById(controlsId);
                if (controls) controls.style.display = 'none';

                return ["Rendering visualization...", {"display": "none"}, false, false];
            }

            // Clear any previous content and prepare for rendering
            renderContainer.innerHTML = '';
            renderContainer.style.display = 'block';
            renderContainer.style.alignItems = 'unset';
            renderContainer.style.justifyContent = 'unset';
            renderContainer.style.color = 'unset';
            renderContainer.style.fontSize = 'unset';
            renderContainer.style.backgroundColor = 'transparent';
            renderContainer.style.border = 'none';
            renderContainer.style.borderRadius = 'unset';
            renderContainer.style.position = 'relative';

            try {
                // Create chart and verify data
                if (!store_data.data || !store_data.data.clusters) {
                    throw new Error('Invalid visualization data: missing clusters');
                }

                const chart = ClusterMap.ClusterMap().config(config);
                const selection = d3.select('#' + containerId)
                    .datum(store_data.data)
                    .call(chart);

                // Verify the selection was created successfully
                if (!selection || selection.empty()) {
                    throw new Error('Failed to create D3 selection for visualization');
                }

                // Center the visualization
                setTimeout(function() {
                    try {
                        const svg = selection.select('svg');
                        if (!svg.empty()) {
                            const svgNode = svg.node();
                            if (svgNode) {
                                // Get the SVG dimensions and center it
                                const bbox = svgNode.getBoundingClientRect();
                                const currentContainer = document.getElementById(containerId);
                                if (currentContainer) {
                                    const containerRect = currentContainer.getBoundingClientRect();

                                    if (bbox.width > 0 && bbox.height > 0) {
                                        // Center horizontally and vertically
                                        const leftOffset = Math.max(0, (containerRect.width - bbox.width) / 2);
                                        const topOffset = Math.max(0, (containerRect.height - bbox.height) / 2);

                                        svg.style('margin-left', leftOffset + 'px');
                                        svg.style('margin-top', topOffset + 'px');
                                    }
                                }
                            }
                        }
                    } catch (centerError) {
                        console.warn('Error centering visualization:', centerError);
                        // Don't fail the entire visualization for centering errors
                    }
                }, 100);

                console.log("ClusterMap created successfully");

                // Show and enable controls after successful rendering
                const controls = document.getElementById(controlsId);
                if (controls) controls.style.display = 'flex';

                const saveBtn = document.getElementById(saveBtnId);
                if (saveBtn) saveBtn.disabled = false;

                const resetBtn = document.getElementById(resetBtnId);
                if (resetBtn) resetBtn.disabled = false;

            } catch (error) {
                console.error("Error creating ClusterMap:", error);
                const errorRenderContainer = document.getElementById(containerId);
                if (errorRenderContainer) {
                    errorRenderContainer.innerHTML = '<div style="color: red; padding: 20px; text-align: center; background: #fff5f5; border: 1px solid #feb2b2; border-radius: 4px;"><strong>Visualization Error:</strong><br>' + error.message + '<br><br>Please try again with a different combination.</div>';
                    errorRenderContainer.style.display = 'flex';
                    errorRenderContainer.style.alignItems = 'center';
                    errorRenderContainer.style.justifyContent = 'center';
                    errorRenderContainer.style.backgroundColor = 'transparent';
                }

                // Hide controls on error
                const controls = document.getElementById(controlsId);
                if (controls) controls.style.display = 'none';
            }

            return ["", {"display": "none"}, false, false];

        } catch (error) {
            console.error("Error in synteny visualization:", error);
            return [`Error: ${error.message}`, {"display": "none"}, true, true];
        }
    }
    """,
    [
        Output("synteny-viz-message", "children", allow_duplicate=True),
        Output("synteny-controls", "style", allow_duplicate=True),
        Output("synteny-save-svg", "disabled", allow_duplicate=True),
        Output("synteny-reset-view", "disabled", allow_duplicate=True)
    ],
    [
        Input("synteny-data-store", "data"),
        Input("synteny-loading-state", "data"),
        Input("synteny-scale-factor", "value"),
        Input("synteny-cluster-spacing", "value"),
        Input("synteny-show-links", "value"),
        Input("synteny-show-gene-labels", "value"),
        Input("synteny-identity-threshold", "value"),
        Input("synteny-clear-viz", "n_clicks"),
    ],
    allow_duplicate=True,
    prevent_initial_call='initial_duplicate',
)

# Callback to trigger loading state when button is clicked
@callback(
    [
        Output("synteny-loading-state", "data", allow_duplicate=True),
        Output("synteny-update-button", "disabled", allow_duplicate=True),
        Output("synteny-update-button", "children", allow_duplicate=True),
    ],
    Input("synteny-update-button", "n_clicks"),
    prevent_initial_call=True
)
def trigger_loading_state(n_clicks):
    """Set loading state to True when processing starts"""
    if n_clicks:
        return True, True, [dmc.Loader(size="sm"), dmc.Text(" Processing...", size="lg", ml="xs")]
    return False, False, dmc.Text("Generate Visualization", size="lg")


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
                // Find the SVG in the visualization container
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
                // Find the SVG in the visualization container
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

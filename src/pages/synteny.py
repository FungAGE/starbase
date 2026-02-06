import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update, clientside_callback, ALL
from dash.dependencies import Output, Input, State
from dash_iconify import DashIconify
import json
import re
import tempfile
from pathlib import Path

from src.components.callbacks import create_modal_callback, handle_callback_error
# Note: search_components utilities available for future enhancements
# Current autocomplete uses Dash Mantine's built-in filtering
from src.config.logging import get_logger
from src.config.settings import GBK_PATH

logger = get_logger(__name__)

dash.register_page(__name__, path="/synteny")

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="synteny-url", refresh=False),
        dcc.Store(id="synteny-data-store"),
        dcc.Store(id="synteny-sequences-count", data=1),
        dcc.Store(id="synteny-available-ships", data=[]),
        
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title("Starship Synteny Viewer", order=1, mb="md"),
                dmc.Text(
                    "Interactive synteny visualization - Compare genomic sequences with detailed GFF annotation data",
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
                                # Dynamic Sequence Selectors
                                dmc.Stack([
                                    dmc.Group([
                                        dmc.Title("Select Sequences", order=2),
                                        dmc.Badge(
                                            id="sequence-count-badge",
                                            children="1 selected",
                                            color="blue",
                                            variant="filled",
                                        ),
                                    ], justify="space-between"),
                                    
                                dmc.Text(
                                    "Add sequences to compare (max 4)",
                                    size="sm",
                                    c="dimmed",
                                ),
                                
                                # Container for dynamic sequence selectors - all 4 created upfront
                                html.Div(
                                    id="sequence-selectors-container",
                                    children=[
                                        # Selector 0
                                        dmc.Paper(
                                            id="sequence-paper-0",
                                            children=dmc.Stack([
                                                dmc.Group([
                                                    dmc.Badge("Sequence 1", color="blue", variant="light"),
                                                    dmc.ActionIcon(
                                                        html.I(className="bi bi-x"),
                                                        id="remove-sequence-0",
                                                        variant="subtle",
                                                        color="red",
                                                        size="sm",
                                                        style={"visibility": "hidden"},
                                                    ),
                                                ], justify="space-between"),
                                                dmc.Autocomplete(
                                                    id="sequence-selector-0",
                                                    placeholder="Type to search taxonomy, family, accession...",
                                                    data=[],
                                                    limit=20,
                                                    leftSection=DashIconify(icon="bi:search"),
                                                    style={"width": "100%"},
                                                ),
                                            ], gap="xs"),
                                            p="sm",
                                            withBorder=True,
                                            radius="sm",
                                            mb="xs",
                                        ),
                                        # Selector 1
                                        dmc.Paper(
                                            id="sequence-paper-1",
                                            children=dmc.Stack([
                                                dmc.Group([
                                                    dmc.Badge("Sequence 2", color="green", variant="light"),
                                                    dmc.ActionIcon(
                                                        html.I(className="bi bi-x"),
                                                        id="remove-sequence-1",
                                                        variant="subtle",
                                                        color="red",
                                                        size="sm",
                                                    ),
                                                ], justify="space-between"),
                                                dmc.Autocomplete(
                                                    id="sequence-selector-1",
                                                    placeholder="Type to search taxonomy, family, accession...",
                                                    data=[],
                                                    limit=20,
                                                    leftSection=DashIconify(icon="bi:search"),
                                                    style={"width": "100%"},
                                                ),
                                            ], gap="xs"),
                                            p="sm",
                                            withBorder=True,
                                            radius="sm",
                                            mb="xs",
                                            style={"display": "none"},
                                        ),
                                        # Selector 2
                                        dmc.Paper(
                                            id="sequence-paper-2",
                                            children=dmc.Stack([
                                                dmc.Group([
                                                    dmc.Badge("Sequence 3", color="orange", variant="light"),
                                                    dmc.ActionIcon(
                                                        html.I(className="bi bi-x"),
                                                        id="remove-sequence-2",
                                                        variant="subtle",
                                                        color="red",
                                                        size="sm",
                                                    ),
                                                ], justify="space-between"),
                                                dmc.Autocomplete(
                                                    id="sequence-selector-2",
                                                    placeholder="Type to search taxonomy, family, accession...",
                                                    data=[],
                                                    limit=20,
                                                    leftSection=DashIconify(icon="bi:search"),
                                                    style={"width": "100%"},
                                                ),
                                            ], gap="xs"),
                                            p="sm",
                                            withBorder=True,
                                            radius="sm",
                                            mb="xs",
                                            style={"display": "none"},
                                        ),
                                        # Selector 3
                                        dmc.Paper(
                                            id="sequence-paper-3",
                                            children=dmc.Stack([
                                                dmc.Group([
                                                    dmc.Badge("Sequence 4", color="purple", variant="light"),
                                                    dmc.ActionIcon(
                                                        html.I(className="bi bi-x"),
                                                        id="remove-sequence-3",
                                                        variant="subtle",
                                                        color="red",
                                                        size="sm",
                                                    ),
                                                ], justify="space-between"),
                                                dmc.Autocomplete(
                                                    id="sequence-selector-3",
                                                    placeholder="Type to search taxonomy, family, accession...",
                                                    data=[],
                                                    limit=20,
                                                    leftSection=DashIconify(icon="bi:search"),
                                                    style={"width": "100%"},
                                                ),
                                            ], gap="xs"),
                                            p="sm",
                                            withBorder=True,
                                            radius="sm",
                                            mb="xs",
                                            style={"display": "none"},
                                        ),
                                    ],
                                ),
                                    
                                    dmc.Group([
                                        dmc.Button(
                                            "Add Sequence",
                                            id="add-sequence-button",
                                            variant="light",
                                            color="blue",
                                            size="sm",
                                            leftSection=html.I(className="bi bi-plus-circle"),
                                        ),
                                        dmc.Button(
                                            "Clear All",
                                            id="clear-sequences-button",
                                            variant="subtle",
                                            color="red",
                                            size="sm",
                                            leftSection=html.I(className="bi bi-trash"),
                                        ),
                                    ], gap="xs"),
                                ], gap="sm"),
                                
                                dmc.Divider(),
                                
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
                                                            description="Minimum similarity to show links",
                                                            id="synteny-identity-threshold",
                                                            value=0.3,
                                                            min=0.1,
                                                            max=1,
                                                            step=0.05,
                                                            decimalScale=2,
                                                        ),
                                                        dmc.NumberInput(
                                                            label="Scale Factor",
                                                            description="Zoom level for genes",
                                                            id="synteny-scale-factor",
                                                            value=10,
                                                            min=1,
                                                            max=50,
                                                        ),
                                                        dmc.NumberInput(
                                                            label="Cluster Spacing",
                                                            description="Vertical spacing between sequences",
                                                            id="synteny-cluster-spacing",
                                                            value=50,
                                                            min=10,
                                                            max=200,
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
                                                        dmc.Switch(
                                                            id="synteny-show-tooltips",
                                                            label="Enable gene tooltips",
                                                            checked=True,
                                                        ),
                                                    ], gap="sm")
                                                ])
                                            ],
                                        ),
                                        dmc.AccordionItem(
                                            value="display-options",
                                            children=[
                                                dmc.AccordionControl(dmc.Title("Display Options", order=3)),
                                                dmc.AccordionPanel([
                                                    dmc.Stack([
                                                        dmc.Select(
                                                            label="Color genes by",
                                                            description="Gene coloring scheme",
                                                            id="synteny-color-by",
                                                            value="category",
                                                            data=[
                                                                {"value": "category", "label": "Functional Category"},
                                                                {"value": "family", "label": "Gene Family (Target_ID)"},
                                                                {"value": "strand", "label": "Strand"},
                                                                {"value": "type", "label": "Feature Type"},
                                                                {"value": "source", "label": "Source"},
                                                            ],
                                                        ),
                                                        dmc.Select(
                                                            label="Gene label field",
                                                            description="Which attribute to show in labels",
                                                            id="synteny-label-field",
                                                            value="Alias",
                                                            data=[
                                                                {"value": "Alias", "label": "Gene Alias"},
                                                                {"value": "Target_ID", "label": "Target ID"},
                                                                {"value": "type", "label": "Feature Type"},
                                                            ],
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
                                        leftSection=html.I(className="bi bi-diagram-3"),
                                        fullWidth=True,
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
                
                # Right Column - Visualization and Details
                dmc.GridCol(
                    span={"base": 12, "md": 8},
                    children=[
                        # Visualization Container
                        dmc.Paper(
                            children=dmc.Stack([
                                html.Div(id="synteny-message"),
                                dcc.Loading(
                                    id="synteny-loading-viz",
                                    type="circle",
                                    children=[
                                        # STATIC CONTAINER
                                        html.Div(
                                            "Select sequences and click 'Generate Visualization'",
                                            id="synteny-static-viz",
                                            style={
                                                "height": "600px",
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
                                        html.Div(
                                            id="synteny-controls",
                                            style={"display": "none"},
                                            children=dmc.Group([
                                                dmc.Button(
                                                    "Save as SVG",
                                                    id="synteny-save-svg",
                                                    variant="light",
                                                    leftSection=html.I(className="bi bi-download"),
                                                    size="sm",
                                                ),
                                                dmc.Button(
                                                    "Reset View",
                                                    id="synteny-reset-view",
                                                    variant="light",
                                                    leftSection=html.I(className="bi bi-arrow-counterclockwise"),
                                                    size="sm",
                                                ),
                                                dmc.Button(
                                                    "Clear",
                                                    id="synteny-clear-viz",
                                                    variant="light",
                                                    color="red",
                                                    leftSection=html.I(className="bi bi-x-circle"),
                                                    size="sm",
                                                ),
                                            ], gap="xs", mt="sm"),
                                        ),
                                    ],
                                ),
                            ], gap="xs"),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            mb="md",
                        ),
                        
                        # Gene Details Panel
                        dmc.Paper(
                            children=dmc.Stack([
                                dmc.Title("Gene Details", order=3),
                                html.Div(
                                    id="gene-details-panel",
                                    children=dmc.Text(
                                        "Click on a gene in the visualization to see detailed information",
                                        c="dimmed",
                                        size="sm",
                                    ),
                                ),
                            ], gap="sm"),
                            p="md",
                            radius="md",
                            withBorder=True,
                        ),
                    ],
                ),
            ],
        ),
    ],
)


# Load available ships on page load (using correct DB access via synteny_queries)
@callback(
    Output("synteny-available-ships", "data"),
    Input("synteny-url", "pathname"),
)
@handle_callback_error
def load_available_ships(pathname):
    """Load list of available ships with GFF data from database using unified search."""
    if pathname != "/synteny":
        return no_update

    try:
        from src.utils.synteny_queries import get_ships_unified_search

        ships_data = get_ships_unified_search()
        
        # Format for Autocomplete: value = ship_id (for easy parsing), label = full searchable text
        options = []
        for ship in ships_data:
            # Value is just the ship ID (what gets returned when selected)
            value = str(ship['id'])
            
            # Label is the full searchable text (what user sees and searches)
            label = ship['label']
            
            options.append({
                "value": value,  # Simple ship ID
                "label": label,  # Full searchable display text
            })
        
        logger.info(f"Loaded {len(options)} ships with GFF data (unified search)")
        return options
    except Exception as e:
        logger.error(f"Error loading ships: {e}", exc_info=True)
        return []


# Callback to manage sequence selector visibility and populate autocomplete data
@callback(
    [Output(f"sequence-paper-{i}", "style") for i in range(4)] +
    [Output(f"sequence-selector-{i}", "data") for i in range(4)] +
    [Output(f"remove-sequence-{i}", "style") for i in range(4)] +
    [Output("sequence-count-badge", "children"),
     Output("synteny-sequences-count", "data")],
    [Input("add-sequence-button", "n_clicks"),
     Input("clear-sequences-button", "n_clicks")] +
     [Input(f"remove-sequence-{i}", "n_clicks") for i in range(4)],
    [State("synteny-sequences-count", "data"),
     State("synteny-available-ships", "data")],
    prevent_initial_call=False
)
@handle_callback_error
def manage_sequence_selectors(add_clicks, clear_clicks, *args):
    """Manage sequence selector visibility and populate their autocomplete data."""
    # Extract remove clicks and state from args
    remove_clicks = args[:4]  # First 4 args are remove button clicks
    current_count = args[4]    # 5th arg is current count
    ships_data = args[5]       # 6th arg is ships data
    
    ctx = dash.callback_context

    if not ctx.triggered:
        count = 1
    else:
        trigger = ctx.triggered[0]["prop_id"]
        if "add-sequence-button" in trigger and current_count < 4:
            count = current_count + 1
        elif "clear-sequences-button" in trigger:
            count = 1
        elif "remove-sequence" in trigger:
            count = max(1, current_count - 1)
        else:
            count = current_count

    if not isinstance(ships_data, list):
        ships_data = []

    # Generate outputs for all 4 selectors
    paper_styles = []
    selector_data = []
    remove_button_styles = []
    
    for i in range(4):
        is_visible = i < count
        
        # Paper visibility
        paper_styles.append({"display": "block" if is_visible else "none"})
        
        # Autocomplete data (populate all visible ones)
        selector_data.append(ships_data if is_visible else [])
        
        # Remove button visibility
        remove_button_styles.append({
            "visibility": "visible" if count > 1 and is_visible else "hidden"
        })

    badge_text = f"{count} selected" if count == 1 else f"{count} sequences"
    
    return paper_styles + selector_data + remove_button_styles + [badge_text, count]


# Main callback to generate visualization data
@callback(
    Output("synteny-data-store", "data"),
    Input("synteny-update-button", "n_clicks"),
    [State(f"sequence-selector-{i}", "value") for i in range(4)] +  # Max 4 sequences
    [State("synteny-color-by", "value"),
     State("synteny-label-field", "value")],
    prevent_initial_call=True
)
@handle_callback_error
def generate_synteny_data(n_clicks, *args):
    """
    Generate synteny visualization data from selected sequences.
    This queries the Gff table and formats data for ClusterMap.js
    """
    # Extract selected sequences and other states from args
    selected_sequences = args[:4]  # First 4 args are sequence values
    color_by = args[4]              # 5th arg is color_by
    label_field = args[5]           # 6th arg is label_field
    from src.utils.synteny_queries import get_gff_by_ship_ids  # Adjust import
    
    if not n_clicks or not any(selected_sequences):
        return no_update
    
    # Filter out None values and convert to integers
    # Autocomplete returns ship IDs as strings
    ship_ids = []
    for s in selected_sequences:
        if s is not None and str(s).strip():
            try:
                ship_ids.append(int(s))
            except ValueError:
                logger.warning(f"Could not parse ship ID from: {s}")
    
    if not ship_ids:
        return {"error": "Please select at least one sequence", "success": False}
    
    try:
        # Get GFF data for selected ships
        gff_data = get_gff_by_ship_ids(ship_ids)
        
        if not gff_data:
            return {"error": "No GFF data found for selected sequences", "success": False}
        
        # Transform data to ClusterMap.js format
        clusters = []
        for ship_id in ship_ids:
            ship_gff = [g for g in gff_data if g["ship_id"] == ship_id]
            
            if not ship_gff:
                continue
            
            # Parse attributes to extract key information
            genes = []
            for idx, gff_entry in enumerate(ship_gff):
                # Parse the attributes string
                attrs = parse_gff_attributes(gff_entry["attributes"])
                
                # Categorize the gene based on functional patterns
                category = categorize_gene(gff_entry)
                
                # Determine color grouping based on color_by option
                if color_by == "category":
                    group = category
                elif color_by == "family":
                    group = attrs.get("Target_ID", category)
                else:
                    group = gff_entry.get(color_by, "unknown")
                
                # Create gene object with full metadata
                gene = {
                    "uid": f"gene_{ship_id}_{idx}",
                    "start": gff_entry["start"],
                    "end": gff_entry["end"],
                    "strand": 1 if gff_entry["strand"] == "+" else -1,
                    "_start": gff_entry["start"],
                    "_end": gff_entry["end"],
                    "_strand": 1 if gff_entry["strand"] == "+" else -1,
                    
                    # Labels and identifiers
                    "name": attrs.get("Alias", attrs.get("Target_ID", f"gene_{idx}")),
                    "label": attrs.get(label_field, attrs.get("Alias", f"gene_{idx}")),
                    "locus_tag": attrs.get("SeqID", ""),
                    
                    # Full metadata for tooltips
                    "metadata": {
                        "id": gff_entry["id"],
                        "source": gff_entry["source"],
                        "type": gff_entry["type"],
                        "score": gff_entry["score"],
                        "phase": gff_entry["phase"],
                        "Target_ID": attrs.get("Target_ID", ""),
                        "Alias": attrs.get("Alias", ""),
                        "SeqID": attrs.get("SeqID", ""),
                        "category": category,
                        "category_label": GENE_CATEGORY_LABELS.get(category, category),
                        **attrs  # Include all other attributes
                    },
                    
                    # Color grouping
                    "group": group,
                    
                    # Set explicit color for category-based coloring
                    "colour": GENE_CATEGORY_COLORS.get(category) if color_by == "category" else None,
                }
                genes.append(gene)
            
            # Create locus (genomic region)
            locus = {
                "uid": f"locus_{ship_id}",
                "name": ship_gff[0].get("accession", f"Locus {ship_id}"),
                "start": min(g["start"] for g in genes),
                "end": max(g["end"] for g in genes),
                "_start": min(g["start"] for g in genes),
                "_end": max(g["end"] for g in genes),
                "genes": genes,
            }
            
            # Create cluster (collection of loci)
            cluster = {
                "uid": f"cluster_{ship_id}",
                "name": f"{ship_gff[0]['ship_accession_display']} - {ship_gff[0]['name']}",
                "loci": [locus],
            }
            clusters.append(cluster)
        
        # Generate links between homologous genes based on Target_ID
        links = generate_synteny_links(clusters)
        logger.info(f"Generated {len(clusters)} clusters with {len(links)} links")

        return {
            "clusters": clusters,
            "links": links,
            "success": True,
        }
    
    except Exception as e:
        logger.error(f"Error generating synteny data: {e}", exc_info=True)
        return {"error": f"Error: {str(e)}", "success": False}


def parse_gff_attributes(attr_string):
    """Parse GFF attributes string into dictionary."""
    if not attr_string:
        return {}
    attrs = {}
    for item in attr_string.split(";"):
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
    return attrs


def categorize_gene(gff_entry):
    """
    Categorize a gene based on functional patterns in its attributes.
    
    Categories:
    - captain: Tyrosine recombinase or capsid proteins (tyr|cap)
    - nlr: NLR genes
    - plp: PLP genes  
    - fre: FRE genes
    - other: Genes not matching any specific pattern
    
    Args:
        gff_entry: Dictionary with GFF fields including attributes
        
    Returns:
        str: Category name ('captain', 'nlr', 'plp', 'fre', or 'other')
    """
    attrs = parse_gff_attributes(gff_entry.get("attributes", ""))
    
    # Collect all searchable text
    searchable_text = " ".join([
        attrs.get("Alias", ""),
        attrs.get("Target_ID", ""),
        attrs.get("product", ""),
        attrs.get("Name", ""),
        attrs.get("gene", ""),
        gff_entry.get("type", ""),
        gff_entry.get("source", ""),
    ]).lower()
    
    # Pattern matching for categories (order matters - most specific first)
    if re.search(r"\b(tyr|cap|tyrosine|capsid|recombinase)\b", searchable_text):
        return "captain"
    elif re.search(r"\bnlr\b", searchable_text):
        return "nlr"
    elif re.search(r"\bplp\b", searchable_text):
        return "plp"
    elif re.search(r"\bfre\b", searchable_text):
        return "fre"
    else:
        return "other"


# Define color scheme for gene categories
GENE_CATEGORY_COLORS = {
    "captain": "#e74c3c",  # Red
    "nlr": "#3498db",      # Blue
    "plp": "#2ecc71",      # Green
    "fre": "#f39c12",      # Orange
    "other": "#95a5a6",    # Gray
}

GENE_CATEGORY_LABELS = {
    "captain": "Captain (Tyr/Cap)",
    "nlr": "NLR",
    "plp": "PLP",
    "fre": "FRE",
    "other": "Other",
}


def generate_synteny_links(clusters):
    """
    Generate synteny links between genes based on Target_IDs and functional categories.
    Links are created when genes share the same Target_ID or functional category.
    """
    links = []
    link_id = 0
    
    # Build index of genes by Target_ID and category
    target_index = {}
    category_index = {}
    
    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                # Index by Target_ID
                target_id = gene["metadata"].get("Target_ID")
                if target_id and target_id != "":
                    if target_id not in target_index:
                        target_index[target_id] = []
                    target_index[target_id].append(gene)
                
                # Index by category
                category = gene["metadata"].get("category", "other")
                if category not in category_index:
                    category_index[category] = []
                category_index[category].append(gene)
    
    # Create links for genes with same Target_ID (highest confidence)
    for target_id, genes in target_index.items():
        if len(genes) > 1:
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    gene_i_cat = genes[i]["metadata"].get("category", "other")
                    gene_j_cat = genes[j]["metadata"].get("category", "other")
                    
                    # Same Target_ID gets high identity
                    # Matching category boosts identity slightly
                    identity = 0.95 if gene_i_cat == gene_j_cat else 0.85
                    
                    link = {
                        "uid": f"link_{link_id}",
                        "query": {"uid": genes[i]["uid"]},
                        "target": {"uid": genes[j]["uid"]},
                        "identity": identity,
                    }
                    links.append(link)
                    link_id += 1
    
    # Create additional links for genes in same functional category but different Target_IDs
    # Only for important categories (captain, nlr, plp, fre) and only if they don't already have Target_ID links
    linked_genes = set()
    for link in links:
        linked_genes.add(link["query"]["uid"])
        linked_genes.add(link["target"]["uid"])
    
    for category in ["captain", "nlr", "plp", "fre"]:
        if category in category_index:
            genes = [g for g in category_index[category] if g["uid"] not in linked_genes]
            if len(genes) > 1:
                # Create links between genes in same category
                for i in range(len(genes)):
                    for j in range(i + 1, len(genes)):
                        link = {
                            "uid": f"link_{link_id}",
                            "query": {"uid": genes[i]["uid"]},
                            "target": {"uid": genes[j]["uid"]},
                            "identity": 0.6,  # Lower identity for category-only matches
                        }
                        links.append(link)
                        link_id += 1
    
    return links


# Enhanced clientside callback for visualization with tooltip support
clientside_callback(
    """
    function(store_data, scale_factor, cluster_spacing, show_links, show_gene_labels, 
             identity_threshold, show_tooltips, clear_clicks) {
        
        const containerId = 'synteny-static-viz';
        
        if (clear_clicks) {
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = 'Select sequences and click "Generate Visualization"';
                container.style.display = 'flex';
                container.style.alignItems = 'center';
                container.style.justifyContent = 'center';
                container.style.color = '#6c757d';
                container.style.fontSize = '16px';
            }
            return ["", {"display": "none"}];
        }
        
        if (!store_data || !store_data.success) {
            const container = document.getElementById(containerId);
            if (container && store_data && store_data.error) {
                container.innerHTML = '<div style="color: red; padding: 20px; text-align: center;">' +
                    store_data.error + '</div>';
            }
            return ["", {"display": "none"}];
        }
        
        try {
            const config = {
                plot: {
                    scaleFactor: scale_factor || 10,
                    scaleGenes: true,
                },
                cluster: {
                    spacing: cluster_spacing || 50,
                    alignLabels: true,
                },
                gene: {
                    label: {
                        show: show_gene_labels || false,
                    },
                    shape: {
                        onClick: function(event, gene) {
                            if (gene && gene.metadata) {
                                const panel = document.getElementById('gene-details-panel');
                                if (panel) {
                                    let html = '<div style="font-family: system-ui; font-size: 13px;">';
                                    
                                    // Header with gene name
                                    html += '<div style="margin-bottom: 16px;">';
                                    html += '<div style="font-size: 18px; font-weight: bold; color: #228be6; margin-bottom: 4px;">';
                                    html += (gene.label || gene.name || 'Unknown Gene');
                                    html += '</div>';
                                    
                                    // Add category badge
                                    if (gene.metadata.category_label) {
                                        const categoryColors = {
                                            'captain': '#e74c3c',
                                            'nlr': '#3498db',
                                            'plp': '#2ecc71',
                                            'fre': '#f39c12',
                                            'other': '#95a5a6'
                                        };
                                        const bgColor = categoryColors[gene.metadata.category] || '#95a5a6';
                                        html += '<span style="display: inline-block; background: ' + bgColor + '; color: white; ';
                                        html += 'padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: 600;">';
                                        html += gene.metadata.category_label + '</span>';
                                    }
                                    html += '</div>';
                                    
                                    // Genomic Location
                                    html += '<div style="background: #f8f9fa; padding: 12px; border-radius: 6px; margin-bottom: 12px;">';
                                    html += '<div style="font-weight: bold; margin-bottom: 8px; color: #495057;">📍 Genomic Location</div>';
                                    html += '<div style="font-size: 12px;">';
                                    html += '<div style="margin-bottom: 4px;"><strong>Position:</strong> ' +
                                        gene.start.toLocaleString() + ' - ' + gene.end.toLocaleString() + '</div>';
                                    html += '<div style="margin-bottom: 4px;"><strong>Length:</strong> ' +
                                        (gene.end - gene.start).toLocaleString() + ' bp</div>';
                                    html += '<div><strong>Strand:</strong> ' + (gene.strand > 0 ? 'Forward (+)' : 'Reverse (-)') + '</div>';
                                    html += '</div></div>';
                                    
                                    // Primary Annotation
                                    if (gene.metadata) {
                                        html += '<div style="background: #f8f9fa; padding: 12px; border-radius: 6px; margin-bottom: 12px;">';
                                        html += '<div style="font-weight: bold; margin-bottom: 8px; color: #495057;">🏷️ Primary Annotation</div>';
                                        html += '<div style="font-size: 12px;">';
                                        
                                        const primaryFields = ['Target_ID', 'Alias', 'SeqID', 'product', 'Name', 'type', 'source'];
                                        for (const key of primaryFields) {
                                            if (gene.metadata[key] && gene.metadata[key] !== '.' && gene.metadata[key] !== '') {
                                                html += '<div style="margin-bottom: 4px;"><strong>' + key + ':</strong> ';
                                                html += '<span style="font-family: monospace; background: #e9ecef; padding: 1px 4px; border-radius: 3px;">';
                                                html += gene.metadata[key] + '</span></div>';
                                            }
                                        }
                                        html += '</div></div>';
                                        
                                        // Additional Attributes
                                        const excludedKeys = ['id', 'category', 'category_label', ...primaryFields];
                                        const additionalAttrs = Object.keys(gene.metadata).filter(k => 
                                            !excludedKeys.includes(k) && 
                                            gene.metadata[k] && 
                                            gene.metadata[k] !== '.' && 
                                            gene.metadata[k] !== ''
                                        );
                                        
                                        if (additionalAttrs.length > 0) {
                                            html += '<div style="background: #f8f9fa; padding: 12px; border-radius: 6px;">';
                                            html += '<div style="font-weight: bold; margin-bottom: 8px; color: #495057;">📋 Additional Attributes</div>';
                                            html += '<div style="font-size: 12px;">';
                                            for (const key of additionalAttrs) {
                                                html += '<div style="margin-bottom: 4px;"><strong>' + key + ':</strong> ' + gene.metadata[key] + '</div>';
                                            }
                                            html += '</div></div>';
                                        }
                                    }
                                    html += '</div>';
                                    panel.innerHTML = html;
                                }
                            }
                        }
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
                return ["Error: Container not found", {"display": "none"}];
            }
            
            // Clear container
            container.innerHTML = '';
            container.style.display = 'block';
            container.style.alignItems = 'unset';
            container.style.justifyContent = 'unset';
            
            if (typeof d3 !== 'undefined') {
                d3.selectAll('#' + containerId + ' *').remove();
            }
            
            // Create visualization
            setTimeout(function() {
                try {
                    const chart = ClusterMap.ClusterMap().config(config);
                    d3.select('#' + containerId)
                        .datum(store_data)
                        .call(chart);
                    
                    if (show_tooltips && typeof d3 !== 'undefined') {
                        const tooltip = d3.select('body').selectAll('.gene-tooltip').data([null])
                            .join('div')
                            .attr('class', 'gene-tooltip')
                            .style('position', 'absolute')
                            .style('background', 'rgba(255, 255, 255, 0.98)')
                            .style('color', '#212529')
                            .style('padding', '12px')
                            .style('border-radius', '6px')
                            .style('box-shadow', '0 4px 12px rgba(0,0,0,0.15)')
                            .style('border', '1px solid #dee2e6')
                            .style('pointer-events', 'none')
                            .style('opacity', 0)
                            .style('z-index', 10000)
                            .style('max-width', '350px')
                            .style('font-family', 'system-ui, -apple-system, sans-serif');
                        d3.selectAll('#' + containerId + ' .genePolygon')
                            .on('mouseenter', function(event, d) {
                                if (d && d.metadata) {
                                    let html = '<div style="font-size: 12px; line-height: 1.5;">';
                                    
                                    // Gene name
                                    html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 6px; color: #228be6;">';
                                    html += (d.label || d.name || 'Unknown Gene') + '</div>';
                                    
                                    // Category badge
                                    if (d.metadata.category_label) {
                                        const categoryColors = {
                                            'captain': '#e74c3c',
                                            'nlr': '#3498db',
                                            'plp': '#2ecc71',
                                            'fre': '#f39c12',
                                            'other': '#95a5a6'
                                        };
                                        const bgColor = categoryColors[d.metadata.category] || '#95a5a6';
                                        html += '<div style="margin-bottom: 6px;">';
                                        html += '<span style="background: ' + bgColor + '; color: white; padding: 2px 6px; ';
                                        html += 'border-radius: 10px; font-size: 10px; font-weight: 600;">';
                                        html += d.metadata.category_label + '</span></div>';
                                    }
                                    
                                    // Position and basic info
                                    html += '<div style="border-top: 1px solid #dee2e6; padding-top: 6px; margin-top: 6px;">';
                                    html += '<strong>Position:</strong> ' + d.start.toLocaleString() + ' - ' + d.end.toLocaleString();
                                    html += ' (' + (d.end - d.start).toLocaleString() + ' bp)<br/>';
                                    html += '<strong>Strand:</strong> ' + (d.strand > 0 ? 'Forward (+)' : 'Reverse (-)')  + '<br/>';
                                    
                                    // Key metadata fields
                                    if (d.metadata.Target_ID && d.metadata.Target_ID !== '') {
                                        html += '<strong>Target ID:</strong> ' + d.metadata.Target_ID + '<br/>';
                                    }
                                    if (d.metadata.Alias && d.metadata.Alias !== '') {
                                        html += '<strong>Alias:</strong> ' + d.metadata.Alias + '<br/>';
                                    }
                                    if (d.metadata.product && d.metadata.product !== '') {
                                        html += '<strong>Product:</strong> ' + d.metadata.product + '<br/>';
                                    }
                                    if (d.metadata.type && d.metadata.type !== '') {
                                        html += '<strong>Type:</strong> ' + d.metadata.type + '<br/>';
                                    }
                                    html += '</div>';
                                    html += '<div style="margin-top: 6px; font-size: 10px; color: #6c757d; font-style: italic;">Click for full details</div>';
                                    html += '</div>';
                                    
                                    tooltip.html(html)
                                        .style('left', (event.pageX + 15) + 'px')
                                        .style('top', (event.pageY - 10) + 'px')
                                        .style('opacity', 1);
                                }
                            })
                            .on('mousemove', function(event) {
                                tooltip
                                    .style('left', (event.pageX + 15) + 'px')
                                    .style('top', (event.pageY - 10) + 'px');
                            })
                            .on('mouseleave', function() {
                                tooltip.style('opacity', 0);
                            });
                    }
                    
                } catch (error) {
                    console.error("Error creating ClusterMap:", error);
                    container.innerHTML = '<div style="color: red; padding: 20px;">Error: ' + error.message + '</div>';
                }
            }, 200);
            
            return ["", {"display": "block"}];
            
        } catch (error) {
            console.error("Error in synteny visualization:", error);
            return ['Error: ' + error.message, {"display": "none"}];
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
        Input("synteny-show-tooltips", "checked"),
        Input("synteny-clear-viz", "n_clicks"),
    ]
)


# SVG export callback
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) {
            try {
                const svg = document.querySelector('#synteny-static-viz svg');
                if (svg) {
                    const serializer = new XMLSerializer();
                    const svgString = serializer.serializeToString(svg);
                    const blob = new Blob([svgString], {type: 'image/svg+xml'});
                    
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement('a');
                    link.href = url;
                    link.download = 'starship-synteny-' + Date.now() + '.svg';
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                    URL.revokeObjectURL(url);
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


# Reset view callback
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) {
            try {
                const svg = document.querySelector('#synteny-static-viz svg');
                if (svg) {
                    const g = svg.querySelector('g');
                    if (g) {
                        g.setAttribute('transform', 'translate(20, 50) scale(1.2)');
                    }
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


# Clear visualization callback
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) {
            return null;
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("synteny-data-store", "data", allow_duplicate=True),
    Input("synteny-clear-viz", "n_clicks"),
    prevent_initial_call=True
)

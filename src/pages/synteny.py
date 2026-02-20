import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback, no_update, clientside_callback, ALL
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
from dash_iconify import DashIconify
import json
import re

from src.components.callbacks import handle_callback_error
from src.config.logging import get_logger

logger = get_logger(__name__)

dash.register_page(__name__, path="/synteny")

SEQUENCE_COLORS = ["blue", "green", "orange", "grape"]

GENE_CATEGORY_COLORS = {
    "captain": "#e74c3c",
    "nlr": "#3498db",
    "plp": "#2ecc71",
    "fre": "#f39c12",
    "other": "#95a5a6",
}

GENE_CATEGORY_LABELS = {
    "captain": "Captain (Tyr)",
    "nlr": "NLR",
    "plp": "PLP",
    "fre": "FRE",
    "other": "Other",
}

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="synteny-url", refresh=False),
        dcc.Store(id="synteny-data-store"),
        dcc.Store(id="synteny-available-ships", data=[]),
        dcc.Store(id="synteny-filtered-ships", data=None),
        dcc.Store(id="synteny-selected-ships", data=[]),
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title("Starship Synteny Viewer", order=1, mb="md"),
                dmc.Text(
                    "Interactive synteny visualization — compare genomic sequences with detailed GFF annotation data",
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
                # Left Column — Search, Results, Selection, Settings
                dmc.GridCol(
                    span={"base": 12, "md": 4},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack(
                                [
                                    # ── Search Section ────────────────────────
                                    dmc.Stack(
                                        [
                                            dmc.Title("Search Sequences", order=2),
                                            dmc.Text(
                                                "Filter available sequences by taxonomy, family, or SSB accession",
                                                size="sm",
                                                c="dimmed",
                                            ),
                                            dmc.Autocomplete(
                                                id="synteny-taxa-search",
                                                label="Taxonomy",
                                                placeholder="Search species, genus, etc.",
                                                data=[],
                                                limit=20,
                                                style={"width": "100%"},
                                            ),
                                            dmc.Autocomplete(
                                                id="synteny-family-search",
                                                label="Starship Family",
                                                placeholder="Search by family...",
                                                data=[],
                                                limit=20,
                                                style={"width": "100%"},
                                            ),
                                            dmc.Autocomplete(
                                                id="synteny-ssb-search",
                                                label="SSB Accession",
                                                placeholder="Search by SSB accession...",
                                                data=[],
                                                limit=20,
                                                style={"width": "100%"},
                                            ),
                                            dmc.Group(
                                                [
                                                    dmc.Button(
                                                        "Reset",
                                                        id="synteny-reset-search",
                                                        variant="outline",
                                                        leftSection=DashIconify(
                                                            icon="tabler:refresh"
                                                        ),
                                                        size="sm",
                                                    ),
                                                    dmc.Button(
                                                        "Search",
                                                        id="synteny-apply-search",
                                                        variant="filled",
                                                        leftSection=DashIconify(
                                                            icon="tabler:search"
                                                        ),
                                                        size="sm",
                                                    ),
                                                ],
                                                justify="flex-end",
                                            ),
                                        ],
                                        gap="sm",
                                    ),
                                    dmc.Divider(),
                                    # ── Results Section ───────────────────────
                                    dmc.Stack(
                                        [
                                            dmc.Title("Results", order=3),
                                            dcc.Loading(
                                                type="circle",
                                                children=html.Div(
                                                    id="synteny-search-results",
                                                    children=dmc.Text(
                                                        "Use the search fields above to find sequences",
                                                        c="dimmed",
                                                        size="sm",
                                                    ),
                                                ),
                                            ),
                                        ],
                                        gap="xs",
                                    ),
                                    dmc.Divider(),
                                    # ── Selected Sequences Section ────────────
                                    dmc.Stack(
                                        [
                                            dmc.Group(
                                                [
                                                    dmc.Title(
                                                        "Selected for Comparison",
                                                        order=3,
                                                    ),
                                                    dmc.Badge(
                                                        id="synteny-selected-count",
                                                        children="0 / 4",
                                                        color="blue",
                                                        variant="filled",
                                                    ),
                                                ],
                                                justify="space-between",
                                            ),
                                            html.Div(
                                                id="synteny-selected-display",
                                                children=dmc.Text(
                                                    "No sequences selected yet",
                                                    c="dimmed",
                                                    size="sm",
                                                ),
                                            ),
                                            dmc.Button(
                                                "Clear All",
                                                id="synteny-clear-selected",
                                                variant="subtle",
                                                color="red",
                                                size="sm",
                                                leftSection=DashIconify(
                                                    icon="tabler:trash"
                                                ),
                                            ),
                                        ],
                                        gap="sm",
                                    ),
                                    dmc.Divider(),
                                    # ── Visualization Settings ────────────────
                                    dmc.Accordion(
                                        variant="filled",
                                        chevronPosition="right",
                                        chevronSize=16,
                                        value="visualization-settings",
                                        children=[
                                            dmc.AccordionItem(
                                                value="visualization-settings",
                                                children=[
                                                    dmc.AccordionControl(
                                                        dmc.Title(
                                                            "Visualization Settings",
                                                            order=3,
                                                        )
                                                    ),
                                                    dmc.AccordionPanel(
                                                        [
                                                            dmc.Stack(
                                                                [
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
                                                                ],
                                                                gap="sm",
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                            ),
                                            dmc.AccordionItem(
                                                value="display-options",
                                                children=[
                                                    dmc.AccordionControl(
                                                        dmc.Title(
                                                            "Display Options", order=3
                                                        )
                                                    ),
                                                    dmc.AccordionPanel(
                                                        [
                                                            dmc.Stack(
                                                                [
                                                                    dmc.Select(
                                                                        label="Color genes by",
                                                                        description="Gene coloring scheme",
                                                                        id="synteny-color-by",
                                                                        value="category",
                                                                        data=[
                                                                            {
                                                                                "value": "category",
                                                                                "label": "Functional Category",
                                                                            },
                                                                            {
                                                                                "value": "family",
                                                                                "label": "Gene Family (Target_ID)",
                                                                            },
                                                                            {
                                                                                "value": "strand",
                                                                                "label": "Strand",
                                                                            },
                                                                            {
                                                                                "value": "type",
                                                                                "label": "Feature Type",
                                                                            },
                                                                            {
                                                                                "value": "source",
                                                                                "label": "Source",
                                                                            },
                                                                        ],
                                                                    ),
                                                                    dmc.Select(
                                                                        label="Gene label field",
                                                                        description="Which attribute to show in labels",
                                                                        id="synteny-label-field",
                                                                        value="Alias",
                                                                        data=[
                                                                            {
                                                                                "value": "Alias",
                                                                                "label": "Gene Alias",
                                                                            },
                                                                            {
                                                                                "value": "Target_ID",
                                                                                "label": "Target ID",
                                                                            },
                                                                            {
                                                                                "value": "type",
                                                                                "label": "Feature Type",
                                                                            },
                                                                        ],
                                                                    ),
                                                                ],
                                                                gap="sm",
                                                            ),
                                                        ]
                                                    ),
                                                ],
                                            ),
                                        ],
                                    ),
                                    dmc.Button(
                                        dmc.Text("Generate Visualization", size="lg"),
                                        id="synteny-update-button",
                                        variant="gradient",
                                        gradient={"from": "indigo", "to": "cyan"},
                                        leftSection=html.I(className="bi bi-diagram-3"),
                                        fullWidth=True,
                                    ),
                                ],
                                gap="md",
                            ),
                            p="xl",
                            radius="md",
                            withBorder=True,
                        ),
                    ],
                ),
                # Right Column — Visualization and Gene Details
                dmc.GridCol(
                    span={"base": 12, "md": 8},
                    children=[
                        dmc.Paper(
                            children=dmc.Stack(
                                [
                                    html.Div(id="synteny-message"),
                                    dcc.Loading(
                                        id="synteny-loading-viz",
                                        type="circle",
                                        children=[
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
                                            html.Div(
                                                id="synteny-controls",
                                                style={"display": "none"},
                                                children=dmc.Group(
                                                    [
                                                        dmc.Button(
                                                            "Save as SVG",
                                                            id="synteny-save-svg",
                                                            variant="light",
                                                            leftSection=html.I(
                                                                className="bi bi-download"
                                                            ),
                                                            size="sm",
                                                        ),
                                                        dmc.Button(
                                                            "Reset View",
                                                            id="synteny-reset-view",
                                                            variant="light",
                                                            leftSection=html.I(
                                                                className="bi bi-arrow-counterclockwise"
                                                            ),
                                                            size="sm",
                                                        ),
                                                        dmc.Button(
                                                            "Clear",
                                                            id="synteny-clear-viz",
                                                            variant="light",
                                                            color="red",
                                                            leftSection=html.I(
                                                                className="bi bi-x-circle"
                                                            ),
                                                            size="sm",
                                                        ),
                                                    ],
                                                    gap="xs",
                                                    mt="sm",
                                                ),
                                            ),
                                        ],
                                    ),
                                    dmc.Divider(mt="xs", mb="xs"),
                                    dmc.Group(
                                        [
                                            dmc.Group(
                                                [
                                                    html.Div(
                                                        style={
                                                            "width": "14px",
                                                            "height": "14px",
                                                            "borderRadius": "3px",
                                                            "backgroundColor": GENE_CATEGORY_COLORS[
                                                                cat
                                                            ],
                                                            "flexShrink": "0",
                                                        }
                                                    ),
                                                    dmc.Text(
                                                        GENE_CATEGORY_LABELS[cat],
                                                        size="xs",
                                                    ),
                                                ],
                                                gap="xs",
                                            )
                                            for cat in GENE_CATEGORY_COLORS
                                        ],
                                        gap="md",
                                    ),
                                ],
                                gap="xs",
                            ),
                            p="xl",
                            radius="md",
                            withBorder=True,
                            mb="md",
                        ),
                        dmc.Paper(
                            children=dmc.Stack(
                                [
                                    dmc.Title("Gene Details", order=3),
                                    html.Div(
                                        id="gene-details-panel",
                                        children=dmc.Text(
                                            "Click on a gene in the visualization to see detailed information",
                                            c="dimmed",
                                            size="sm",
                                        ),
                                    ),
                                ],
                                gap="sm",
                            ),
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


# ── Load all available ships on page load ──────────────────────────────────────
@callback(
    Output("synteny-available-ships", "data"),
    Input("synteny-url", "pathname"),
)
@handle_callback_error
def load_available_ships(pathname):
    if pathname != "/synteny":
        return no_update

    try:
        from src.utils.synteny_queries import get_ships_unified_search

        ships_data = get_ships_unified_search()

        options = []
        for ship in ships_data:
            options.append(
                {
                    "value": str(ship["id"]),
                    "label": ship["label"],
                    "taxonomy_name": ship.get("taxonomy_name", ""),
                    "familyName": ship.get("familyName", ""),
                    "ship_accession_display": ship.get("ship_accession_display", ""),
                }
            )

        logger.info(f"Loaded {len(options)} ships with GFF data")
        return options
    except Exception as e:
        logger.error(f"Error loading ships: {e}", exc_info=True)
        return []


# ── Populate autocomplete options from available ships ─────────────────────────
@callback(
    [
        Output("synteny-taxa-search", "data"),
        Output("synteny-family-search", "data"),
        Output("synteny-ssb-search", "data"),
    ],
    Input("synteny-available-ships", "data"),
)
@handle_callback_error
def populate_search_options(ships_data):
    if not ships_data:
        return [], [], []

    taxa_values = set()
    family_values = set()
    ssb_values = set()

    for ship in ships_data:
        if ship.get("taxonomy_name"):
            taxa_values.add(ship["taxonomy_name"])
        if ship.get("familyName"):
            family_values.add(ship["familyName"])
        if ship.get("ship_accession_display"):
            ssb_values.add(ship["ship_accession_display"])

    def _to_opts(values):
        return sorted(
            [{"value": v, "label": v} for v in values],
            key=lambda x: len(x["label"]),
        )

    return _to_opts(taxa_values), _to_opts(family_values), _to_opts(ssb_values)


# ── Search / Reset ─────────────────────────────────────────────────────────────
@callback(
    [
        Output("synteny-filtered-ships", "data"),
        Output("synteny-taxa-search", "value"),
        Output("synteny-family-search", "value"),
        Output("synteny-ssb-search", "value"),
    ],
    [
        Input("synteny-apply-search", "n_clicks"),
        Input("synteny-reset-search", "n_clicks"),
    ],
    [
        State("synteny-taxa-search", "value"),
        State("synteny-family-search", "value"),
        State("synteny-ssb-search", "value"),
        State("synteny-available-ships", "data"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def search_synteny_ships(
    search_clicks, reset_clicks, taxa_val, family_val, ssb_val, ships_data
):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if triggered_id == "synteny-reset-search":
        return None, "", "", ""

    if not search_clicks or not ships_data:
        raise PreventUpdate

    filtered = list(ships_data)

    if taxa_val and taxa_val.strip():
        filtered = [
            s
            for s in filtered
            if taxa_val.lower() in s.get("taxonomy_name", "").lower()
        ]

    if family_val and family_val.strip():
        filtered = [
            s for s in filtered if family_val.lower() == s.get("familyName", "").lower()
        ]

    if ssb_val and ssb_val.strip():
        filtered = [
            s
            for s in filtered
            if ssb_val.lower() in s.get("ship_accession_display", "").lower()
        ]

    return filtered, no_update, no_update, no_update


# ── Display search results ─────────────────────────────────────────────────────
@callback(
    Output("synteny-search-results", "children"),
    [
        Input("synteny-filtered-ships", "data"),
        Input("synteny-selected-ships", "data"),
    ],
)
@handle_callback_error
def display_search_results(filtered_ships, selected_ships):
    if filtered_ships is None:
        return dmc.Text(
            "Use the search fields above to find sequences",
            c="dimmed",
            size="sm",
        )

    if not filtered_ships:
        return dmc.Alert(
            "No sequences found matching your search criteria.",
            color="blue",
            variant="light",
        )

    selected_ids = {str(s["id"]) for s in (selected_ships or [])}
    at_capacity = len(selected_ids) >= 4

    items = []
    for ship in filtered_ships:
        ship_id = str(ship["value"])
        label = ship.get("label", ship_id)
        already_added = ship_id in selected_ids
        disabled = already_added or (at_capacity and not already_added)

        items.append(
            dmc.Paper(
                dmc.Group(
                    [
                        dmc.Text(
                            label,
                            size="xs",
                            style={"flex": 1, "minWidth": 0, "wordBreak": "break-word"},
                        ),
                        dmc.ActionIcon(
                            DashIconify(
                                icon="tabler:circle-check"
                                if already_added
                                else "tabler:circle-plus"
                            ),
                            id={"type": "synteny-add-ship", "index": ship_id},
                            variant="subtle",
                            color="green" if already_added else "blue",
                            disabled=disabled,
                            size="sm",
                        ),
                    ],
                    justify="space-between",
                    gap="xs",
                    wrap="nowrap",
                ),
                p="xs",
                withBorder=True,
                radius="sm",
                style={"opacity": "0.55" if already_added else "1"},
            )
        )

    return dmc.Stack(
        [
            dmc.Text(
                f"{len(filtered_ships)} sequence{'s' if len(filtered_ships) != 1 else ''} found",
                size="xs",
                c="dimmed",
            ),
            html.Div(
                items,
                style={"maxHeight": "280px", "overflowY": "auto"},
            ),
        ],
        gap="xs",
    )


# ── Add / remove / clear selected ships ───────────────────────────────────────
@callback(
    Output("synteny-selected-ships", "data"),
    [
        Input({"type": "synteny-add-ship", "index": ALL}, "n_clicks"),
        Input({"type": "synteny-remove-ship", "index": ALL}, "n_clicks"),
        Input("synteny-clear-selected", "n_clicks"),
    ],
    [
        State("synteny-selected-ships", "data"),
        State("synteny-available-ships", "data"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def manage_selected_ships(
    add_clicks, remove_clicks, clear_clicks, current_selected, all_ships
):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    triggered_prop = ctx.triggered[0]["prop_id"]
    triggered_value = ctx.triggered[0]["value"]

    if not triggered_value:
        raise PreventUpdate

    if current_selected is None:
        current_selected = []

    if "synteny-clear-selected" in triggered_prop:
        return []

    try:
        triggered_id = json.loads(triggered_prop.split(".")[0])
    except (json.JSONDecodeError, ValueError):
        raise PreventUpdate

    action_type = triggered_id.get("type")
    ship_id = str(triggered_id.get("index"))

    if action_type == "synteny-add-ship":
        if len(current_selected) >= 4:
            raise PreventUpdate
        if any(str(s["id"]) == ship_id for s in current_selected):
            raise PreventUpdate
        ship = next((s for s in (all_ships or []) if str(s["value"]) == ship_id), None)
        if ship:
            return current_selected + [{"id": int(ship_id), "label": ship["label"]}]

    elif action_type == "synteny-remove-ship":
        return [s for s in current_selected if str(s["id"]) != ship_id]

    return current_selected


# ── Render selected sequences list ─────────────────────────────────────────────
@callback(
    [
        Output("synteny-selected-display", "children"),
        Output("synteny-selected-count", "children"),
    ],
    Input("synteny-selected-ships", "data"),
)
@handle_callback_error
def display_selected_ships(selected_ships):
    if not selected_ships:
        return dmc.Text("No sequences selected yet", c="dimmed", size="sm"), "0 / 4"

    items = []
    for i, ship in enumerate(selected_ships):
        color = SEQUENCE_COLORS[i % len(SEQUENCE_COLORS)]
        items.append(
            dmc.Group(
                [
                    dmc.Badge(str(i + 1), color=color, variant="filled", size="sm"),
                    dmc.Text(
                        ship["label"],
                        size="xs",
                        style={"flex": 1, "minWidth": 0, "wordBreak": "break-word"},
                    ),
                    dmc.ActionIcon(
                        DashIconify(icon="tabler:x"),
                        id={"type": "synteny-remove-ship", "index": str(ship["id"])},
                        variant="subtle",
                        color="red",
                        size="sm",
                    ),
                ],
                justify="space-between",
                gap="xs",
                wrap="nowrap",
            )
        )

    return dmc.Stack(items, gap="xs"), f"{len(selected_ships)} / 4"


# ── Generate visualization data ────────────────────────────────────────────────
@callback(
    Output("synteny-data-store", "data"),
    Input("synteny-update-button", "n_clicks"),
    [
        State("synteny-selected-ships", "data"),
        State("synteny-color-by", "value"),
        State("synteny-label-field", "value"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def generate_synteny_data(n_clicks, selected_ships, color_by, label_field):
    from src.utils.synteny_queries import get_gff_by_ship_ids

    if not n_clicks or not selected_ships:
        return no_update

    ship_ids = [s["id"] for s in selected_ships]

    if not ship_ids:
        return {"error": "Please select at least one sequence", "success": False}

    try:
        gff_data = get_gff_by_ship_ids(ship_ids)

        if not gff_data:
            return {
                "error": "No GFF data found for selected sequences",
                "success": False,
            }

        clusters = []
        for ship_id in ship_ids:
            ship_gff = [g for g in gff_data if g["ship_id"] == ship_id]

            if not ship_gff:
                continue

            genes = []
            for idx, gff_entry in enumerate(ship_gff):
                attrs = parse_gff_attributes(gff_entry["attributes"])
                category = categorize_gene(gff_entry)

                if color_by == "category":
                    group = category
                elif color_by == "family":
                    group = attrs.get("Target_ID", category)
                else:
                    group = gff_entry.get(color_by, "unknown")

                gene = {
                    "uid": f"gene_{ship_id}_{idx}",
                    "start": gff_entry["start"],
                    "end": gff_entry["end"],
                    "strand": 1 if gff_entry["strand"] == "+" else -1,
                    "_start": gff_entry["start"],
                    "_end": gff_entry["end"],
                    "_strand": 1 if gff_entry["strand"] == "+" else -1,
                    "name": attrs.get("Alias", attrs.get("Target_ID", f"gene_{idx}")),
                    "label": attrs.get(label_field, attrs.get("Alias", f"gene_{idx}")),
                    "locus_tag": attrs.get("SeqID", ""),
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
                        **attrs,
                    },
                    "group": group,
                    "colour": GENE_CATEGORY_COLORS.get(category)
                    if color_by == "category"
                    else None,
                }
                genes.append(gene)

            locus = {
                "uid": f"locus_{ship_id}",
                "name": ship_gff[0].get("ship_accession_display", f"Locus {ship_id}"),
                "start": min(g["start"] for g in genes),
                "end": max(g["end"] for g in genes),
                "_start": min(g["start"] for g in genes),
                "_end": max(g["end"] for g in genes),
                "genes": genes,
            }

            cluster = {
                "uid": f"cluster_{ship_id}",
                "name": f"{ship_gff[0]['ship_accession_display']}",
                "loci": [locus],
            }
            clusters.append(cluster)

        links = generate_synteny_links(clusters)
        logger.info(f"Generated {len(clusters)} clusters with {len(links)} links")

        return {
            "clusters": clusters,
            "links": links,
            "config": {"updateGroups": False},
            "success": True,
        }

    except Exception as e:
        logger.error(f"Error generating synteny data: {e}", exc_info=True)
        return {"error": f"Error: {str(e)}", "success": False}


# ── Helper functions ───────────────────────────────────────────────────────────


def parse_gff_attributes(attr_string):
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
    attrs = parse_gff_attributes(gff_entry.get("attributes", ""))
    searchable_text = " ".join(
        [
            attrs.get("Alias", ""),
            attrs.get("Target_ID", ""),
            attrs.get("product", ""),
            attrs.get("Name", ""),
            attrs.get("gene", ""),
            gff_entry.get("type", ""),
            gff_entry.get("source", ""),
        ]
    ).lower()

    # Use (?<![a-z])...(?![a-z]) instead of \b so that abbreviations embedded in
    # underscore-delimited gene names (e.g. "altals1_tyr65") still match.
    # Full words like "tyrosine" / "recombinase" are matched without boundary guards.
    if re.search(
        r"tyrosine|recombinase|capsid|(?<![a-z])(tyr|cap)(?![a-z])", searchable_text
    ):
        return "captain"
    elif re.search(r"(?<![a-z])nlr(?![a-z])", searchable_text):
        return "nlr"
    elif re.search(r"(?<![a-z])plp(?![a-z])", searchable_text):
        return "plp"
    elif re.search(r"(?<![a-z])fre(?![a-z])", searchable_text):
        return "fre"
    else:
        return "other"


def generate_synteny_links(clusters):
    links = []
    link_id = 0

    target_index = {}
    category_index = {}

    for cluster in clusters:
        for locus in cluster["loci"]:
            for gene in locus["genes"]:
                target_id = gene["metadata"].get("Target_ID")
                if target_id and target_id != "":
                    if target_id not in target_index:
                        target_index[target_id] = []
                    target_index[target_id].append(gene)

                category = gene["metadata"].get("category", "other")
                if category not in category_index:
                    category_index[category] = []
                category_index[category].append(gene)

    for target_id, genes in target_index.items():
        if len(genes) > 1:
            for i in range(len(genes)):
                for j in range(i + 1, len(genes)):
                    gene_i_cat = genes[i]["metadata"].get("category", "other")
                    gene_j_cat = genes[j]["metadata"].get("category", "other")
                    identity = 0.95 if gene_i_cat == gene_j_cat else 0.85
                    links.append(
                        {
                            "uid": f"link_{link_id}",
                            "query": {"uid": genes[i]["uid"]},
                            "target": {"uid": genes[j]["uid"]},
                            "identity": identity,
                        }
                    )
                    link_id += 1

    linked_genes = set()
    for link in links:
        linked_genes.add(link["query"]["uid"])
        linked_genes.add(link["target"]["uid"])

    for category in ["captain", "nlr", "plp", "fre"]:
        if category in category_index:
            genes = [
                g for g in category_index[category] if g["uid"] not in linked_genes
            ]
            if len(genes) > 1:
                for i in range(len(genes)):
                    for j in range(i + 1, len(genes)):
                        links.append(
                            {
                                "uid": f"link_{link_id}",
                                "query": {"uid": genes[i]["uid"]},
                                "target": {"uid": genes[j]["uid"]},
                                "identity": 0.6,
                            }
                        )
                        link_id += 1

    return links


# ── Clientside: render ClusterMap visualization ────────────────────────────────
clientside_callback(
    """
    function(store_data, scale_factor, cluster_spacing, show_links, 
             identity_threshold, clear_clicks) {
        
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
                    alignLabels: false,
                },
                gene: {
                    label: {
                        show: false,
                    },
                    shape: {
                        onClick: function(event, gene) {
                            if (gene && gene.metadata) {
                                const panel = document.getElementById('gene-details-panel');
                                if (panel) {
                                    let html = '<div style="font-family: system-ui; font-size: 13px;">';
                                    html += '<div style="margin-bottom: 16px;">';
                                    html += '<div style="font-size: 18px; font-weight: bold; color: #228be6; margin-bottom: 4px;">';
                                    html += (gene.label || gene.name || 'Unknown Gene');
                                    html += '</div>';
                                    if (gene.metadata.category_label) {
                                        const categoryColors = {
                                            'captain': '#e74c3c', 'nlr': '#3498db',
                                            'plp': '#2ecc71', 'fre': '#f39c12', 'other': '#95a5a6'
                                        };
                                        const bgColor = categoryColors[gene.metadata.category] || '#95a5a6';
                                        html += '<span style="display: inline-block; background: ' + bgColor + '; color: white; ';
                                        html += 'padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: 600;">';
                                        html += gene.metadata.category_label + '</span>';
                                    }
                                    html += '</div>';
                                    html += '<div style="background: #f8f9fa; padding: 12px; border-radius: 6px; margin-bottom: 12px;">';
                                    html += '<div style="font-weight: bold; margin-bottom: 8px; color: #495057;">📍 Genomic Location</div>';
                                    html += '<div style="font-size: 12px;">';
                                    html += '<div style="margin-bottom: 4px;"><strong>Position:</strong> ' +
                                        gene.start.toLocaleString() + ' - ' + gene.end.toLocaleString() + '</div>';
                                    html += '<div style="margin-bottom: 4px;"><strong>Length:</strong> ' +
                                        (gene.end - gene.start).toLocaleString() + ' bp</div>';
                                    html += '<div><strong>Strand:</strong> ' + (gene.strand > 0 ? 'Forward (+)' : 'Reverse (-)') + '</div>';
                                    html += '</div></div>';
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
                legend: { show: false },
                link: {
                    show: show_links !== false,
                    threshold: identity_threshold || 0.3,
                }
            };
            
            const container = document.getElementById(containerId);
            if (!container) return ["Error: Container not found", {"display": "none"}];
            
            container.innerHTML = '';
            container.style.display = 'block';
            container.style.alignItems = 'unset';
            container.style.justifyContent = 'unset';
            
            if (typeof d3 !== 'undefined') {
                d3.selectAll('#' + containerId + ' *').remove();
            }
            
            setTimeout(function() {
                try {
                    const chart = ClusterMap.ClusterMap().config(config);
                    d3.select('#' + containerId)
                        .datum(store_data)
                        .call(chart);
                    
                    if (typeof d3 !== 'undefined') {
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
                                    html += '<div style="font-weight: bold; font-size: 13px; margin-bottom: 6px; color: #228be6;">';
                                    html += (d.label || d.name || 'Unknown Gene') + '</div>';
                                    if (d.metadata.category_label) {
                                        const categoryColors = {
                                            'captain': '#e74c3c', 'nlr': '#3498db',
                                            'plp': '#2ecc71', 'fre': '#f39c12', 'other': '#95a5a6'
                                        };
                                        const bgColor = categoryColors[d.metadata.category] || '#95a5a6';
                                        html += '<div style="margin-bottom: 6px;">';
                                        html += '<span style="background: ' + bgColor + '; color: white; padding: 2px 6px; ';
                                        html += 'border-radius: 10px; font-size: 10px; font-weight: 600;">';
                                        html += d.metadata.category_label + '</span></div>';
                                    }
                                    html += '<div style="border-top: 1px solid #dee2e6; padding-top: 6px; margin-top: 6px;">';
                                    html += '<strong>Position:</strong> ' + d.start.toLocaleString() + ' - ' + d.end.toLocaleString();
                                    html += ' (' + (d.end - d.start).toLocaleString() + ' bp)<br/>';
                                    html += '<strong>Strand:</strong> ' + (d.strand > 0 ? 'Forward (+)' : 'Reverse (-)') + '<br/>';
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
        Input("synteny-identity-threshold", "value"),
        Input("synteny-clear-viz", "n_clicks"),
    ],
)


# ── SVG export ─────────────────────────────────────────────────────────────────
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
    prevent_initial_call=True,
)


# ── Reset view ─────────────────────────────────────────────────────────────────
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
    prevent_initial_call=True,
)


# ── Clear visualization ────────────────────────────────────────────────────────
clientside_callback(
    """
    function(n_clicks) {
        if (n_clicks) { return null; }
        return window.dash_clientside.no_update;
    }
    """,
    Output("synteny-data-store", "data", allow_duplicate=True),
    Input("synteny-clear-viz", "n_clicks"),
    prevent_initial_call=True,
)

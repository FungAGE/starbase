from datetime import datetime
from urllib.parse import parse_qs

import dash
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash import (
    ALL,
    Input,
    Output,
    State,
    callback,
    clientside_callback,
    dcc,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from src.components.data import create_ship_accession_modal_data, render_ship_accession_modal
from src.config.logging import get_logger
from src.config.settings import ADMIN_TOKEN
from src.database import curation_manager
from src.database.curation_constants import FLAG_COLORS, FLAG_LABELS
from src.database.sql_manager import fetch_ships
from src.utils.seq_utils import create_ncbi_style_header

logger = get_logger(__name__)

dash.register_page(__name__, path="/curation", title="Curation", name="Curation")

FLAG_FILTER_OPTIONS = [{"value": "all", "label": "All"}] + [
    {"value": str(k), "label": v} for k, v in FLAG_LABELS.items()
]
FLAG_EDIT_OPTIONS = [{"value": str(k), "label": v} for k, v in FLAG_LABELS.items()]

RESULT_TOOLS = [
    ("blastp", "Blastp"),
    ("rpsblast", "RPS-BLAST"),
    ("hhsearch", "HHsearch"),
    ("interpro", "InterProScan"),
]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def _queue_rows(flag_value):
    flag = None if not flag_value or flag_value == "all" else int(flag_value)
    rows = curation_manager.fetch_annotation_queue(flag=flag, limit=200)
    for r in rows:
        r["annotation"] = r["annotation"] or ""
        r["assigned_to"] = r["assigned_to"] or ""
    return rows


def _flag_badge(flag):
    label = FLAG_LABELS.get(flag, "UNKNOWN")
    color = FLAG_COLORS.get(flag, "gray")
    return dmc.Badge(label, color=color, variant="light")


def _history_table(history):
    if not history:
        return dmc.Text("No edits yet.", size="sm", c="dimmed")
    rows = []
    for h in reversed(history):
        changes = []
        if h["new_flag"] is not None:
            changes.append(
                f"flag: {FLAG_LABELS.get(h['old_flag'], h['old_flag'])} → "
                f"{FLAG_LABELS.get(h['new_flag'], h['new_flag'])}"
            )
        if h["new_annotation"] is not None:
            changes.append(
                f"name: {h['old_annotation'] or ''!r} → {h['new_annotation']!r}"
            )
        if h["new_public_notes"] is not None:
            changes.append("public notes edited")
        if h["new_private_notes"] is not None:
            changes.append("private notes edited")
        rows.append(
            html.Tr(
                [
                    html.Td(h["changed_at"], style={"whiteSpace": "nowrap"}),
                    html.Td(h["changed_by"] or ""),
                    html.Td("; ".join(changes)),
                ]
            )
        )
    return dmc.Table(
        [
            html.Thead(html.Tr([html.Th("When"), html.Th("Who"), html.Th("What")])),
            html.Tbody(rows),
        ],
        striped=True,
        highlightOnHover=True,
        fz="sm",
    )


def _gene_features_table(features):
    if not features:
        return dmc.Text("No linked gene features.", size="sm", c="dimmed")
    rows = [
        html.Tr(
            [
                html.Td(f["joined_ship_id"]),
                html.Td(f["type"] or ""),
                html.Td(f"{f['start']}–{f['stop']}"),
                html.Td(f["strand"] or ""),
            ]
        )
        for f in features
    ]
    return dmc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Ship (joined_ship_id)"),
                        html.Th("Type"),
                        html.Th("Range"),
                        html.Th("Strand"),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        striped=True,
        highlightOnHover=True,
        fz="sm",
    )


def _ship_gene_features_panel(features):
    """Linked gene features for a ship-detail panel -- each row jumps into the
    Annotation queue tab if it has a linked Annotation, plain text otherwise."""
    if not features:
        return dmc.Text("No gene features called on this ship yet.", size="sm", c="dimmed")
    rows = []
    for f in features:
        if f["annotation_id"] is not None:
            name_cell = dmc.Anchor(
                f["annotation"] or f"Annotation #{f['annotation_id']}",
                id={"type": "curation-jump-to-annotation", "index": f["annotation_id"]},
                size="sm",
            )
        else:
            name_cell = dmc.Text("—", size="sm", c="dimmed")
        rows.append(
            html.Tr(
                [
                    html.Td(f["type"] or ""),
                    html.Td(f"{f['start']}–{f['stop']}"),
                    html.Td(f["strand"] or ""),
                    html.Td(_flag_badge(f["flag"]) if f["flag"] is not None else ""),
                    html.Td(name_cell),
                ]
            )
        )
    return dmc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Type"),
                        html.Th("Range"),
                        html.Th("Strand"),
                        html.Th("Flag"),
                        html.Th("Annotation"),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        striped=True,
        highlightOnHover=True,
        fz="sm",
    )


def _results_panel(results):
    panels = []
    for tool_key, tool_label in RESULT_TOOLS:
        rows = (results or {}).get(tool_key) or []
        if not rows:
            body = dmc.Text("No results yet.", size="sm", c="dimmed")
        else:
            body = dmc.Stack(
                [
                    dmc.Text(
                        f"{r['database']} — {r['status']} ({r['run_date'] or 'no run date'})",
                        size="sm",
                    )
                    for r in rows
                ]
            )
        panels.append(
            dmc.Paper(
                dmc.Stack([dmc.Text(tool_label, fw=600, size="sm"), body], gap="xs"),
                p="sm",
                withBorder=True,
                radius="sm",
            )
        )
    return dmc.SimpleGrid(panels, cols=2, spacing="sm")


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def _build_curation_layout():
    return html.Div(
        dmc.Tabs(
            [
                dmc.TabsList(
                    [
                        dmc.TabsTab("Ships overview", value="ships"),
                        dmc.TabsTab("Annotation queue", value="annotations"),
                    ]
                ),
                dmc.TabsPanel(_build_ships_tab(), value="ships"),
                dmc.TabsPanel(_build_annotations_tab(), value="annotations"),
            ],
            id="curation-tabs",
            value="ships",
        ),
        style={"padding": "1rem 2rem"},
    )


def _build_annotations_tab():
    return html.Div(
        [
            dmc.Group(
                [
                    dmc.Title("Annotation Curation", order=3),
                    dmc.Select(
                        id="curation-flag-filter",
                        data=FLAG_FILTER_OPTIONS,
                        value="all",
                        w=200,
                        allowDeselect=False,
                    ),
                ],
                justify="space-between",
                my="md",
            ),
            html.Div(id="curation-queue-section"),
            html.Div(id="curation-detail-section", style={"display": "none"}),
        ]
    )


def _build_ships_tab():
    return html.Div(
        [
            dmc.Group(
                [dmc.Title("Ships Overview", order=3)],
                justify="space-between",
                my="md",
            ),
            html.Div(id="curation-ships-queue-section"),
            html.Div(id="curation-ship-detail-section", style={"display": "none"}),
            dcc.Download(id="curation-ships-dl"),
        ]
    )


def _build_queue_section():
    return dag.AgGrid(
        id="curation-queue-grid",
        columnDefs=[
            {"field": "id", "headerName": "ID", "width": 90, "pinned": "left"},
            {"field": "annotation", "headerName": "Name", "flex": 2},
            {"field": "flag_label", "headerName": "Flag", "width": 140},
            {"field": "assigned_to", "headerName": "Assigned to", "width": 160},
            {"field": "updated_at", "headerName": "Updated", "width": 170},
        ],
        rowData=[],
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": 25,
            "rowSelection": "single",
            "rowHeight": 40,
        },
        style={"height": "60vh"},
        className="ag-theme-alpine",
    )


def _build_detail_section():
    return dmc.Stack(
        [
            dmc.Group(
                [
                    dmc.Button(
                        "← Back to queue",
                        id="curation-back-btn",
                        variant="subtle",
                        size="sm",
                    ),
                    dmc.Group(
                        [
                            dmc.Button(
                                "Save", id="curation-save-btn", color="blue", size="sm"
                            ),
                            dmc.Button(
                                "Next unreviewed →",
                                id="curation-next-btn",
                                variant="outline",
                                size="sm",
                            ),
                        ],
                        gap="xs",
                    ),
                ],
                justify="space-between",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Group(
                            [
                                dmc.Text(id="curation-detail-id", fw=700),
                                html.Div(id="curation-detail-flag-badge"),
                            ],
                            gap="sm",
                        ),
                        dmc.Textarea(
                            id="curation-sequence-display",
                            label="Sequence",
                            autosize=True,
                            minRows=2,
                            maxRows=6,
                            readOnly=True,
                        ),
                        dmc.Group(
                            [
                                dmc.TextInput(
                                    id="curation-annotation-input",
                                    label="Name / product",
                                    w="60%",
                                ),
                                dmc.Select(
                                    id="curation-flag-select",
                                    label="Flag",
                                    data=FLAG_EDIT_OPTIONS,
                                    w="35%",
                                    allowDeselect=False,
                                ),
                            ],
                            grow=True,
                        ),
                        dmc.TextInput(
                            id="curation-assigned-to-input", label="Assigned to"
                        ),
                        dmc.Textarea(
                            id="curation-public-notes",
                            label="Public notes",
                            autosize=True,
                            minRows=2,
                            maxRows=8,
                        ),
                        dmc.Textarea(
                            id="curation-private-notes",
                            label="Private notes (curator-only)",
                            autosize=True,
                            minRows=2,
                            maxRows=8,
                        ),
                    ],
                    gap="sm",
                ),
                p="md",
                withBorder=True,
                radius="sm",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Text("Gene map", fw=600, size="sm"),
                        html.Div(id="curation-ship-viz", style={"minHeight": "260px"}),
                    ],
                    gap="xs",
                ),
                p="sm",
                withBorder=True,
                radius="sm",
            ),
            dmc.Accordion(
                [
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Linked gene features"),
                            dmc.AccordionPanel(html.Div(id="curation-gene-features")),
                        ],
                        value="features",
                    ),
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Search results"),
                            dmc.AccordionPanel(html.Div(id="curation-results")),
                        ],
                        value="results",
                    ),
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Edit history"),
                            dmc.AccordionPanel(html.Div(id="curation-history")),
                        ],
                        value="history",
                    ),
                ],
                multiple=True,
                value=["features"],
            ),
        ],
        gap="md",
    )


def _build_ships_queue_section():
    return html.Div(
        [
            dmc.Group(
                dmc.Button(
                    "Download FASTA",
                    id="curation-ships-download-btn",
                    variant="outline",
                    size="sm",
                ),
                justify="flex-end",
                mb="xs",
            ),
            dag.AgGrid(
                id="curation-ships-grid",
                columnDefs=[
                    {"field": "id", "headerName": "ID", "width": 90, "pinned": "left"},
                    {"field": "starshipID", "headerName": "Starship ID", "flex": 1},
                    {"field": "accession_tag", "headerName": "Accession", "width": 140},
                    {"field": "familyName", "headerName": "Family", "width": 120},
                    {"field": "navis_name", "headerName": "Navis", "width": 120},
                    {"field": "haplotype_name", "headerName": "Haplotype", "width": 120},
                    {"field": "taxonomy_name", "headerName": "Organism", "flex": 1},
                    {"field": "curated_status", "headerName": "Curated status", "width": 140},
                    {"field": "evidence", "headerName": "Evidence", "width": 110},
                ],
                rowData=[],
                dashGridOptions={
                    "pagination": True,
                    "paginationPageSize": 25,
                    "rowSelection": "single",
                    "rowHeight": 40,
                },
                style={"height": "60vh"},
                className="ag-theme-alpine",
            ),
        ]
    )


def _build_ship_detail_section():
    return dmc.Stack(
        [
            dmc.Group(
                [
                    dmc.Button(
                        "← Back to ships",
                        id="curation-ships-back-btn",
                        variant="subtle",
                        size="sm",
                    ),
                ],
                justify="space-between",
            ),
            dmc.Paper(
                html.Div(id="curation-ship-meta"),
                p="md",
                withBorder=True,
                radius="sm",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Text("Gene map", fw=600, size="sm"),
                        html.Div(
                            id="curation-ship-detail-viz", style={"minHeight": "260px"}
                        ),
                    ],
                    gap="xs",
                ),
                p="sm",
                withBorder=True,
                radius="sm",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Text("Linked gene features", fw=600, size="sm"),
                        html.Div(id="curation-ship-gene-features-list"),
                    ],
                    gap="xs",
                ),
                p="sm",
                withBorder=True,
                radius="sm",
            ),
        ],
        gap="md",
    )


layout = html.Div(
    [
        dcc.Location(id="curation-url", refresh=False),
        dcc.Store(id="curation-selected-id", data=None),
        dcc.Store(id="curation-selected-ship-id", data=None),
        dcc.Store(id="curation-authed", data=False),
        dcc.Store(id="curation-ship-viz-data", data=None),
        dcc.Store(id="curation-ship-detail-viz-data", data=None),
        html.Div(id="curation-content"),
    ]
)


# ---------------------------------------------------------------------------
# Callbacks — token gate
# ---------------------------------------------------------------------------


@callback(
    [Output("curation-content", "children"), Output("curation-authed", "data")],
    Input("curation-url", "search"),
)
def render_curation_content(search):
    params = parse_qs((search or "").lstrip("?"))
    token = params.get("token", [None])[0]
    if not ADMIN_TOKEN or not token or token != ADMIN_TOKEN:
        return (
            dmc.Container(
                dmc.Alert(
                    "A valid ?token= query parameter is required.",
                    title="Unauthorized",
                    color="red",
                    mt="xl",
                ),
                size="sm",
            ),
            False,
        )
    return _build_curation_layout(), True


# ---------------------------------------------------------------------------
# Callbacks — queue
# ---------------------------------------------------------------------------


@callback(
    Output("curation-queue-section", "children"),
    Input("curation-authed", "data"),
)
def init_queue_section(authed):
    if not authed:
        raise PreventUpdate
    return _build_queue_section()


@callback(
    Output("curation-queue-grid", "rowData"),
    Input("curation-flag-filter", "value"),
    Input("curation-authed", "data"),
    prevent_initial_call=True,
)
def refresh_queue(flag_value, authed):
    if not authed:
        raise PreventUpdate
    return _queue_rows(flag_value)


clientside_callback(
    """
    function(selectedRows) {
        if (!selectedRows || selectedRows.length === 0) {
            return window.dash_clientside.no_update;
        }
        return selectedRows[0].id;
    }
    """,
    Output("curation-selected-id", "data"),
    Input("curation-queue-grid", "selectedRows"),
    prevent_initial_call=True,
)


# assets/js/starship_visualization.js defines window.renderStarshipViz.
# Output is a no-op style re-set on the same container -- clientside
# callbacks need a real Output; the render itself is the side effect.
clientside_callback(
    """
    function(vizData) {
        if (!vizData) {
            return {"minHeight": "260px"};
        }
        if (window.renderStarshipViz) {
            window.renderStarshipViz("curation-ship-viz", vizData);
        }
        return {"minHeight": "260px"};
    }
    """,
    Output("curation-ship-viz", "style"),
    Input("curation-ship-viz-data", "data"),
    prevent_initial_call=True,
)


clientside_callback(
    """
    function(vizData) {
        if (!vizData) {
            return {"minHeight": "260px"};
        }
        if (window.renderStarshipViz) {
            window.renderStarshipViz("curation-ship-detail-viz", vizData);
        }
        return {"minHeight": "260px"};
    }
    """,
    Output("curation-ship-detail-viz", "style"),
    Input("curation-ship-detail-viz-data", "data"),
    prevent_initial_call=True,
)


@callback(
    Output("curation-selected-id", "data", allow_duplicate=True),
    Input("curation-back-btn", "n_clicks"),
    prevent_initial_call=True,
)
def go_back_to_queue(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return None


# ---------------------------------------------------------------------------
# Callbacks — section toggling
# ---------------------------------------------------------------------------


@callback(
    [
        Output("curation-queue-section", "style"),
        Output("curation-detail-section", "style"),
        Output("curation-detail-section", "children"),
    ],
    Input("curation-selected-id", "data"),
    prevent_initial_call=True,
)
def toggle_sections(selected_id):
    if selected_id is None:
        return {"display": "block"}, {"display": "none"}, no_update
    return {"display": "none"}, {"display": "block"}, _build_detail_section()


# ---------------------------------------------------------------------------
# Callbacks — detail load
# ---------------------------------------------------------------------------


@callback(
    [
        Output("curation-detail-id", "children"),
        Output("curation-detail-flag-badge", "children"),
        Output("curation-sequence-display", "value"),
        Output("curation-annotation-input", "value"),
        Output("curation-flag-select", "value"),
        Output("curation-assigned-to-input", "value"),
        Output("curation-public-notes", "value"),
        Output("curation-private-notes", "value"),
        Output("curation-gene-features", "children"),
        Output("curation-results", "children"),
        Output("curation-history", "children"),
        Output("curation-ship-viz-data", "data"),
    ],
    Input("curation-selected-id", "data"),
    Input("curation-detail-section", "children"),
    prevent_initial_call=True,
)
def load_detail(selected_id, _section_children):
    if selected_id is None:
        raise PreventUpdate
    a = curation_manager.fetch_annotation(selected_id)

    viz_data = None
    if a["gene_features"]:
        joined_ship_id = a["gene_features"][0]["joined_ship_id"]
        try:
            ship_features = curation_manager.fetch_ship_gene_features(joined_ship_id)
            viz_data = {**ship_features, "current_annotation_id": selected_id}
        except Exception as exc:
            logger.error("Failed to load gene map for ship %s: %s", joined_ship_id, exc)

    return (
        f"Annotation #{a['id']}",
        _flag_badge(a["flag"]),
        a["sequence"],
        a["annotation"] or "",
        str(a["flag"]),
        a["assigned_to"] or "",
        a["public_notes"] or "",
        a["private_notes"] or "",
        _gene_features_table(a["gene_features"]),
        _results_panel(a["results"]),
        _history_table(a["history"]),
        viz_data,
    )


# ---------------------------------------------------------------------------
# Callbacks — save / next
# ---------------------------------------------------------------------------


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("curation-selected-id", "data", allow_duplicate=True),
    Input("curation-save-btn", "n_clicks"),
    State("curation-selected-id", "data"),
    State("curation-annotation-input", "value"),
    State("curation-flag-select", "value"),
    State("curation-assigned-to-input", "value"),
    State("curation-public-notes", "value"),
    State("curation-private-notes", "value"),
    prevent_initial_call=True,
)
def save_annotation(
    n_clicks, selected_id, annotation, flag, assigned_to, public_notes, private_notes
):
    if not n_clicks or selected_id is None:
        raise PreventUpdate
    changes = {
        "annotation": annotation,
        "flag": int(flag),
        "assigned_to": assigned_to,
        "public_notes": public_notes,
        "private_notes": private_notes,
    }
    try:
        curation_manager.update_annotation(selected_id, changes, changed_by="curator")
        notif = dmc.Notification(
            id="curation-save-ok",
            title="Saved",
            message=f"Annotation #{selected_id} updated.",
            color="green",
            action="show",
            autoClose=4000,
        )
    except Exception as exc:
        logger.error("Failed to save annotation %s: %s", selected_id, exc)
        notif = dmc.Notification(
            id="curation-save-err",
            title="Save failed",
            message=str(exc),
            color="red",
            action="show",
            autoClose=8000,
        )
        return notif, no_update
    # re-trigger detail reload with the same id
    return notif, selected_id


@callback(
    Output("curation-selected-id", "data", allow_duplicate=True),
    Input("curation-next-btn", "n_clicks"),
    State("curation-selected-id", "data"),
    State("curation-flag-filter", "value"),
    prevent_initial_call=True,
)
def go_next_unreviewed(n_clicks, selected_id, flag_value):
    if not n_clicks:
        raise PreventUpdate
    rows = _queue_rows(flag_value)
    ids = [r["id"] for r in rows]
    if not ids:
        raise PreventUpdate
    if selected_id in ids:
        idx = ids.index(selected_id)
        next_idx = (idx + 1) % len(ids)
    else:
        next_idx = 0
    return ids[next_idx]


# ---------------------------------------------------------------------------
# Callbacks — ships overview
# ---------------------------------------------------------------------------


@callback(
    Output("curation-ships-queue-section", "children"),
    Input("curation-authed", "data"),
)
def init_ships_queue_section(authed):
    if not authed:
        raise PreventUpdate
    return _build_ships_queue_section()


@callback(
    Output("curation-ships-grid", "rowData"),
    Input("curation-authed", "data"),
    prevent_initial_call=True,
)
def refresh_ships_grid(authed):
    if not authed:
        raise PreventUpdate
    return curation_manager.fetch_ships_overview()


clientside_callback(
    """
    function(selectedRows) {
        if (!selectedRows || selectedRows.length === 0) {
            return window.dash_clientside.no_update;
        }
        return selectedRows[0];
    }
    """,
    Output("curation-selected-ship-id", "data"),
    Input("curation-ships-grid", "selectedRows"),
    prevent_initial_call=True,
)


@callback(
    Output("curation-selected-ship-id", "data", allow_duplicate=True),
    Input("curation-ships-back-btn", "n_clicks"),
    prevent_initial_call=True,
)
def go_back_to_ships(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return None


@callback(
    [
        Output("curation-ships-queue-section", "style"),
        Output("curation-ship-detail-section", "style"),
        Output("curation-ship-detail-section", "children"),
    ],
    Input("curation-selected-ship-id", "data"),
    prevent_initial_call=True,
)
def toggle_ship_sections(selected_ship_id):
    if selected_ship_id is None:
        return {"display": "block"}, {"display": "none"}, no_update
    return {"display": "none"}, {"display": "block"}, _build_ship_detail_section()


@callback(
    [
        Output("curation-ship-meta", "children"),
        Output("curation-ship-detail-viz-data", "data"),
        Output("curation-ship-gene-features-list", "children"),
    ],
    Input("curation-selected-ship-id", "data"),
    Input("curation-ship-detail-section", "children"),
    prevent_initial_call=True,
)
def load_ship_detail(selected_ship, _section_children):
    if selected_ship is None:
        raise PreventUpdate

    joined_ship_id = selected_ship["id"]
    ship_features = curation_manager.fetch_ship_gene_features(joined_ship_id)

    ship_accession_tag = selected_ship.get("ship_accession_tag") or selected_ship.get(
        "accession_tag"
    )

    if ship_accession_tag:
        meta = render_ship_accession_modal(
            create_ship_accession_modal_data(ship_accession_tag)
        )
    else:
        meta = dmc.Alert(
            "No accession tag on this ship -- can't load metadata.",
            title="No metadata",
            color="yellow",
        )

    return (
        meta,
        ship_features,
        _ship_gene_features_panel(ship_features["features"]),
    )


def _ship_download_payload(row):
    """Build a single-sequence FASTA download dict for one ships-overview row,
    reusing the same primitives as wiki.py's generate_download_helper."""
    tag = row.get("ship_accession_tag") or row.get("accession_tag")
    if not tag:
        raise ValueError("Selected ship has no accession tag to download.")

    dl_df = fetch_ships(accessions=[tag], curated=False, dereplicate=True, with_sequence=True)
    if dl_df is None or dl_df.empty:
        raise ValueError(f"No sequence found for {tag}.")

    seq_row = dl_df.iloc[0]
    header = create_ncbi_style_header(seq_row)
    if header is None:
        raise ValueError(f"Could not build a FASTA header for {tag}.")

    fasta_str = f"{header}\n{seq_row['sequence']}"
    return dict(
        content=fasta_str,
        filename=f"{tag}_{datetime.now().strftime('%y%m%d')}.fasta",
        type="text/plain",
    )


@callback(
    Output("curation-ships-dl", "data"),
    Output("notifications-container", "children", allow_duplicate=True),
    Input("curation-ships-download-btn", "n_clicks"),
    State("curation-ships-grid", "selectedRows"),
    prevent_initial_call=True,
)
def download_selected_ship(n_clicks, selected_rows):
    if not n_clicks:
        raise PreventUpdate
    if not selected_rows:
        raise PreventUpdate
    try:
        return _ship_download_payload(selected_rows[0]), no_update
    except Exception as exc:
        logger.error("Ship FASTA download failed: %s", exc)
        notif = dmc.Notification(
            id="curation-download-err",
            title="Download failed",
            message=str(exc),
            color="red",
            action="show",
            autoClose=8000,
        )
        return no_update, notif


@callback(
    Output("curation-selected-id", "data", allow_duplicate=True),
    Output("curation-tabs", "value", allow_duplicate=True),
    Input({"type": "curation-jump-to-annotation", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def jump_to_annotation(n_clicks_list):
    triggered_id = dash.callback_context.triggered_id
    if not triggered_id or not any(n_clicks_list):
        raise PreventUpdate
    return triggered_id["index"], "annotations"

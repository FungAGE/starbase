from urllib.parse import parse_qs

import dash
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash import (
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

from src.config.logging import get_logger
from src.config.settings import ADMIN_TOKEN
from src.database import curation_manager
from src.database.curation_constants import FLAG_COLORS, FLAG_LABELS

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
                mb="md",
            ),
            html.Div(id="curation-queue-section"),
            html.Div(id="curation-detail-section", style={"display": "none"}),
        ],
        style={"padding": "1rem 2rem"},
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


layout = html.Div(
    [
        dcc.Location(id="curation-url", refresh=False),
        dcc.Store(id="curation-selected-id", data=None),
        dcc.Store(id="curation-authed", data=False),
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
    ],
    Input("curation-selected-id", "data"),
    Input("curation-detail-section", "children"),
    prevent_initial_call=True,
)
def load_detail(selected_id, _section_children):
    if selected_id is None:
        raise PreventUpdate
    a = curation_manager.fetch_annotation(selected_id)
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

from urllib.parse import parse_qs
import csv
import io
import time

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
from src.database import starfish_manager

logger = get_logger(__name__)

dash.register_page(__name__, path="/starfish", title="Starfish", name="Starfish")

STATUS_COLORS = {
    "pending": "gray",
    "running": "yellow",
    "completed": "green",
    "failed": "red",
    "cancelled": "gray",
}

MODEL_OPTIONS = [{"value": "tyr", "label": "tyr"}, {"value": "cap", "label": "cap"}]

_SAMPLESHEET_COLS = ["genomeID", "taxID", "fna", "gff3", "emapper", "cds", "faa"]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def _status_badge(status):
    return dmc.Badge(status, color=STATUS_COLORS.get(status, "gray"), variant="light")


def _list_rows():
    runs = starfish_manager.list_runs(limit=200)
    for r in runs:
        r["status_label"] = r["status"]
        r["duration_display"] = (
            f"{r['duration_seconds']:.0f}s" if r.get("duration_seconds") else ""
        )
    return runs


def _parse_genomes_csv(text):
    """Parse the pasted samplesheet CSV into the dict shape create_run
    expects. genomeID/fna/gff3 required, matching backend validation."""
    if not text or not text.strip():
        raise ValueError("Paste at least one genome row")
    reader = csv.DictReader(io.StringIO(text.strip()))
    genomes = []
    for row in reader:
        genome_id = (row.get("genomeID") or "").strip()
        fna = (row.get("fna") or "").strip()
        gff3 = (row.get("gff3") or "").strip()
        if not (genome_id and fna and gff3):
            continue
        genomes.append(
            {
                "genome_id": genome_id,
                "fna_path": fna,
                "gff3_path": gff3,
                "tax_id": (row.get("taxID") or "").strip() or None,
                "emapper_path": (row.get("emapper") or "").strip() or None,
                "cds_path": (row.get("cds") or "").strip() or None,
                "faa_path": (row.get("faa") or "").strip() or None,
            }
        )
    if not genomes:
        raise ValueError(
            f"No valid genome rows found -- header must include {','.join(_SAMPLESHEET_COLS)}"
        )
    return genomes


def _genomes_table(genomes):
    if not genomes:
        return dmc.Text("No genomes on this run.", size="sm", c="dimmed")
    rows = [
        html.Tr(
            [
                html.Td(g["genome_id"]),
                html.Td(g.get("tax_id") or ""),
                html.Td(_status_badge(g["status"])),
                html.Td(
                    g.get("num_elements") if g.get("num_elements") is not None else ""
                ),
                html.Td(
                    g.get("error_message") or "",
                    style={"color": "var(--mantine-color-red-6)"},
                ),
            ]
        )
        for g in genomes
    ]
    return dmc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Genome ID"),
                        html.Th("Tax ID"),
                        html.Th("Status"),
                        html.Th("Elements"),
                        html.Th("Error"),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        striped=True,
        highlightOnHover=True,
        fz="sm",
    )


def _elements_table(elements):
    if not elements:
        return dmc.Text(
            "No elements found yet -- run the pipeline to detect starship elements.",
            size="sm",
            c="dimmed",
        )
    rows = [
        html.Tr(
            [
                html.Td(e["element_id"]),
                html.Td(e["contig_id"] or ""),
                html.Td(f"{e['start']}-{e['end']}"),
                html.Td(e["length"]),
                html.Td(e.get("family") or ""),
                html.Td(
                    dmc.Badge("imported", color="green", variant="light")
                    if e.get("imported_submission_id")
                    else dmc.Badge("not imported", color="gray", variant="light")
                ),
            ]
        )
        for e in elements
    ]
    return dmc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th("Element ID"),
                        html.Th("Contig"),
                        html.Th("Position"),
                        html.Th("Length"),
                        html.Th("Family"),
                        html.Th("Import status"),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        striped=True,
        highlightOnHover=True,
        fz="sm",
    )


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


def _build_starfish_layout():
    return html.Div(
        [
            dmc.Group(
                [
                    dmc.Title("Starfish Runs", order=3),
                    dmc.Button("New Run", id="starfish-new-run-btn", color="blue"),
                ],
                justify="space-between",
                mb="md",
            ),
            html.Div(id="starfish-list-section"),
            html.Div(id="starfish-detail-section", style={"display": "none"}),
            _create_run_modal(),
        ],
        style={"padding": "1rem 2rem"},
    )


def _create_run_modal():
    return dmc.Modal(
        id="starfish-create-modal",
        title="New Starfish Run",
        size="lg",
        opened=False,
        children=[
            dmc.Stack(
                [
                    dmc.TextInput(
                        id="starfish-create-run-name", label="Run name", required=True
                    ),
                    dmc.Textarea(
                        id="starfish-create-description",
                        label="Description",
                        autosize=True,
                        minRows=1,
                    ),
                    dmc.Group(
                        [
                            dmc.Select(
                                id="starfish-create-model",
                                label="Model",
                                data=MODEL_OPTIONS,
                                value="tyr",
                                w=120,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-threads",
                                label="Threads",
                                value=20,
                                min=1,
                                w=110,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-missing",
                                label="Missing",
                                value=1,
                                min=0,
                                w=110,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-maxcopy",
                                label="Max copy",
                                value=5,
                                min=1,
                                w=110,
                            ),
                        ],
                        grow=True,
                    ),
                    dmc.Group(
                        [
                            dmc.NumberInput(
                                id="starfish-create-pid",
                                label="PID %",
                                value=90,
                                min=0,
                                max=100,
                                w=110,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-hsp",
                                label="HSP",
                                value=1000,
                                min=1,
                                w=110,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-flank",
                                label="Flank",
                                value=6,
                                min=0,
                                w=110,
                            ),
                            dmc.NumberInput(
                                id="starfish-create-neighbourhood",
                                label="Neighbourhood",
                                value=10000,
                                min=0,
                                w=140,
                            ),
                        ],
                        grow=True,
                    ),
                    dmc.Textarea(
                        id="starfish-create-genomes-csv",
                        label="Genomes (CSV, paste)",
                        description=f"Header row required: {','.join(_SAMPLESHEET_COLS)}",
                        placeholder="genomeID,taxID,fna,gff3,emapper,cds,faa\ngenome1,5599,/data/genome1.fna,/data/genome1.gff3,,,",
                        autosize=True,
                        minRows=4,
                        maxRows=10,
                    ),
                    dmc.Group(
                        [
                            dmc.Button(
                                "Cancel",
                                id="starfish-create-cancel-btn",
                                variant="subtle",
                            ),
                            dmc.Button(
                                "Create", id="starfish-create-submit-btn", color="blue"
                            ),
                        ],
                        justify="flex-end",
                    ),
                ],
                gap="sm",
            ),
        ],
    )


def _build_list_section():
    return dag.AgGrid(
        id="starfish-runs-grid",
        columnDefs=[
            {"field": "id", "headerName": "ID", "width": 80, "pinned": "left"},
            {"field": "run_name", "headerName": "Run name", "flex": 2},
            {"field": "status_label", "headerName": "Status", "width": 120},
            {"field": "num_genomes", "headerName": "Genomes", "width": 100},
            {"field": "num_elements_found", "headerName": "Elements", "width": 100},
            {"field": "created_at", "headerName": "Created", "width": 170},
            {"field": "duration_display", "headerName": "Duration", "width": 100},
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


def _action_buttons(status):
    buttons = []
    if status == "pending":
        buttons.append(dmc.Button("Start", id="starfish-start-btn", color="green"))
    elif status == "running":
        buttons.append(dmc.Button("Cancel", id="starfish-cancel-btn", color="red"))
    elif status == "completed":
        buttons.append(
            dmc.Button(
                "Import all found elements", id="starfish-import-btn", color="blue"
            )
        )
        buttons.append(
            dmc.Button(
                "Re-run", id="starfish-rerun-btn", variant="outline", color="orange"
            )
        )
    elif status in ("failed", "cancelled"):
        buttons.append(dmc.Button("Resume", id="starfish-resume-btn", color="teal"))
        buttons.append(
            dmc.Button(
                "Re-run", id="starfish-rerun-btn", variant="outline", color="orange"
            )
        )
    return dmc.Group(buttons, gap="xs")


def _build_detail_section():
    return dmc.Stack(
        [
            dmc.Group(
                [
                    dmc.Button(
                        "← Back to runs",
                        id="starfish-back-btn",
                        variant="subtle",
                        size="sm",
                    ),
                ],
                justify="space-between",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Group(
                            [
                                dmc.Text(id="starfish-detail-name", fw=700, size="lg"),
                                html.Div(id="starfish-detail-status-badge"),
                            ],
                            gap="sm",
                        ),
                        dmc.Text(
                            id="starfish-detail-description", size="sm", c="dimmed"
                        ),
                        html.Div(id="starfish-detail-actions"),
                        html.Div(id="starfish-detail-error"),
                        dmc.Group(
                            [
                                dmc.Text(
                                    id="starfish-detail-timing", size="sm", c="dimmed"
                                ),
                            ]
                        ),
                        dmc.SimpleGrid(
                            [
                                dmc.Paper(
                                    dmc.Stack(
                                        [
                                            dmc.Text(
                                                id="starfish-stat-genomes",
                                                size="xl",
                                                fw=700,
                                                ta="center",
                                            ),
                                            dmc.Text(
                                                "Genomes",
                                                size="xs",
                                                c="dimmed",
                                                ta="center",
                                            ),
                                        ],
                                        gap=0,
                                    ),
                                    p="sm",
                                    withBorder=True,
                                ),
                                dmc.Paper(
                                    dmc.Stack(
                                        [
                                            dmc.Text(
                                                id="starfish-stat-elements",
                                                size="xl",
                                                fw=700,
                                                ta="center",
                                            ),
                                            dmc.Text(
                                                "Elements found",
                                                size="xs",
                                                c="dimmed",
                                                ta="center",
                                            ),
                                        ],
                                        gap=0,
                                    ),
                                    p="sm",
                                    withBorder=True,
                                ),
                            ],
                            cols=2,
                        ),
                        dmc.Text(id="starfish-detail-config", size="sm", c="dimmed"),
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
                            dmc.AccordionControl("Genomes"),
                            dmc.AccordionPanel(html.Div(id="starfish-detail-genomes")),
                        ],
                        value="genomes",
                    ),
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Elements"),
                            dmc.AccordionPanel(html.Div(id="starfish-detail-elements")),
                        ],
                        value="elements",
                    ),
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Logs"),
                            dmc.AccordionPanel(
                                html.Pre(
                                    id="starfish-detail-logs",
                                    style={
                                        "maxHeight": "400px",
                                        "overflowY": "auto",
                                        "backgroundColor": "#222",
                                        "color": "#eee",
                                        "padding": "0.75rem",
                                        "fontSize": "0.8rem",
                                    },
                                )
                            ),
                        ],
                        value="logs",
                    ),
                ],
                multiple=True,
                value=["genomes"],
            ),
            dcc.Interval(id="starfish-refresh-interval", interval=5000, disabled=True),
        ],
        gap="md",
    )


layout = html.Div(
    [
        dcc.Location(id="starfish-url", refresh=False),
        dcc.Store(id="starfish-selected-run-id", data=None),
        dcc.Store(id="starfish-authed", data=False),
        dcc.Store(id="starfish-list-refresh", data=0),
        html.Div(id="starfish-content"),
    ]
)


# ---------------------------------------------------------------------------
# Callbacks — token gate
# ---------------------------------------------------------------------------


@callback(
    [Output("starfish-content", "children"), Output("starfish-authed", "data")],
    Input("starfish-url", "search"),
)
def render_starfish_content(search):
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
    return _build_starfish_layout(), True


# ---------------------------------------------------------------------------
# Callbacks — list
# ---------------------------------------------------------------------------


@callback(
    Output("starfish-list-section", "children"),
    Input("starfish-authed", "data"),
)
def init_list_section(authed):
    if not authed:
        raise PreventUpdate
    return _build_list_section()


@callback(
    Output("starfish-runs-grid", "rowData"),
    Input("starfish-authed", "data"),
    Input("starfish-list-refresh", "data"),
    prevent_initial_call=True,
)
def refresh_runs_list(authed, _refresh):
    if not authed:
        raise PreventUpdate
    return _list_rows()


clientside_callback(
    """
    function(selectedRows) {
        if (!selectedRows || selectedRows.length === 0) {
            return window.dash_clientside.no_update;
        }
        return selectedRows[0].id;
    }
    """,
    Output("starfish-selected-run-id", "data"),
    Input("starfish-runs-grid", "selectedRows"),
    prevent_initial_call=True,
)


@callback(
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-back-btn", "n_clicks"),
    prevent_initial_call=True,
)
def go_back_to_list(n_clicks):
    if not n_clicks:
        raise PreventUpdate
    return None


@callback(
    [
        Output("starfish-list-section", "style"),
        Output("starfish-detail-section", "style"),
        Output("starfish-detail-section", "children"),
    ],
    Input("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def toggle_sections(selected_id):
    if selected_id is None:
        return {"display": "block"}, {"display": "none"}, no_update
    return {"display": "none"}, {"display": "block"}, _build_detail_section()


# ---------------------------------------------------------------------------
# Callbacks — create run
# ---------------------------------------------------------------------------


@callback(
    Output("starfish-create-modal", "opened"),
    Input("starfish-new-run-btn", "n_clicks"),
    Input("starfish-create-cancel-btn", "n_clicks"),
    prevent_initial_call=True,
)
def toggle_create_modal(open_clicks, cancel_clicks):
    trigger = dash.ctx.triggered_id
    return trigger == "starfish-new-run-btn"


@callback(
    Output("starfish-create-modal", "opened", allow_duplicate=True),
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-list-refresh", "data", allow_duplicate=True),
    Input("starfish-create-submit-btn", "n_clicks"),
    State("starfish-create-run-name", "value"),
    State("starfish-create-description", "value"),
    State("starfish-create-model", "value"),
    State("starfish-create-threads", "value"),
    State("starfish-create-missing", "value"),
    State("starfish-create-maxcopy", "value"),
    State("starfish-create-pid", "value"),
    State("starfish-create-hsp", "value"),
    State("starfish-create-flank", "value"),
    State("starfish-create-neighbourhood", "value"),
    State("starfish-create-genomes-csv", "value"),
    State("starfish-list-refresh", "data"),
    prevent_initial_call=True,
)
def submit_create_run(
    n_clicks,
    run_name,
    description,
    model,
    threads,
    missing,
    maxcopy,
    pid_threshold,
    hsp,
    flank,
    neighbourhood,
    genomes_csv,
    refresh_count,
):
    if not n_clicks:
        raise PreventUpdate
    try:
        genomes = _parse_genomes_csv(genomes_csv)
        run = starfish_manager.create_run(
            run_name,
            genomes,
            description=description or "",
            model=model,
            threads=threads,
            missing=missing,
            maxcopy=maxcopy,
            pid_threshold=pid_threshold,
            hsp=hsp,
            flank=flank,
            neighbourhood=neighbourhood,
        )
        notif = dmc.Notification(
            id=f"starfish-create-ok-{time.time()}",
            title="Run created",
            message=f"{run['run_name']} ({run['num_genomes']} genome(s))",
            color="green",
            action="show",
            autoClose=5000,
        )
        return False, notif, (refresh_count or 0) + 1
    except ValueError as exc:
        notif = dmc.Notification(
            id=f"starfish-create-err-{time.time()}",
            title="Could not create run",
            message=str(exc),
            color="red",
            action="show",
            autoClose=8000,
        )
        return True, notif, no_update


# ---------------------------------------------------------------------------
# Callbacks — detail load
# ---------------------------------------------------------------------------


@callback(
    [
        Output("starfish-detail-name", "children"),
        Output("starfish-detail-status-badge", "children"),
        Output("starfish-detail-description", "children"),
        Output("starfish-detail-actions", "children"),
        Output("starfish-detail-error", "children"),
        Output("starfish-detail-timing", "children"),
        Output("starfish-stat-genomes", "children"),
        Output("starfish-stat-elements", "children"),
        Output("starfish-detail-config", "children"),
        Output("starfish-detail-genomes", "children"),
        Output("starfish-detail-elements", "children"),
        Output("starfish-detail-logs", "children"),
        Output("starfish-refresh-interval", "disabled"),
    ],
    Input("starfish-selected-run-id", "data"),
    Input("starfish-detail-section", "children"),
    Input("starfish-refresh-interval", "n_intervals"),
    prevent_initial_call=True,
)
def load_detail(selected_id, _section_children, _n_intervals):
    if selected_id is None:
        raise PreventUpdate
    run = starfish_manager.get_run(selected_id)

    timing_parts = [f"Created {run['created_at']}"]
    if run.get("started_at"):
        timing_parts.append(f"Started {run['started_at']}")
    if run.get("completed_at"):
        timing_parts.append(f"Completed {run['completed_at']}")
    if run.get("duration_seconds"):
        timing_parts.append(f"Duration {run['duration_seconds']:.0f}s")

    error_children = (
        dmc.Alert(run["error_message"], title="Error", color="red", mt="sm")
        if run.get("error_message")
        else None
    )

    config = (
        f"model={run['model']}  threads={run['threads']}  missing={run['missing']}  "
        f"maxcopy={run['maxcopy']}  pid={run['pid_threshold']}%  hsp={run['hsp']}  "
        f"flank={run['flank']}  neighbourhood={run['neighbourhood']}"
    )

    try:
        log_content = (
            starfish_manager.get_run_log(selected_id) or "No log file available yet."
        )
    except Exception as exc:
        log_content = f"Error loading log: {exc}"

    return (
        run["run_name"],
        _status_badge(run["status"]),
        run.get("description") or "",
        _action_buttons(run["status"]),
        error_children,
        "  •  ".join(timing_parts),
        str(run.get("num_genomes") or 0),
        str(run.get("num_elements_found") or 0),
        config,
        _genomes_table(run["genomes"]),
        _elements_table(run["elements"]),
        log_content,
        run["status"] != "running",
    )


# ---------------------------------------------------------------------------
# Callbacks — run actions
# ---------------------------------------------------------------------------


def _run_action(action_fn, run_id, ok_title):
    try:
        action_fn(run_id)
        return dmc.Notification(
            id=f"starfish-action-ok-{time.time()}",
            title=ok_title,
            message="",
            color="green",
            action="show",
            autoClose=4000,
        )
    except ValueError as exc:
        return dmc.Notification(
            id=f"starfish-action-err-{time.time()}",
            title="Action failed",
            message=str(exc),
            color="red",
            action="show",
            autoClose=8000,
        )


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-start-btn", "n_clicks"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def start_run(n_clicks, run_id):
    if not n_clicks:
        raise PreventUpdate
    notif = _run_action(starfish_manager.start_run, run_id, "Run started")
    return notif, run_id


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-cancel-btn", "n_clicks"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def cancel_run(n_clicks, run_id):
    if not n_clicks:
        raise PreventUpdate
    notif = _run_action(starfish_manager.cancel_run, run_id, "Run cancelled")
    return notif, run_id


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-rerun-btn", "n_clicks"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def rerun_run(n_clicks, run_id):
    if not n_clicks:
        raise PreventUpdate
    notif = _run_action(starfish_manager.rerun_run, run_id, "Run restarted")
    return notif, run_id


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-resume-btn", "n_clicks"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def resume_run(n_clicks, run_id):
    if not n_clicks:
        raise PreventUpdate
    notif = _run_action(starfish_manager.resume_run, run_id, "Run resumed")
    return notif, run_id


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("starfish-selected-run-id", "data", allow_duplicate=True),
    Input("starfish-import-btn", "n_clicks"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def import_all_elements(n_clicks, run_id):
    if not n_clicks:
        raise PreventUpdate
    run = starfish_manager.get_run(run_id)
    to_import = [e for e in run["elements"] if not e.get("imported_submission_id")]
    imported, errors = 0, 0
    for e in to_import:
        try:
            starfish_manager.import_element_to_submission(e["id"])
            imported += 1
        except ValueError as exc:
            logger.warning("Import failed for element %s: %s", e["id"], exc)
            errors += 1
    color = "red" if errors and not imported else ("yellow" if errors else "green")
    notif = dmc.Notification(
        id=f"starfish-import-all-{time.time()}",
        title="Import complete",
        message=f"{imported} element(s) imported to the submissions queue"
        + (f", {errors} error(s)" if errors else "")
        + (" -- none were unimported" if not to_import else ""),
        color=color,
        action="show",
        autoClose=6000,
    )
    return notif, run_id

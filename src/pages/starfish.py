from urllib.parse import parse_qs
import base64
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


def _elements_grid(elements):
    rows = [
        {
            **e,
            "position_display": f"{e['start']}-{e['end']}",
            "import_display": "imported"
            if e.get("imported_submission_id")
            else "not imported",
        }
        for e in elements
    ]
    return dag.AgGrid(
        id="starfish-elements-grid",
        columnDefs=[
            {"field": "element_id", "headerName": "Element ID", "flex": 1},
            {"field": "genome_id", "headerName": "Genome", "width": 110},
            {"field": "contig_id", "headerName": "Contig", "width": 110},
            {"field": "position_display", "headerName": "Position", "width": 130},
            {"field": "length", "headerName": "Length", "width": 80},
            {"field": "family", "headerName": "Family", "width": 110},
            {"field": "navis", "headerName": "Navis", "width": 110},
            {"field": "haplotype", "headerName": "Haplotype", "width": 110},
            {"field": "confidence", "headerName": "Confidence", "width": 100},
            {"field": "import_display", "headerName": "Import", "width": 100},
        ],
        rowData=rows,
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": 25,
            "rowSelection": "single",
            "rowHeight": 40,
        },
        style={"height": "40vh"},
        className="ag-theme-alpine",
    )


def _element_actions(has_elements):
    if not has_elements:
        return html.Div()
    return dmc.Group(
        [
            dmc.Button(
                "Import selected",
                id="starfish-element-import-btn",
                color="blue",
                size="sm",
            ),
            dmc.Button(
                "Edit selected",
                id="starfish-element-edit-btn",
                variant="light",
                size="sm",
            ),
            dmc.Button(
                "Delete selected",
                id="starfish-element-delete-btn",
                variant="light",
                color="red",
                size="sm",
            ),
        ],
        gap="xs",
    )


def _viz_file_panel(title, section, file_list):
    if not file_list:
        return dmc.Paper(
            dmc.Text(
                f"No {title.lower()} outputs yet -- they appear once the pipeline "
                "reaches the viz steps.",
                size="sm",
                c="dimmed",
            ),
            p="sm",
            withBorder=True,
        )
    buttons = [
        dmc.Button(
            f,
            id={"type": "starfish-viz-open", "section": section, "index": i},
            variant="subtle",
            size="xs",
        )
        for i, f in enumerate(file_list)
    ]
    return dmc.Paper(
        dmc.Stack(
            [
                dmc.Text(title, fw=700, size="sm"),
                dmc.Group(buttons, wrap=True, gap="xs"),
            ],
            gap="xs",
        ),
        p="sm",
        withBorder=True,
    )


def _viz_section(files):
    return dmc.SimpleGrid(
        [
            _viz_file_panel(
                "Locus visualizations", "locusViz", files.get("locusViz_files") or []
            ),
            _viz_file_panel(
                "Pair visualizations", "pairViz", files.get("pairViz_files") or []
            ),
        ],
        cols=2,
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
            _viz_modal(),
            _edit_element_modal(),
            _delete_element_modal(),
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


def _viz_modal():
    return dmc.Modal(
        id="starfish-viz-modal",
        title="Visualization",
        size="lg",
        withCloseButton=True,
        children=[
            dmc.Stack(
                [
                    dmc.Text(id="starfish-viz-name", fw=700, size="sm"),
                    dmc.Image(
                        id="starfish-viz-image",
                        style={"maxHeight": "70vh", "display": "none"},
                    ),
                    dmc.Group(
                        [
                            dmc.Button(
                                "Download",
                                id="starfish-viz-download-btn",
                                size="sm",
                            ),
                        ],
                        justify="flex-end",
                    ),
                ],
                gap="sm",
            ),
        ],
    )


def _edit_element_modal():
    return dmc.Modal(
        id="starfish-edit-modal",
        title="Edit element",
        size="md",
        withCloseButton=True,
        children=[
            dmc.Stack(
                [
                    dmc.Text(id="starfish-edit-element-label", size="sm", c="dimmed"),
                    dmc.TextInput(id="starfish-edit-family", label="Family"),
                    dmc.TextInput(id="starfish-edit-navis", label="Navis"),
                    dmc.TextInput(id="starfish-edit-haplotype", label="Haplotype"),
                    dmc.TextInput(id="starfish-edit-confidence", label="Confidence"),
                    dmc.Textarea(
                        id="starfish-edit-notes",
                        label="Notes",
                        autosize=True,
                        minRows=2,
                        maxRows=6,
                    ),
                    dmc.Group(
                        [
                            dmc.Button(
                                "Cancel",
                                id="starfish-edit-cancel-btn",
                                variant="subtle",
                            ),
                            dmc.Button(
                                "Save", id="starfish-edit-save-btn", color="blue"
                            ),
                        ],
                        justify="flex-end",
                    ),
                ],
                gap="sm",
            ),
        ],
    )


def _delete_element_modal():
    return dmc.Modal(
        id="starfish-delete-modal",
        title="Delete element?",
        size="sm",
        withCloseButton=True,
        children=[
            dmc.Stack(
                [
                    dmc.Text(id="starfish-delete-element-label", size="sm"),
                    dmc.Text(
                        "Removes the element row from this run. Imported elements "
                        "cannot be deleted.",
                        size="sm",
                        c="dimmed",
                    ),
                    dmc.Group(
                        [
                            dmc.Button(
                                "Cancel",
                                id="starfish-delete-cancel-btn",
                                variant="subtle",
                            ),
                            dmc.Button(
                                "Delete",
                                id="starfish-delete-confirm-btn",
                                color="red",
                            ),
                        ],
                        justify="flex-end",
                    ),
                ],
                gap="sm",
            ),
        ],
    )


def _build_list_section(rows=None):
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
        rowData=rows if rows is not None else [],
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
                            dmc.AccordionPanel(
                                dmc.Stack(
                                    [
                                        html.Div(id="starfish-detail-elements-hint"),
                                        html.Div(id="starfish-detail-elements-actions"),
                                        html.Div(id="starfish-detail-elements"),
                                    ],
                                    gap="xs",
                                )
                            ),
                        ],
                        value="elements",
                    ),
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl("Visualizations"),
                            dmc.AccordionPanel(html.Div(id="starfish-detail-viz")),
                        ],
                        value="viz",
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
                id="starfish-detail-accordion",
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
        dcc.Store(id="starfish-detail-refresh", data=0),
        dcc.Store(
            id="starfish-viz-files", data={"locusViz_files": [], "pairViz_files": []}
        ),
        dcc.Store(id="starfish-viz-payload", data=None),
        dcc.Store(id="starfish-edit-element-id", data=None),
        dcc.Store(id="starfish-delete-element-id", data=None),
        dcc.Download(id="starfish-viz-download"),
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
    return _build_list_section(_list_rows())


@callback(
    Output("starfish-runs-grid", "rowData"),
    Input("starfish-list-refresh", "data"),
    prevent_initial_call=True,
)
def refresh_runs_list(_refresh):
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
        Output("starfish-detail-elements-hint", "children"),
        Output("starfish-detail-elements-actions", "children"),
        Output("starfish-detail-elements", "children"),
        Output("starfish-detail-logs", "children"),
        Output("starfish-refresh-interval", "disabled"),
    ],
    Input("starfish-selected-run-id", "data"),
    Input("starfish-detail-section", "children"),
    Input("starfish-refresh-interval", "n_intervals"),
    Input("starfish-detail-refresh", "data"),
    prevent_initial_call=True,
)
def load_detail(selected_id, _section_children, _n_intervals, _detail_refresh):
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

    elements = run["elements"]
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
        dmc.Text(
            "No elements found yet -- run the pipeline to detect starship elements.",
            size="sm",
            c="dimmed",
        )
        if not elements
        else "",
        _element_actions(bool(elements)),
        _elements_grid(elements),
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


# ---------------------------------------------------------------------------
# Callbacks — per-element actions (import / edit / delete)
# ---------------------------------------------------------------------------


def _notif(title, message, color):
    return dmc.Notification(
        id=f"starfish-element-{color}-{time.time()}",
        title=title,
        message=message,
        color=color,
        action="show",
        autoClose=6000,
    )


@callback(
    [
        Output("notifications-container", "children", allow_duplicate=True),
        Output("starfish-detail-refresh", "data", allow_duplicate=True),
    ],
    Input("starfish-element-import-btn", "n_clicks"),
    State("starfish-elements-grid", "selectedRows"),
    State("starfish-detail-refresh", "data"),
    prevent_initial_call=True,
)
def import_selected_element(n_clicks, selected_rows, refresh_count):
    if not n_clicks:
        raise PreventUpdate
    if not selected_rows:
        raise PreventUpdate
    row = selected_rows[0]
    try:
        result = starfish_manager.import_element_to_submission(row["id"])
        notif = _notif(
            "Element imported",
            f"{result['element_id']} added to the submissions queue as "
            f"submission {result['submission_id']}",
            "green",
        )
    except ValueError as exc:
        notif = _notif("Import failed", str(exc), "red")
    return notif, (refresh_count or 0) + 1


@callback(
    [
        Output("starfish-edit-modal", "opened"),
        Output("starfish-edit-element-id", "data"),
        Output("starfish-edit-element-label", "children"),
        Output("starfish-edit-family", "value"),
        Output("starfish-edit-navis", "value"),
        Output("starfish-edit-haplotype", "value"),
        Output("starfish-edit-confidence", "value"),
        Output("starfish-edit-notes", "value"),
    ],
    Input("starfish-element-edit-btn", "n_clicks"),
    Input("starfish-edit-cancel-btn", "n_clicks"),
    State("starfish-elements-grid", "selectedRows"),
    prevent_initial_call=True,
)
def open_edit_modal(edit_clicks, cancel_clicks, selected_rows):
    if dash.ctx.triggered_id == "starfish-edit-cancel-btn":
        return (
            False,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
        )
    if not edit_clicks:
        raise PreventUpdate
    if not selected_rows:
        raise PreventUpdate
    row = selected_rows[0]
    return (
        True,
        row["id"],
        f"{row['element_id']}  ({row.get('contig_id') or '?'}:{row['start']}-{row['end']})",
        row.get("family") or "",
        row.get("navis") or "",
        row.get("haplotype") or "",
        row.get("confidence") or "",
        row.get("notes") or "",
    )


@callback(
    [
        Output("starfish-edit-modal", "opened", allow_duplicate=True),
        Output("notifications-container", "children", allow_duplicate=True),
        Output("starfish-detail-refresh", "data", allow_duplicate=True),
    ],
    Input("starfish-edit-save-btn", "n_clicks"),
    State("starfish-edit-element-id", "data"),
    State("starfish-edit-family", "value"),
    State("starfish-edit-navis", "value"),
    State("starfish-edit-haplotype", "value"),
    State("starfish-edit-confidence", "value"),
    State("starfish-edit-notes", "value"),
    State("starfish-detail-refresh", "data"),
    prevent_initial_call=True,
)
def save_edit(
    n_clicks,
    element_id,
    family,
    navis,
    haplotype,
    confidence,
    notes,
    refresh_count,
):
    if not n_clicks:
        raise PreventUpdate
    if element_id is None:
        raise PreventUpdate
    try:
        updated = starfish_manager.update_element(
            element_id,
            family=family,
            navis=navis,
            haplotype=haplotype,
            confidence=confidence,
            notes=notes,
        )
        notif = _notif(
            "Element updated",
            f"{updated['element_id']} saved",
            "green",
        )
    except ValueError as exc:
        notif = _notif("Could not save element", str(exc), "red")
    return False, notif, (refresh_count or 0) + 1


@callback(
    [
        Output("starfish-delete-modal", "opened"),
        Output("starfish-delete-element-id", "data"),
        Output("starfish-delete-element-label", "children"),
    ],
    Input("starfish-element-delete-btn", "n_clicks"),
    Input("starfish-delete-cancel-btn", "n_clicks"),
    State("starfish-elements-grid", "selectedRows"),
    prevent_initial_call=True,
)
def open_delete_modal(delete_clicks, cancel_clicks, selected_rows):
    if dash.ctx.triggered_id == "starfish-delete-cancel-btn":
        return False, no_update, no_update
    if not delete_clicks:
        raise PreventUpdate
    if not selected_rows:
        raise PreventUpdate
    row = selected_rows[0]
    return (
        True,
        row["id"],
        f"{row['element_id']}  ({row.get('contig_id') or '?'}:{row['start']}-{row['end']})",
    )


@callback(
    [
        Output("starfish-delete-modal", "opened", allow_duplicate=True),
        Output("notifications-container", "children", allow_duplicate=True),
        Output("starfish-detail-refresh", "data", allow_duplicate=True),
    ],
    Input("starfish-delete-confirm-btn", "n_clicks"),
    State("starfish-delete-element-id", "data"),
    State("starfish-detail-refresh", "data"),
    prevent_initial_call=True,
)
def confirm_delete(n_clicks, element_id, refresh_count):
    if not n_clicks:
        raise PreventUpdate
    if element_id is None:
        raise PreventUpdate
    try:
        result = starfish_manager.delete_element(element_id)
        notif = _notif(
            "Element deleted",
            f"{result['deleted']} removed from this run",
            "green",
        )
    except ValueError as exc:
        notif = _notif("Could not delete element", str(exc), "red")
    return False, notif, (refresh_count or 0) + 1


# ---------------------------------------------------------------------------
# Callbacks — visualizations
# ---------------------------------------------------------------------------


@callback(
    [
        Output("starfish-detail-viz", "children"),
        Output("starfish-viz-files", "data"),
    ],
    Input("starfish-selected-run-id", "data"),
    Input("starfish-detail-accordion", "value"),
    Input("starfish-detail-refresh", "data"),
    prevent_initial_call=True,
)
def load_visualizations(selected_id, accordion_value, _detail_refresh):
    if selected_id is None or "viz" not in (accordion_value or []):
        raise PreventUpdate
    try:
        files = starfish_manager.list_visualizations(selected_id)
    except ValueError:
        return no_update, no_update
    return _viz_section(files), files


@callback(
    [
        Output("starfish-viz-modal", "opened"),
        Output("starfish-viz-name", "children"),
        Output("starfish-viz-image", "src"),
        Output("starfish-viz-image", "style"),
        Output("starfish-viz-payload", "data"),
    ],
    Input({"type": "starfish-viz-open"}, "n_clicks"),
    State({"type": "starfish-viz-open"}, "id"),
    State("starfish-viz-files", "data"),
    State("starfish-selected-run-id", "data"),
    prevent_initial_call=True,
)
def open_viz_modal(_n_clicks, button_id, files, run_id):
    if not isinstance(button_id, dict) or "index" not in button_id:
        raise PreventUpdate
    if not run_id:
        raise PreventUpdate
    section = button_id["section"]
    index = button_id["index"]
    file_list = (files or {}).get(f"{section}_files") or []
    if index >= len(file_list):
        raise PreventUpdate
    filename = file_list[index]
    try:
        data, media_type = starfish_manager.get_visualization_file(
            run_id, section, filename
        )
    except ValueError as exc:
        logger.warning("Could not load viz file %s: %s", filename, exc)
        raise PreventUpdate
    b64 = base64.b64encode(data).decode()
    is_image = media_type.startswith("image/")
    return (
        True,
        filename,
        f"data:{media_type};base64,{b64}" if is_image else None,
        None if is_image else {"display": "none"},
        {"content": b64, "filename": filename, "type": "base64"},
    )


@callback(
    Output("starfish-viz-download", "data"),
    Input("starfish-viz-download-btn", "n_clicks"),
    State("starfish-viz-payload", "data"),
    prevent_initial_call=True,
)
def trigger_viz_download(n_clicks, payload):
    if not n_clicks or not payload:
        raise PreventUpdate
    return payload

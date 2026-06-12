import uuid
from datetime import datetime
from urllib.parse import parse_qs

import dash
import dash_ag_grid as dag
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import pandas as pd
from dash import Input, Output, State, callback, dcc, html, no_update
from sqlalchemy import text

from src.config.logging import get_logger
from src.config.settings import ADMIN_TOKEN
from src.database.sql_engine import get_starbase_session, get_submissions_session
from src.database.sql_manager import get_database_version, set_database_version

logger = get_logger(__name__)

dash.register_page(__name__, path="/admin", title="Admin", name="Admin")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EDITABLE_COLS = {
    "joined_ships": {"curated_status", "evidence", "source"},
    "submissions": {
        "needs_review",
        "classification_family",
        "classification_navis",
        "classification_haplotype",
        "classification_confidence",
        "comment",
    },
    "taxonomy": {"name", "taxID", "genus", "species", "strain"},
    "papers": {
        "Title",
        "Author",
        "PublicationYear",
        "DOI",
        "Url",
        "shortCitation",
        "starshipMentioned",
        "typePaper",
    },
}

TABLE_CONFIG = {
    "joined_ships": {"sql_table": "joined_ships", "db": "starbase", "pk": "id"},
    "submissions": {"sql_table": "submissions", "db": "submissions", "pk": "id"},
    "taxonomy": {"sql_table": "taxonomy", "db": "starbase", "pk": "id"},
    "papers": {"sql_table": "papers", "db": "starbase", "pk": "id"},
}

# Columns that should be stored as integers (SQLite booleans)
BOOLEAN_COLS = {"needs_review"}

# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def _fetch_joined_ships():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT j.id, j.starshipID, j.curated_status, j.evidence, j.source,
                   a.accession_tag,
                   f.familyName,
                   n.navis_name,
                   h.haplotype_name,
                   t.name AS taxonomy_name, t.genus, t.species
            FROM joined_ships j
            LEFT JOIN accessions a ON j.accession_id = a.id
            LEFT JOIN family_names f ON j.ship_family_id = f.id
            LEFT JOIN navis_names n ON j.ship_navis_id = n.id
            LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
            LEFT JOIN taxonomy t ON j.tax_id = t.id
            ORDER BY j.id DESC
            """,
            session.bind,
        )


def _fetch_submissions():
    with get_submissions_session() as session:
        return pd.read_sql_query(
            """
            SELECT id, seq_filename, uploader, seq_date,
                   needs_review, comment,
                   classification_family, classification_navis,
                   classification_haplotype, classification_confidence
            FROM submissions
            ORDER BY id DESC
            """,
            session.bind,
        )


def _fetch_taxonomy():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT id, name, taxID, genus, species, strain, kingdom, phylum
            FROM taxonomy
            ORDER BY id
            """,
            session.bind,
        )


def _fetch_papers():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT id, Title, Author, PublicationYear, DOI, Url,
                   shortCitation, starshipMentioned, typePaper
            FROM papers
            ORDER BY id
            """,
            session.bind,
        )


def _bump_version(version_str, bump_type):
    """Increment semantic version string by minor or patch."""
    try:
        parts = [int(x) for x in str(version_str).split(".")]
        while len(parts) < 3:
            parts.append(0)
        major, minor, patch = parts[0], parts[1], parts[2]
        if bump_type == "minor":
            return f"{major}.{minor + 1}.0"
        if bump_type == "patch":
            return f"{major}.{minor}.{patch + 1}"
    except Exception:
        pass
    return version_str


def _do_update(table_key, row_id, col_id, new_value):
    """Run a whitelisted UPDATE. Returns (success: bool, error: str|None)."""
    allowed = EDITABLE_COLS.get(table_key, set())
    if col_id not in allowed:
        return False, f"Column '{col_id}' is not editable."

    config = TABLE_CONFIG[table_key]

    if col_id in BOOLEAN_COLS:
        new_value = 1 if str(new_value).lower() in ("true", "1", "yes") else 0

    sql = text(
        f"UPDATE {config['sql_table']} SET {col_id} = :val WHERE {config['pk']} = :pk"
    )

    try:
        ctx = (
            get_starbase_session
            if config["db"] == "starbase"
            else get_submissions_session
        )
        with ctx() as session:
            session.execute(sql, {"val": new_value, "pk": row_id})
            session.commit()
        logger.info(
            "admin UPDATE %s.%s id=%s  %r",
            config["sql_table"],
            col_id,
            row_id,
            new_value,
        )
        return True, None
    except Exception as exc:
        logger.error(
            "admin UPDATE error (%s.%s id=%s): %s", table_key, col_id, row_id, exc
        )
        return False, str(exc)


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------


def _make_grid(df, grid_id, editable_cols):
    col_defs = []
    for col in df.columns:
        is_editable = col in editable_cols
        col_def = {
            "field": col,
            "headerName": col,
            "editable": is_editable,
            "filter": True,
            "sortable": True,
            "resizable": True,
            "minWidth": 90,
        }
        if col == "id":
            col_def.update(
                {
                    "width": 72,
                    "pinned": "left",
                    "editable": False,
                    "filter": "agNumberColumnFilter",
                }
            )
        elif is_editable:
            col_def["cellStyle"] = {"backgroundColor": "var(--mantine-color-yellow-0)"}
        col_defs.append(col_def)

    return dag.AgGrid(
        id=grid_id,
        columnDefs=col_defs,
        rowData=df.fillna("").to_dict("records"),
        defaultColDef={"resizable": True, "minWidth": 80},
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": 50,
            "suppressPropertyNamesCheck": True,
            "rowHeight": 40,
            "headerHeight": 44,
        },
        getRowId="params.data.id",
        className="ag-theme-alpine",
        style={"width": "100%", "height": "600px"},
    )


def _version_bump_modal():
    return dmc.Modal(
        id="admin-version-modal",
        title="Bump database content version?",
        centered=True,
        size="md",
        opened=False,
        children=dmc.Stack(
            [
                html.Div(id="admin-version-modal-info"),
                dmc.SegmentedControl(
                    id="admin-version-bump-type",
                    value="minor",
                    data=[
                        {"label": "Minor — new/updated content", "value": "minor"},
                        {"label": "Patch — small correction", "value": "patch"},
                        {"label": "Skip", "value": "skip"},
                    ],
                    fullWidth=True,
                ),
                dmc.Textarea(
                    id="admin-version-description",
                    label="Description (recorded with the version entry)",
                    placeholder="e.g. corrected curated_status for SSB123 via admin",
                    autosize=True,
                    minRows=2,
                ),
                dmc.Group(
                    [
                        dmc.Button("Confirm", id="admin-version-confirm", color="blue"),
                        dmc.Button(
                            "Cancel",
                            id="admin-version-cancel",
                            variant="subtle",
                            color="gray",
                        ),
                    ],
                    justify="flex-end",
                ),
            ],
            gap="md",
        ),
    )


def _build_admin_layout():
    try:
        js_df = _fetch_joined_ships()
        sub_df = _fetch_submissions()
        tax_df = _fetch_taxonomy()
        pap_df = _fetch_papers()
    except Exception as exc:
        logger.error("Admin data load failed: %s", exc)
        return dmc.Alert(str(exc), title="Database Error", color="red", mt="xl")

    current_version = get_database_version()

    return dmc.Container(
        [
            dmc.Group(
                [
                    dmc.Title("Admin Panel", order=2),
                    dmc.Badge(
                        f"Content v{current_version}",
                        color="blue",
                        variant="light",
                        size="lg",
                    ),
                ],
                justify="space-between",
                align="center",
                mt="lg",
                mb="xs",
            ),
            dmc.Text(
                "Yellow-highlighted cells are editable. "
                "Press Enter or click away to save. All changes are logged.",
                size="sm",
                c="dimmed",
                mb="md",
            ),
            dbc.Tabs(
                [
                    dbc.Tab(
                        _make_grid(
                            js_df,
                            "admin-joined-ships-grid",
                            EDITABLE_COLS["joined_ships"],
                        ),
                        label=f"Joined Ships ({len(js_df):,})",
                        tab_id="joined_ships",
                    ),
                    dbc.Tab(
                        _make_grid(
                            sub_df,
                            "admin-submissions-grid",
                            EDITABLE_COLS["submissions"],
                        ),
                        label=f"Submissions ({len(sub_df):,})",
                        tab_id="submissions",
                    ),
                    dbc.Tab(
                        _make_grid(
                            tax_df,
                            "admin-taxonomy-grid",
                            EDITABLE_COLS["taxonomy"],
                        ),
                        label=f"Taxonomy ({len(tax_df):,})",
                        tab_id="taxonomy",
                    ),
                    dbc.Tab(
                        _make_grid(
                            pap_df,
                            "admin-papers-grid",
                            EDITABLE_COLS["papers"],
                        ),
                        label=f"Papers ({len(pap_df):,})",
                        tab_id="papers",
                    ),
                ],
                active_tab="joined_ships",
            ),
        ],
        size="xl",
        style={"paddingBottom": "4rem"},
    )


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

layout = html.Div(
    [
        dcc.Location(id="admin-url", refresh=False),
        dcc.Store(id="admin-save-store"),
        html.Div(id="admin-content"),
        _version_bump_modal(),
    ]
)

# ---------------------------------------------------------------------------
# Callbacks — token gate
# ---------------------------------------------------------------------------


@callback(
    Output("admin-content", "children"),
    Input("admin-url", "search"),
)
def render_admin_content(search):
    params = parse_qs((search or "").lstrip("?"))
    token = params.get("token", [None])[0]
    if not ADMIN_TOKEN or not token or token != ADMIN_TOKEN:
        return dmc.Container(
            dmc.Alert(
                "A valid ?token= query parameter is required.",
                title="Unauthorized",
                color="red",
                mt="xl",
            ),
            size="sm",
        )
    return _build_admin_layout()


# ---------------------------------------------------------------------------
# Callbacks — cell edits → DB update → store
# ---------------------------------------------------------------------------


def _handle_cell_change(table_key, cell_changed):
    if not cell_changed:
        return no_update
    row_data = cell_changed.get("data", {})
    row_id = row_data.get("id")
    col_id = cell_changed.get("colId")
    new_value = cell_changed.get("value")
    old_value = cell_changed.get("oldValue")

    if row_id is None or col_id is None:
        return no_update

    success, error = _do_update(table_key, row_id, col_id, new_value)
    return {
        "table": table_key,
        "row_id": row_id,
        "col_id": col_id,
        "old_value": old_value,
        "new_value": new_value,
        "success": success,
        "error": error,
        "ts": datetime.now().isoformat(),
    }


@callback(
    Output("admin-save-store", "data"),
    Input("admin-joined-ships-grid", "cellValueChanged"),
    prevent_initial_call=True,
)
def save_joined_ships(cell_changed):
    return _handle_cell_change("joined_ships", cell_changed)


@callback(
    Output("admin-save-store", "data", allow_duplicate=True),
    Input("admin-submissions-grid", "cellValueChanged"),
    prevent_initial_call=True,
)
def save_submissions(cell_changed):
    return _handle_cell_change("submissions", cell_changed)


@callback(
    Output("admin-save-store", "data", allow_duplicate=True),
    Input("admin-taxonomy-grid", "cellValueChanged"),
    prevent_initial_call=True,
)
def save_taxonomy(cell_changed):
    return _handle_cell_change("taxonomy", cell_changed)


@callback(
    Output("admin-save-store", "data", allow_duplicate=True),
    Input("admin-papers-grid", "cellValueChanged"),
    prevent_initial_call=True,
)
def save_papers(cell_changed):
    return _handle_cell_change("papers", cell_changed)


# ---------------------------------------------------------------------------
# Callbacks — store → open version modal (or show error notification)
# ---------------------------------------------------------------------------


@callback(
    [
        Output("admin-version-modal", "opened"),
        Output("admin-version-modal-info", "children"),
        Output("admin-version-description", "value"),
        Output("notifications-container", "children"),
    ],
    Input("admin-save-store", "data"),
    prevent_initial_call=True,
)
def on_save_store_updated(store_data):
    if not store_data:
        return no_update, no_update, no_update, no_update

    if not store_data.get("success"):
        err = store_data.get("error", "Unknown error")
        notification = dmc.Notification(
            id=f"admin-err-{uuid.uuid4().hex[:6]}",
            title="Save failed",
            message=err,
            color="red",
            action="show",
            autoClose=7000,
        )
        return False, no_update, no_update, notification

    current_ver = get_database_version()
    table = store_data.get("table", "")
    col_id = store_data.get("col_id", "")
    row_id = store_data.get("row_id", "")
    old_val = store_data.get("old_value", "")
    new_val = store_data.get("new_value", "")

    minor_ver = _bump_version(current_ver, "minor")
    patch_ver = _bump_version(current_ver, "patch")

    modal_info = dmc.Stack(
        [
            dmc.Text(f"Saved: {table}.{col_id}  (row {row_id})", size="sm", fw=500),
            dmc.Code(f'"{old_val}"  →  "{new_val}"', block=False),
            dmc.Divider(),
            dmc.Group(
                [
                    dmc.Text("Current version:", size="sm"),
                    dmc.Badge(current_ver, variant="outline", color="gray"),
                ],
                gap="xs",
            ),
            dmc.Group(
                [
                    dmc.Text("Minor bump →", size="sm"),
                    dmc.Badge(minor_ver, variant="light", color="blue"),
                    dmc.Text("Patch bump →", size="sm"),
                    dmc.Badge(patch_ver, variant="light", color="teal"),
                ],
                gap="xs",
            ),
        ],
        gap="xs",
    )

    default_description = f"edited {table}.{col_id} row {row_id} via admin"
    return True, modal_info, default_description, no_update


# ---------------------------------------------------------------------------
# Callbacks — version modal confirm / cancel
# ---------------------------------------------------------------------------


@callback(
    [
        Output("admin-version-modal", "opened", allow_duplicate=True),
        Output("notifications-container", "children", allow_duplicate=True),
    ],
    Input("admin-version-confirm", "n_clicks"),
    [
        State("admin-version-bump-type", "value"),
        State("admin-version-description", "value"),
    ],
    prevent_initial_call=True,
)
def confirm_version_bump(n_clicks, bump_type, description):
    if not n_clicks:
        return no_update, no_update

    if bump_type == "skip":
        notification = dmc.Notification(
            id=f"admin-ok-{uuid.uuid4().hex[:6]}",
            title="Saved",
            message="Change saved. Version not bumped.",
            color="green",
            action="show",
            autoClose=4000,
        )
        return False, notification

    current_ver = get_database_version()
    new_ver = _bump_version(current_ver, bump_type)

    try:
        set_database_version(new_ver, description or "", created_by="admin")
        msg = f"Content version: {current_ver} → {new_ver}"
        color = "green"
    except Exception as exc:
        msg = f"Saved but version bump failed: {exc}"
        color = "orange"

    notification = dmc.Notification(
        id=f"admin-ok-{uuid.uuid4().hex[:6]}",
        title="Saved",
        message=msg,
        color=color,
        action="show",
        autoClose=5000,
    )
    return False, notification


@callback(
    Output("admin-version-modal", "opened", allow_duplicate=True),
    Input("admin-version-cancel", "n_clicks"),
    prevent_initial_call=True,
)
def cancel_version_bump(n_clicks):
    return False

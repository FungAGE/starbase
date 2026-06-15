import base64
import uuid
from io import StringIO
from urllib.parse import parse_qs

from Bio import SeqIO
import dash
import dash_ag_grid as dag
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import pandas as pd
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
    "family_names": {
        "familyName",
        "notes",
        "type_element_reference",
        "clade",
        "longFamilyID",
        "oldFamilyID",
        "otherFamilyID",
    },
    "navis_names": {"navis_name", "previous_navis_name", "activity"},
    "haplotype_names": {"haplotype_name", "previous_haplotype_name", "activity"},
    "accessions": {"version_tag"},
    "ship_accessions": {"ship_accession_display", "ship_version_tag"},
    "genomes": {
        "ome",
        "version",
        "genomeSource",
        "citation",
        "biosample",
        "acquisition_date",
        "assembly_accession",
    },
}

TABLE_CONFIG = {
    "joined_ships": {"sql_table": "joined_ships", "db": "starbase", "pk": "id"},
    "submissions": {"sql_table": "submissions", "db": "submissions", "pk": "id"},
    "taxonomy": {"sql_table": "taxonomy", "db": "starbase", "pk": "id"},
    "papers": {"sql_table": "papers", "db": "starbase", "pk": "id"},
    "family_names": {"sql_table": "family_names", "db": "starbase", "pk": "id"},
    "navis_names": {"sql_table": "navis_names", "db": "starbase", "pk": "id"},
    "haplotype_names": {"sql_table": "haplotype_names", "db": "starbase", "pk": "id"},
    "accessions": {"sql_table": "accessions", "db": "starbase", "pk": "id"},
    "ship_accessions": {"sql_table": "ship_accessions", "db": "starbase", "pk": "id"},
    "genomes": {"sql_table": "genomes", "db": "starbase", "pk": "id"},
}

# Columns that should be stored as integers (SQLite booleans)
BOOLEAN_COLS = {"needs_review"}

# Grid IDs keyed by table_key — used for batch save/discard
GRID_IDS = {
    "joined_ships": "admin-joined-ships-grid",
    "submissions": "admin-submissions-grid",
    "taxonomy": "admin-taxonomy-grid",
    "papers": "admin-papers-grid",
    "family_names": "admin-family-names-grid",
    "navis_names": "admin-navis-names-grid",
    "haplotype_names": "admin-haplotype-names-grid",
    "accessions": "admin-accessions-grid",
    "ship_accessions": "admin-ship-accessions-grid",
    "genomes": "admin-genomes-grid",
}

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


def _fetch_family_names():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT id, familyName, longFamilyID, oldFamilyID, otherFamilyID,
                   clade, newFamilyID, type_element_reference, notes
            FROM family_names
            ORDER BY id
            """,
            session.bind,
        )


def _fetch_navis_names():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT n.id, n.navis_name, n.previous_navis_name, n.activity,
                   f.familyName
            FROM navis_names n
            LEFT JOIN family_names f ON n.ship_family_id = f.id
            ORDER BY n.id
            """,
            session.bind,
        )


def _fetch_haplotype_names():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT h.id, h.haplotype_name, h.previous_haplotype_name, h.activity,
                   n.navis_name, f.familyName
            FROM haplotype_names h
            LEFT JOIN navis_names n ON h.navis_id = n.id
            LEFT JOIN family_names f ON h.ship_family_id = f.id
            ORDER BY h.id
            """,
            session.bind,
        )


def _fetch_accessions():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT a.id, a.accession_tag, a.version_tag
            FROM accessions a
            ORDER BY a.id
            """,
            session.bind,
        )


def _fetch_ship_accessions():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT sa.id, sa.ship_accession_tag, sa.ship_accession_display,
                   sa.ship_version_tag, sa.ship_id
            FROM ship_accessions sa
            ORDER BY sa.id
            """,
            session.bind,
        )


def _fetch_genomes():
    with get_starbase_session() as session:
        return pd.read_sql_query(
            """
            SELECT g.id, g.ome, g.version, g.genomeSource, g.citation,
                   g.biosample, g.acquisition_date, g.assembly_accession,
                   t.name AS taxonomy_name
            FROM genomes g
            LEFT JOIN taxonomy t ON g.taxonomy_id = t.id
            ORDER BY g.id
            """,
            session.bind,
        )


def _refetch_rowdata(table_key):
    """Re-fetch clean rowData (no _dirty) for a grid after save/discard."""
    fetchers = {
        "joined_ships": _fetch_joined_ships,
        "submissions": _fetch_submissions,
        "taxonomy": _fetch_taxonomy,
        "papers": _fetch_papers,
        "family_names": _fetch_family_names,
        "navis_names": _fetch_navis_names,
        "haplotype_names": _fetch_haplotype_names,
        "accessions": _fetch_accessions,
        "ship_accessions": _fetch_ship_accessions,
        "genomes": _fetch_genomes,
    }
    try:
        return fetchers[table_key]().fillna("").to_dict("records")
    except Exception as exc:
        logger.error("Re-fetch failed for %s: %s", table_key, exc)
        return no_update


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
            result = session.execute(sql, {"val": new_value, "pk": row_id})
            session.commit()
            rc = result.rowcount
        if rc == 0:
            logger.warning(
                "admin UPDATE matched 0 rows: %s.%s id=%r (pk type=%s, value=%r)",
                config["sql_table"],
                col_id,
                row_id,
                type(row_id).__name__,
                new_value,
            )
            return False, f"Row id={row_id!r} not found in {config['sql_table']}."
        logger.info(
            "admin UPDATE %s.%s id=%s  %r  (rowcount=%d)",
            config["sql_table"],
            col_id,
            row_id,
            new_value,
            rc,
        )
        return True, None
    except Exception as exc:
        logger.error(
            "admin UPDATE error (%s.%s id=%s): %s", table_key, col_id, row_id, exc
        )
        return False, str(exc)


# ---------------------------------------------------------------------------
# Promote helper
# ---------------------------------------------------------------------------


def _promote_submission(sub_id: int):
    """
    Load a staging submission and insert it into the main starbase DB via
    perform_database_insertion (SubmissionProcessor).

    Returns (success: bool, accession: str|None, error: str|None).
    """
    from src.utils.web_submission_adapter import perform_database_insertion
    from src.utils.submission_utils import check_sequence_duplicate

    try:
        with get_submissions_session() as session:
            row = session.execute(
                text("SELECT * FROM submissions WHERE id = :id"), {"id": sub_id}
            ).fetchone()

        if not row:
            return False, None, f"Submission {sub_id} not found"

        row = dict(row._mapping)

        if not row.get("seq_contents"):
            return False, None, "No sequence data stored in this submission"

        content_type, content_string = row["seq_contents"].split(",", 1)
        decoded = base64.b64decode(content_string).decode("utf-8")
        records = list(SeqIO.parse(StringIO(decoded), "fasta"))
        if not records:
            return False, None, "No valid FASTA sequences in submission"
        sequence = str(records[0].seq)
        seq_id = records[0].id

        duplicate_info = check_sequence_duplicate(
            sequence, row.get("genus") or "", row.get("species") or ""
        )

        strand = row.get("shipstrand") or "+"
        processed_data = {
            "sequence": sequence,
            "starshipID": seq_id,
            "evidence": row.get("evidence") or "",
            "source": "web_submission_promoted",
            "genus": row.get("genus") or "",
            "species": row.get("species") or "",
            "strain": None,
            "contig_id": row.get("hostchr") or "",
            "element_start": row.get("shipstart"),
            "element_end": row.get("shipend"),
            "element_strand": strand,
            "curator": row.get("uploader") or "",
            "curated_status": "curated",
            "notes": row.get("comment") or "",
            "duplicate_info": duplicate_info,
        }

        if row.get("classification_family") or row.get("classification_navis"):
            processed_data["classification"] = {
                "family": row.get("classification_family"),
                "navis": row.get("classification_navis"),
                "haplotype": row.get("classification_haplotype"),
                "source": row.get("classification_source"),
                "closest_match": row.get("closest_match"),
                "confidence": row.get("classification_confidence"),
            }
            if row.get("classification_family"):
                processed_data["ship_family"] = row["classification_family"]
            if row.get("classification_navis"):
                processed_data["ship_navis"] = row["classification_navis"]
            if row.get("classification_haplotype"):
                processed_data["ship_haplotype"] = row["classification_haplotype"]

        result = perform_database_insertion(
            processed_data,
            anno_contents=row.get("anno_contents"),
            anno_filename=row.get("anno_filename"),
            anno_date=None,
            seq_date=0,
        )

        with get_submissions_session() as session:
            session.execute(
                text("UPDATE submissions SET needs_review = 0 WHERE id = :id"),
                {"id": sub_id},
            )
            session.commit()

        logger.info(
            "Promoted submission %s → accession %s", sub_id, result["accession"]
        )
        return True, result["accession"], None

    except Exception as exc:
        logger.error("Promotion failed for submission %s: %s", sub_id, exc)
        return False, None, str(exc)


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

    rows = [{**r, "_dirty": False} for r in df.fillna("").to_dict("records")]

    return dag.AgGrid(
        id=grid_id,
        columnDefs=col_defs,
        rowData=rows,
        defaultColDef={"resizable": True, "minWidth": 80},
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": 50,
            "suppressPropertyNamesCheck": True,
            "rowHeight": 40,
            "headerHeight": 44,
            "rowClassRules": {"admin-unsaved-row": "params.data._dirty === true"},
        },
        getRowId="params.data.id",
        className="ag-theme-alpine",
        style={"width": "100%", "height": "600px"},
    )


def _promote_modal():
    return dmc.Modal(
        id="admin-promote-modal",
        title="Promote submission to main database?",
        centered=True,
        size="lg",
        opened=False,
        children=dmc.Stack(
            [
                dmc.Alert(
                    "This inserts the submission into the live starbase database. "
                    "The action cannot be undone — a duplicate check will still run.",
                    color="yellow",
                    variant="light",
                ),
                html.Div(id="admin-promote-modal-info"),
                dmc.Group(
                    [
                        dmc.Button(
                            "Promote",
                            id="admin-promote-confirm",
                            color="green",
                        ),
                        dmc.Button(
                            "Cancel",
                            id="admin-promote-cancel",
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
        fam_df = _fetch_family_names()
        nav_df = _fetch_navis_names()
        hap_df = _fetch_haplotype_names()
        acc_df = _fetch_accessions()
        sacc_df = _fetch_ship_accessions()
        gen_df = _fetch_genomes()
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
            dmc.Paper(
                dmc.Group(
                    [
                        dmc.Text(
                            id="admin-pending-label",
                            size="sm",
                            c="dimmed",
                            children="No pending changes. Yellow cells are editable — edit freely, then save.",
                        ),
                        dmc.Group(
                            [
                                dmc.Button(
                                    "Discard",
                                    id="admin-discard-btn",
                                    variant="subtle",
                                    color="gray",
                                    size="sm",
                                    disabled=True,
                                ),
                                dmc.Button(
                                    "Save changes",
                                    id="admin-save-btn",
                                    color="blue",
                                    size="sm",
                                    disabled=True,
                                ),
                            ],
                            gap="xs",
                        ),
                    ],
                    justify="space-between",
                    align="center",
                ),
                p="sm",
                mb="md",
                radius="sm",
                withBorder=True,
                style={"borderColor": "var(--mantine-color-blue-3)"},
            ),
            dbc.Tabs(
                [
                    dbc.Tab(
                        _make_grid(
                            js_df,
                            "admin-joined-ships-grid",
                            EDITABLE_COLS["joined_ships"],
                        ),
                        label=f"Starships ({len(js_df):,})",
                        tab_id="joined_ships",
                    ),
                    dbc.Tab(
                        html.Div(
                            [
                                dmc.Group(
                                    [
                                        dmc.Text(
                                            "Select row(s) to promote into the main starbase database.",
                                            size="sm",
                                            c="dimmed",
                                        ),
                                        dmc.Button(
                                            "Promote to Main DB",
                                            id="admin-promote-btn",
                                            color="green",
                                            size="sm",
                                            disabled=True,
                                        ),
                                    ],
                                    justify="space-between",
                                    align="center",
                                    mb="xs",
                                    mt="xs",
                                ),
                                _make_grid(
                                    sub_df,
                                    "admin-submissions-grid",
                                    EDITABLE_COLS["submissions"],
                                ),
                            ]
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
                    dbc.Tab(
                        _make_grid(
                            fam_df,
                            "admin-family-names-grid",
                            EDITABLE_COLS["family_names"],
                        ),
                        label=f"Families ({len(fam_df):,})",
                        tab_id="family_names",
                    ),
                    dbc.Tab(
                        _make_grid(
                            nav_df,
                            "admin-navis-names-grid",
                            EDITABLE_COLS["navis_names"],
                        ),
                        label=f"Navis ({len(nav_df):,})",
                        tab_id="navis_names",
                    ),
                    dbc.Tab(
                        _make_grid(
                            hap_df,
                            "admin-haplotype-names-grid",
                            EDITABLE_COLS["haplotype_names"],
                        ),
                        label=f"Haplotypes ({len(hap_df):,})",
                        tab_id="haplotype_names",
                    ),
                    dbc.Tab(
                        _make_grid(
                            acc_df,
                            "admin-accessions-grid",
                            EDITABLE_COLS["accessions"],
                        ),
                        label=f"Accessions ({len(acc_df):,})",
                        tab_id="accessions",
                    ),
                    dbc.Tab(
                        _make_grid(
                            sacc_df,
                            "admin-ship-accessions-grid",
                            EDITABLE_COLS["ship_accessions"],
                        ),
                        label=f"Ship Accessions ({len(sacc_df):,})",
                        tab_id="ship_accessions",
                    ),
                    dbc.Tab(
                        _make_grid(
                            gen_df,
                            "admin-genomes-grid",
                            EDITABLE_COLS["genomes"],
                        ),
                        label=f"Genomes ({len(gen_df):,})",
                        tab_id="genomes",
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
        dcc.Store(id="admin-pending-changes", data=[]),
        dcc.Store(id="admin-selected-submissions", data=[]),
        html.Div(id="admin-content"),
        _version_bump_modal(),
        _promote_modal(),
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
# Clientside callbacks — cell edits accumulate into pending store + _dirty
# ---------------------------------------------------------------------------

_CELL_JS = """
function(ev, pending, rowData) {{
    if (!ev || ev.colId === undefined) return [window.dash_clientside.no_update, window.dash_clientside.no_update];
    var rid = ev.data.id, cid = ev.colId, ov = ev.oldValue, tk = "{tk}";
    var nv = (ev.data && ev.data[cid] !== undefined) ? ev.data[cid] : ev.value;
    var p = JSON.parse(JSON.stringify(pending || []));
    var idx = -1;
    for (var i = 0; i < p.length; i++) {{ if (p[i].table===tk && p[i].row_id===rid && p[i].col_id===cid) {{ idx=i; break; }} }}
    if (idx >= 0) {{ p[idx].new_value = nv; }} else {{ p.push({{table:tk,row_id:rid,col_id:cid,old_value:ov,new_value:nv}}); }}
    var rd = (rowData||[]).map(function(r) {{ return r.id===rid ? Object.assign({{}},r,{{_dirty:true}}) : r; }});
    return [rd, p];
}}
"""

for _tk, _gid in GRID_IDS.items():
    clientside_callback(
        _CELL_JS.format(tk=_tk),
        Output(_gid, "rowData"),
        Output("admin-pending-changes", "data", allow_duplicate=True),
        Input(_gid, "cellValueChanged"),
        State("admin-pending-changes", "data"),
        State(_gid, "rowData"),
        prevent_initial_call=True,
    )

# ---------------------------------------------------------------------------
# Clientside callback — relay selectedRows into static-layout store
# (avoids refErr for State("admin-submissions-grid", "selectedRows") in Python
# callbacks that are validated against the initial layout before auth)
# ---------------------------------------------------------------------------

clientside_callback(
    "function(r) { return r || []; }",
    Output("admin-selected-submissions", "data"),
    Input("admin-submissions-grid", "selectedRows"),
    prevent_initial_call=True,
)

# ---------------------------------------------------------------------------
# Clientside callback — pending count drives toolbar labels / disabled state
# ---------------------------------------------------------------------------

clientside_callback(
    """
    function(p) {
        var n = (p||[]).length, dis = n===0;
        var lbl = n > 0 ? ('Save '+n+' change'+(n!==1?'s':'')) : 'Save changes';
        var hint = n > 0 ? (n+' unsaved change'+(n!==1?'s':'')+' \u2014 click Save or Discard.')
                         : 'No pending changes. Yellow cells are editable \u2014 edit freely, then save.';
        return [dis, dis, lbl, hint];
    }
    """,
    Output("admin-save-btn", "disabled"),
    Output("admin-discard-btn", "disabled"),
    Output("admin-save-btn", "children"),
    Output("admin-pending-label", "children"),
    Input("admin-pending-changes", "data"),
    prevent_initial_call=True,
)


# ---------------------------------------------------------------------------
# Callbacks — batch save / discard
# ---------------------------------------------------------------------------


def _make_version_modal_content(successes):
    tables = sorted({c["table"] for c in successes})
    n = len(successes)
    cur = get_database_version()
    info = dmc.Stack(
        [
            dmc.Text(
                f"Saved {n} change{'s' if n != 1 else ''} across: {', '.join(tables)}",
                size="sm",
                fw=500,
            ),
            dmc.Divider(),
            dmc.Group(
                [
                    dmc.Text("Current:", size="sm"),
                    dmc.Badge(cur, variant="outline", color="gray"),
                ],
                gap="xs",
            ),
            dmc.Group(
                [
                    dmc.Text("Minor \u2192", size="sm"),
                    dmc.Badge(
                        _bump_version(cur, "minor"), variant="light", color="blue"
                    ),
                    dmc.Text("Patch \u2192", size="sm"),
                    dmc.Badge(
                        _bump_version(cur, "patch"), variant="light", color="teal"
                    ),
                ],
                gap="xs",
            ),
        ],
        gap="xs",
    )
    desc = f"admin: {n} edit{'s' if n != 1 else ''} across {', '.join(tables)}"
    return info, desc


@callback(
    [
        Output("admin-pending-changes", "data"),
        *[Output(gid, "rowData", allow_duplicate=True) for gid in GRID_IDS.values()],
        Output("admin-version-modal", "opened"),
        Output("admin-version-modal-info", "children"),
        Output("admin-version-description", "value"),
        Output("notifications-container", "children"),
    ],
    Input("admin-save-btn", "n_clicks"),
    State("admin-pending-changes", "data"),
    prevent_initial_call=True,
)
def save_all_pending(n_clicks, pending):
    n_out = (
        len(GRID_IDS) + 5
    )  # pending-store + 10 grids + modal-opened + modal-info + description + notif
    if not n_clicks or not pending:
        return [no_update] * n_out

    logger.info("save_all_pending: %d change(s): %s", len(pending), pending)
    errors, successes = [], []
    for ch in pending:
        ok, err = _do_update(
            ch["table"], ch["row_id"], ch["col_id"], ch.get("new_value")
        )
        (successes if ok else errors).append((ch, err))

    changed_tables = {ch["table"] for ch in pending}
    fresh = [
        _refetch_rowdata(k) if k in changed_tables else no_update for k in GRID_IDS
    ]

    if errors:
        msgs = "; ".join(
            f"{c['table']}.{c['col_id']} row {c['row_id']}: {e}" for c, e in errors[:3]
        )
        notif = dmc.Notification(
            id=f"admin-err-{uuid.uuid4().hex[:6]}",
            title="Some saves failed",
            message=msgs,
            color="red",
            action="show",
            autoClose=10000,
        )
        return [[], *fresh, False, no_update, no_update, notif]

    info, desc = _make_version_modal_content([c for c, _ in successes])
    return [[], *fresh, True, info, desc, no_update]


@callback(
    [
        Output("admin-pending-changes", "data", allow_duplicate=True),
        *[Output(gid, "rowData", allow_duplicate=True) for gid in GRID_IDS.values()],
    ],
    Input("admin-discard-btn", "n_clicks"),
    prevent_initial_call=True,
)
def discard_pending(n_clicks):
    if not n_clicks:
        return [no_update] * (len(GRID_IDS) + 1)
    return [[], *[_refetch_rowdata(k) for k in GRID_IDS]]


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


# ---------------------------------------------------------------------------
# Callbacks — promote submission
# ---------------------------------------------------------------------------


@callback(
    Output("admin-promote-btn", "disabled"),
    Input("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)
def toggle_promote_button(selected):
    return not bool(selected)


@callback(
    [
        Output("admin-promote-modal", "opened"),
        Output("admin-promote-modal-info", "children"),
    ],
    Input("admin-promote-btn", "n_clicks"),
    State("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)
def open_promote_modal(n_clicks, selected_rows):
    if not n_clicks or not selected_rows:
        return False, no_update

    cards = []
    for row in selected_rows:
        family = row.get("classification_family") or "unclassified"
        navis = row.get("classification_navis") or ""
        taxon = f"{row.get('genus', '')} {row.get('species', '')}".strip()
        cards.append(
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Group(
                            [
                                dmc.Text(f"ID {row.get('id')}", fw=700, size="sm"),
                                dmc.Text(
                                    row.get("seq_filename", ""), size="sm", c="dimmed"
                                ),
                                dmc.Badge(
                                    "needs review"
                                    if row.get("needs_review")
                                    else "reviewed",
                                    color="orange"
                                    if row.get("needs_review")
                                    else "green",
                                    variant="light",
                                    size="xs",
                                ),
                            ],
                            gap="xs",
                        ),
                        dmc.Group(
                            [
                                dmc.Badge(
                                    taxon or "unknown taxon",
                                    variant="light",
                                    color="gray",
                                ),
                                dmc.Badge(family, variant="light", color="blue"),
                                dmc.Badge(navis, variant="light", color="teal")
                                if navis
                                else None,
                            ],
                            gap="xs",
                        ),
                    ],
                    gap="xs",
                ),
                p="sm",
                withBorder=True,
                radius="sm",
            )
        )

    return True, dmc.Stack(cards, gap="xs")


@callback(
    [
        Output("admin-promote-modal", "opened", allow_duplicate=True),
        Output("admin-version-modal", "opened", allow_duplicate=True),
        Output("admin-version-modal-info", "children", allow_duplicate=True),
        Output("admin-version-description", "value", allow_duplicate=True),
        Output("notifications-container", "children", allow_duplicate=True),
    ],
    Input("admin-promote-confirm", "n_clicks"),
    State("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)
def run_promotion(n_clicks, selected_rows):
    if not n_clicks or not selected_rows:
        return no_update, no_update, no_update, no_update, no_update

    results = []
    for row in selected_rows:
        sub_id = row.get("id")
        success, accession, error = _promote_submission(sub_id)
        results.append((sub_id, success, accession, error))

    failures = [r for r in results if not r[1]]

    if failures:
        errors = "; ".join(f"sub {r[0]}: {r[3]}" for r in failures)
        notif = dmc.Notification(
            id=f"admin-err-{uuid.uuid4().hex[:6]}",
            title="Promotion failed",
            message=errors,
            color="red",
            action="show",
            autoClose=10000,
        )
        return False, False, no_update, no_update, notif

    successes = [r for r in results if r[1]]
    accessions = ", ".join(r[2] for r in successes if r[2])
    _pseudo_changes = [
        {"table": "submissions", "col_id": f"promote→{a}", "row_id": sid}
        for sid, _, a, _ in successes
    ]
    n = len(successes)
    cur = get_database_version()
    info = dmc.Stack(
        [
            dmc.Text(
                f"Promoted {n} submission{'s' if n != 1 else ''}: {accessions}",
                size="sm",
                fw=500,
            ),
            dmc.Divider(),
            dmc.Group(
                [
                    dmc.Text("Current:", size="sm"),
                    dmc.Badge(cur, variant="outline", color="gray"),
                ],
                gap="xs",
            ),
            dmc.Group(
                [
                    dmc.Text("Minor \u2192", size="sm"),
                    dmc.Badge(
                        _bump_version(cur, "minor"), variant="light", color="blue"
                    ),
                    dmc.Text("Patch \u2192", size="sm"),
                    dmc.Badge(
                        _bump_version(cur, "patch"), variant="light", color="teal"
                    ),
                ],
                gap="xs",
            ),
        ],
        gap="xs",
    )
    desc = f"admin: promoted submission(s) {accessions} to main DB"
    return False, True, info, desc, no_update


@callback(
    Output("admin-promote-modal", "opened", allow_duplicate=True),
    Input("admin-promote-cancel", "n_clicks"),
    prevent_initial_call=True,
)
def cancel_promote(n_clicks):
    return False

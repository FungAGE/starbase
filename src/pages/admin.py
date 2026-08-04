import time
import uuid
from urllib.parse import parse_qs

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
from dash.exceptions import PreventUpdate
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
    "joined_ships": {"curated_status", "evidence", "source", "tax_id", "accession_id"},
    "submissions": {
        "needs_review",
        "classification_family",
        "classification_navis",
        "classification_haplotype",
        "classification_confidence",
        "comment",
    },
    "taxonomy": {
        "name",
        "taxID",
        "genus",
        "species",
        "strain",
        "superkingdom",
        "clade",
        "kingdom",
        "subkingdom",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        "order",
        "suborder",
        "family",
        "section",
        "species_group",
        "subgenus",
    },
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
        "taxonomy_id",
    },
    "ships": {"md5", "rev_comp_md5"},
    "starship_features": {
        "contigID",
        "elementBegin",
        "elementEnd",
        "elementLength",
        "strand",
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
    "ships": {"sql_table": "ships", "db": "starbase", "pk": "id"},
    "starship_features": {
        "sql_table": "starship_features",
        "db": "starbase",
        "pk": "id",
    },
}

# SQLite reserved / quoted column names in UPDATE statements
SQL_QUOTED_COLS = {"class", "order", "source"}


def _sql_col_ref(col_id):
    if col_id in SQL_QUOTED_COLS:
        return f"`{col_id}`"
    return col_id


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
            SELECT j.id, j.starshipID, j.curated_status, j.evidence, j.source, j.tax_id,
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
            SELECT id, submission_group_id, seq_filename, uploader, seq_date,
                   genus, species, strain, hostchr, shipstart, shipend, shipstrand,
                   assembly_accession, evidence,
                   processing_status, accession_tag,
                   needs_review, comment,
                   classification_family, classification_navis,
                   classification_haplotype, classification_confidence
            FROM submissions
            ORDER BY submission_group_id DESC, id ASC
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
                   g.biosample, g.acquisition_date, g.assembly_accession, g.taxonomy_id,
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


def _do_insert(table_key, col_values):
    """Run a whitelisted INSERT for a new row. Returns (success, error, new_id)."""
    allowed = EDITABLE_COLS.get(table_key, set())
    config = TABLE_CONFIG[table_key]

    cols = [c for c in col_values if c in allowed and col_values[c] not in (None, "")]
    if not cols:
        return False, "No editable field was filled in for the new row.", None

    values = {}
    for c in cols:
        v = col_values[c]
        if c in BOOLEAN_COLS:
            v = 1 if str(v).lower() in ("true", "1", "yes") else 0
        values[c] = v

    col_refs = ", ".join(_sql_col_ref(c) for c in cols)
    placeholders = ", ".join(f":{c}" for c in cols)
    sql = text(
        f"INSERT INTO {config['sql_table']} ({col_refs}) VALUES ({placeholders})"
    )

    try:
        ctx = (
            get_starbase_session
            if config["db"] == "starbase"
            else get_submissions_session
        )
        with ctx() as session:
            result = session.execute(sql, values)
            session.commit()
            new_id = result.lastrowid
        logger.info(
            "admin INSERT %s: %r (id=%s)", config["sql_table"], values, new_id
        )
        return True, None, new_id
    except Exception as exc:
        logger.error("admin INSERT error (%s): %s", table_key, exc)
        return False, str(exc), None


def _do_update(table_key, row_id, col_id, new_value):
    """Run a whitelisted UPDATE. Returns (success: bool, error: str|None)."""
    allowed = EDITABLE_COLS.get(table_key, set())
    if col_id not in allowed:
        return False, f"Column '{col_id}' is not editable."

    config = TABLE_CONFIG[table_key]

    if col_id in BOOLEAN_COLS:
        new_value = 1 if str(new_value).lower() in ("true", "1", "yes") else 0

    sql = text(
        f"UPDATE {config['sql_table']} SET {_sql_col_ref(col_id)} = :val WHERE {config['pk']} = :pk"
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
# Admin jobs — consistency / cleanup utilities
# ---------------------------------------------------------------------------

_REPORT_COLS = [
    "ship_id",
    "accession_tag",
    "md5",
    "rev_comp_md5",
    "issue",
]

_STAGE_COLS = [
    "table",
    "row_id",
    "col_id",
    "old_value",
    "new_value",
]

_TAX_LINK_COLS = [
    "table",
    "row_id",
    "col_id",
    "old_value",
    "new_value",
    "strategy",
    "starshipID",
    "note",
]

_TAX_NCBI_COLS = [
    "genome_id",
    "assembly_accession",
    "biosample",
    "ncbi_taxid",
    "organism",
    "taxonomy_id",
    "joined_ships_updated",
    "action",
]

_TAX_VALIDATE_COLS = [
    "issue_type",
    "taxonomy_id",
    "name",
    "detail",
]

_ACCESSION_MISSING_COLS = [
    "ship_id",
    "sequence_preview",
    "sequence_length",
    "md5",
    "rev_comp_md5",
]

_ACCESSION_ASSIGN_COLS = [
    "ship_id",
    "new_accession",
    "needs_review",
    "sequence_length",
    "note",
]

_ACCESSION_SYNC_COLS = [
    "table",
    "row_id",
    "col_id",
    "old_value",
    "new_value",
    "ship_id",
    "starshipID",
]

_ACCESSION_DIAG_COLS = [
    "issue_type",
    "table",
    "row_id",
    "ship_id",
    "detail",
]

_ACCESSION_CLEANUP_COLS = [
    "type",
    "primary_accession",
    "secondary_accessions",
    "reason",
]

_GENOME_COORD_COLS = [
    "issue_type",
    "source",
    "joined_ship_id",
    "submission_id",
    "starshipID",
    "ship_id",
    "assembly_accession",
    "genome_source",
    "contig_id",
    "coordinates",
    "detail",
]

_GENOME_COORD_FIX_COLS = [
    "status",
    "starshipID",
    "ship_id",
    "starship_feature_id",
    "assembly_accession",
    "contig_id",
    "old_coordinates",
    "new_coordinates",
    "strand",
    "coverage",
    "identity",
    "detail",
]

# Cap assign-missing preview in admin (full pipeline per ship is expensive)
_ACCESSION_ASSIGN_PREVIEW_LIMIT = 25


def _taxonomy_orm_field_to_col(field):
    if field == "class_":
        return "class"
    return field


def _consolidator_clean_to_changes(report):
    """Convert TaxonomyConsolidator.clean_taxonomy_data report to pending changes."""
    changes = []
    for item in report.get("whitespace_fixes", []):
        col = _taxonomy_orm_field_to_col(item.get("field", ""))
        nv = item.get("new_value")
        if nv == "NULL":
            nv = None
        ov = item.get("old_value")
        if ov == '""':
            ov = ""
        changes.append(
            {
                "table": "taxonomy",
                "row_id": item["taxonomy_id"],
                "col_id": col,
                "old_value": ov if ov is not None else "",
                "new_value": nv,
                "source": "job",
            }
        )
    for item in report.get("species_cleanup", []):
        changes.append(
            {
                "table": "taxonomy",
                "row_id": item["taxonomy_id"],
                "col_id": "species",
                "old_value": item.get("old_species", ""),
                "new_value": item.get("new_species"),
                "source": "job",
            }
        )
    for item in report.get("strain_cleanup", []):
        changes.append(
            {
                "table": "taxonomy",
                "row_id": item["taxonomy_id"],
                "col_id": "strain",
                "old_value": item.get("old_strain", ""),
                "new_value": item.get("new_strain"),
                "source": "job",
            }
        )
    for item in report.get("name_cleanup", []):
        changes.append(
            {
                "table": "taxonomy",
                "row_id": item["taxonomy_id"],
                "col_id": "name",
                "old_value": item.get("old_name", ""),
                "new_value": item.get("new_name"),
                "source": "job",
            }
        )
    return changes


def _consolidator_backfill_to_changes(report):
    changes = []
    for item in report.get("fields_backfilled", []):
        col = _taxonomy_orm_field_to_col(item.get("field", ""))
        changes.append(
            {
                "table": "taxonomy",
                "row_id": item["taxonomy_id"],
                "col_id": col,
                "old_value": "",
                "new_value": item.get("inferred_value"),
                "source": "job",
            }
        )
    return changes


def _flatten_validation_report(report):
    rows = []
    for item in report.get("duplicate_entries", []):
        ids = ", ".join(str(e.get("id")) for e in item.get("entries", []))
        rows.append(
            {
                "issue_type": "duplicate",
                "taxonomy_id": ids,
                "name": item.get("name") or "",
                "detail": (
                    f"{item.get('count')} rows — "
                    f"{item.get('genus')} {item.get('species')} "
                    f"strain={item.get('strain') or ''}"
                ).strip(),
            }
        )
    for item in report.get("orphaned_entries", []):
        rows.append(
            {
                "issue_type": "orphaned",
                "taxonomy_id": item.get("taxonomy_id"),
                "name": item.get("name") or "",
                "detail": f"taxID={item.get('taxID') or 'none'}",
            }
        )
    for item in report.get("family_inconsistencies", []):
        rows.append(
            {
                "issue_type": "family_inconsistency",
                "taxonomy_id": "",
                "name": item.get("family") or "",
                "detail": f"{item.get('field')}: {item.get('values')}",
            }
        )
    for item in report.get("data_quality_issues", []):
        rows.append(
            {
                "issue_type": item.get("type", "data_quality"),
                "taxonomy_id": item.get("count", ""),
                "name": "",
                "detail": item.get("description", str(item)),
            }
        )
    return rows


def _job_taxonomy_validate():
    from src.database.cleanup.utils.taxonomy_consolidator import TaxonomyConsolidator

    with TaxonomyConsolidator() as consolidator:
        report = consolidator.validate_taxonomy_consistency()
    rows = _flatten_validation_report(report)
    s = report.get("summary", {})
    summary = (
        f"{len(rows)} issue(s): "
        f"{s.get('duplicate_entries', 0)} duplicate group(s), "
        f"{s.get('orphaned_entries', 0)} orphaned, "
        f"{s.get('family_inconsistencies', 0)} family inconsistency group(s)"
        if rows
        else "No taxonomy consistency issues found."
    )
    return {
        "job": "taxonomy_validate",
        "mode": "report",
        "columns": _TAX_VALIDATE_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_taxonomy_clean():
    from src.database.cleanup.utils.taxonomy_consolidator import TaxonomyConsolidator

    with TaxonomyConsolidator() as consolidator:
        report = consolidator.clean_taxonomy_data(dry_run=True)
    proposed = _consolidator_clean_to_changes(report)
    preview = [
        {
            "table": c["table"],
            "row_id": c["row_id"],
            "col_id": c["col_id"],
            "old_value": c["old_value"],
            "new_value": c["new_value"],
        }
        for c in proposed
    ]
    s = report.get("summary", {})
    n = len(proposed)
    summary = (
        f"{n} field fix(es) to stage: "
        f"{s.get('whitespace_fixes', 0)} whitespace, "
        f"{s.get('species_cleanup', 0)} species, "
        f"{s.get('strain_cleanup', 0)} strain, "
        f"{s.get('name_cleanup', 0)} name"
        if n
        else "Taxonomy data already clean — nothing to stage."
    )
    return {
        "job": "taxonomy_clean",
        "mode": "stage",
        "columns": _STAGE_COLS,
        "rows": preview,
        "proposed_changes": proposed,
        "summary": summary,
    }


def _job_taxonomy_hierarchy_backfill():
    from src.database.cleanup.utils.taxonomy_consolidator import TaxonomyConsolidator

    with TaxonomyConsolidator() as consolidator:
        report = consolidator.backfill_taxonomy_hierarchy(dry_run=True)
    proposed = _consolidator_backfill_to_changes(report)
    preview = [
        {
            "table": c["table"],
            "row_id": c["row_id"],
            "col_id": c["col_id"],
            "old_value": c["old_value"],
            "new_value": c["new_value"],
        }
        for c in proposed
    ]
    n = len(proposed)
    summary = (
        f"{n} empty hierarchy field(s) can be backfilled from matching family/genus rows"
        if n
        else "No empty hierarchy fields to backfill."
    )
    return {
        "job": "taxonomy_hierarchy_backfill",
        "mode": "stage",
        "columns": _STAGE_COLS,
        "rows": preview,
        "proposed_changes": proposed,
        "summary": summary,
    }


def _job_taxonomy_link_ships():
    from src.database.cleanup.utils.fill_taxonomy_from_ome_and_source import (
        analyze_fill_taxonomy_from_ome_and_source,
    )

    proposed, preview, counts = analyze_fill_taxonomy_from_ome_and_source(
        overwrite=False
    )
    n = len(proposed)
    create_n = counts.get("taxonomy_would_create", 0)
    summary = (
        f"{n} link(s) to stage: genome={counts.get('from_genome', 0)}, "
        f"ome={counts.get('from_ome', 0)}, "
        f"gluck_thaler_2025={counts.get('from_gluck_thaler_2025', 0)}, "
        f"genome taxonomy_id={counts.get('genome_taxonomy_updates', 0)}"
        + (
            f"; {create_n} new taxonomy row(s) need CLI fill (not stageable)"
            if create_n
            else ""
        )
        if n or create_n
        else "All joined_ships already linked to taxonomy."
    )
    return {
        "job": "taxonomy_link_ships",
        "mode": "stage",
        "columns": _TAX_LINK_COLS,
        "rows": preview,
        "proposed_changes": proposed,
        "summary": summary,
    }


def _job_md5_validation():
    """Report ships with duplicate MD5 hashes or missing rev_comp_md5."""
    with get_starbase_session() as session:
        df = pd.read_sql_query(
            """
            WITH dup_md5 AS (
                SELECT md5 FROM ships
                WHERE md5 IS NOT NULL AND md5 != ''
                GROUP BY md5 HAVING COUNT(*) > 1
            ),
            dup_rev AS (
                SELECT rev_comp_md5 FROM ships
                WHERE rev_comp_md5 IS NOT NULL AND rev_comp_md5 != ''
                GROUP BY rev_comp_md5 HAVING COUNT(*) > 1
            )
            SELECT s.id AS ship_id,
                   a.accession_tag,
                   s.md5,
                   s.rev_comp_md5,
                   CASE
                       WHEN s.rev_comp_md5 IS NULL OR s.rev_comp_md5 = ''
                           THEN 'missing rev_comp_md5'
                       WHEN s.md5 IN (SELECT md5 FROM dup_md5)
                           THEN 'duplicate md5'
                       WHEN s.rev_comp_md5 IN (SELECT rev_comp_md5 FROM dup_rev)
                           THEN 'duplicate rev_comp_md5'
                   END AS issue
            FROM ships s
            LEFT JOIN accessions a ON s.accession_id = a.id
            WHERE (s.rev_comp_md5 IS NULL OR s.rev_comp_md5 = '')
               OR s.md5 IN (SELECT md5 FROM dup_md5)
               OR s.rev_comp_md5 IN (SELECT rev_comp_md5 FROM dup_rev)
            ORDER BY issue, s.md5, s.id
            """,
            session.bind,
        )

    rows = df.fillna("").to_dict("records")
    n_dup = sum(1 for r in rows if "duplicate" in r.get("issue", ""))
    n_missing = sum(1 for r in rows if "missing" in r.get("issue", ""))
    summary = (
        f"{len(rows)} issue(s): {n_dup} duplicate hash row(s), "
        f"{n_missing} missing rev_comp_md5"
        if rows
        else "No MD5 issues found."
    )
    return {
        "job": "md5_validation",
        "mode": "report",
        "columns": _REPORT_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_backfill_rev_comp_md5():
    """Stage rev_comp_md5 backfill for ships where it is NULL."""
    from src.utils.classification_utils import generate_md5_hash
    from src.utils.seq_utils import revcomp

    with get_starbase_session() as session:
        df = pd.read_sql_query(
            """
            SELECT s.id, s.sequence, s.rev_comp_md5, a.accession_tag
            FROM ships s
            LEFT JOIN accessions a ON s.accession_id = a.id
            WHERE s.rev_comp_md5 IS NULL OR s.rev_comp_md5 = ''
            ORDER BY s.id
            """,
            session.bind,
        )

    proposed = []
    preview_rows = []
    for _, row in df.iterrows():
        seq = row.get("sequence") or ""
        if not seq:
            preview_rows.append(
                {
                    "table": "ships",
                    "row_id": row["id"],
                    "col_id": "rev_comp_md5",
                    "old_value": "",
                    "new_value": "",
                    "note": "skipped — no sequence",
                }
            )
            continue
        new_hash = generate_md5_hash(revcomp(seq))
        proposed.append(
            {
                "table": "ships",
                "row_id": int(row["id"]),
                "col_id": "rev_comp_md5",
                "old_value": row.get("rev_comp_md5") or "",
                "new_value": new_hash,
                "source": "job",
            }
        )
        preview_rows.append(
            {
                "table": "ships",
                "row_id": int(row["id"]),
                "col_id": "rev_comp_md5",
                "old_value": row.get("rev_comp_md5") or "",
                "new_value": new_hash,
            }
        )

    summary = (
        f"{len(proposed)} change(s) ready to stage"
        if proposed
        else "No ships need rev_comp_md5 backfill."
    )
    return {
        "job": "backfill_rev_comp_md5",
        "mode": "stage",
        "columns": _STAGE_COLS,
        "rows": preview_rows,
        "proposed_changes": proposed,
        "summary": summary,
    }


def _job_accession_status():
    """Report ships that have sequences but no accession_id."""
    from src.database.cleanup.utils.accession_manager import AccessionManager

    manager = AccessionManager()
    missing_df = manager.check_missing_accessions()
    rows = missing_df.fillna("").to_dict("records") if not missing_df.empty else []
    stats = manager.stats
    n = len(rows)
    summary = (
        f"{n} ship(s) missing accessions out of "
        f"{stats.get('ships_checked', 0)} with sequences"
        if n
        else "All sequenced ships have accessions."
    )
    return {
        "job": "accession_status",
        "mode": "report",
        "columns": _ACCESSION_MISSING_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_accession_assign_preview():
    """
    Preview accession assignments for ships missing accessions (dry-run).

    Creates new accession rows + ships.accession_id updates — not stageable via
    admin save; use CLI accession_manager.py --fix-missing --apply to apply.
    """
    from src.database.cleanup.utils.accession_manager import AccessionManager

    manager = AccessionManager()
    results = manager.assign_missing_accessions(
        dry_run=True, limit=_ACCESSION_ASSIGN_PREVIEW_LIMIT
    )
    rows = []
    for a in results.get("assignments") or []:
        rows.append(
            {
                "ship_id": a.get("ship_id"),
                "new_accession": a.get("new_accession", ""),
                "needs_review": a.get("needs_review", False),
                "sequence_length": a.get("sequence_length", ""),
                "note": "CLI apply only — creates accession rows",
            }
        )
    for err in results.get("errors") or []:
        rows.append(
            {
                "ship_id": err.get("ship_id", ""),
                "new_accession": "",
                "needs_review": "",
                "sequence_length": "",
                "note": f"error: {err.get('error', err)}",
            }
        )
    n = len(results.get("assignments") or [])
    err_n = len(results.get("errors") or [])
    summary = (
        f"Preview ({_ACCESSION_ASSIGN_PREVIEW_LIMIT} ship max): "
        f"{n} assignment(s), {err_n} error(s). "
        "Not stageable — run accession_manager.py --fix-missing --apply."
        if n or err_n
        else "No ships missing accessions."
    )
    return {
        "job": "accession_assign_preview",
        "mode": "report",
        "columns": _ACCESSION_ASSIGN_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_accession_sync():
    """Stage joined_ships.accession_id from ships.accession_id where out of sync."""
    from src.database.cleanup.utils.sync_accession_ids import analyze_sync_accession_ids

    proposed, preview, stats = analyze_sync_accession_ids()
    n = len(proposed)
    summary = (
        f"{n} joined_ships.accession_id update(s) to stage "
        f"({stats.get('already_synced', 0)} already synced, "
        f"{stats.get('checked', 0)} checked)"
        if n
        else f"All {stats.get('checked', 0)} joined_ships accession_ids are in sync."
    )
    return {
        "job": "accession_sync",
        "mode": "stage",
        "columns": _ACCESSION_SYNC_COLS,
        "rows": preview,
        "proposed_changes": proposed,
        "summary": summary,
    }


def _job_accession_diagnose():
    """Report accession sync, FK, duplicate, and orphan issues."""
    from src.database.cleanup.utils.diagnose_accession_sync import (
        analyze_accession_sync_issues,
    )

    rows, counts = analyze_accession_sync_issues()
    n = len(rows)
    summary = (
        f"{n} issue(s): "
        f"{counts.get('out_of_sync', 0)} out-of-sync, "
        f"{counts.get('fk_mismatch', 0)} FK mismatch, "
        f"{counts.get('duplicate_ship_id_groups', 0)} duplicate ship_id group(s), "
        f"{counts.get('orphaned_joined_ships', 0)} orphaned"
        if n
        else "No accession sync issues found."
    )
    return {
        "job": "accession_diagnose",
        "mode": "report",
        "columns": _ACCESSION_DIAG_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_accession_fix_duplicates_preview():
    """Preview duplicate joined_ships.ship_id cleanup (report only)."""
    from src.database.cleanup.utils.diagnose_accession_sync import (
        analyze_fix_duplicate_ship_ids,
    )

    rows, stats = analyze_fix_duplicate_ship_ids()
    n_clear = stats.get("would_clear", 0)
    summary = (
        f"{stats.get('checked', 0)} duplicate row(s): "
        f"{stats.get('kept', 0)} keep, {n_clear} would clear ship_id/accession_id. "
        "Not stageable — run diagnose_accession_sync.py --fix-duplicates --apply."
        if stats.get("checked")
        else "No duplicate ship_id groups in joined_ships."
    )
    return {
        "job": "accession_fix_duplicates_preview",
        "mode": "report",
        "columns": _ACCESSION_DIAG_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_accession_cleanup_analyze():
    """
    Report hash-duplicate and reverse-complement accession consolidation groups.

    Full consolidation merges ships/FKs — use cleanup_accessions.py CLI to apply.
    Nested-sequence detection is skipped here (too slow for interactive admin).
    """
    from src.database.cleanup.utils.cleanup_accessions import analyze_accession_cleanup

    _consolidations, preview, stats = analyze_accession_cleanup(include_nested=False)
    n = stats.get("total_groups", 0)
    sec = stats.get("secondary_accessions", 0)
    summary = (
        f"{n} consolidation group(s) affecting {sec} secondary accession(s): "
        f"{stats.get('hash_groups', 0)} hash, "
        f"{stats.get('rev_comp_groups', 0)} rev-comp "
        f"({stats.get('sequences', 0)} sequences scanned). "
        "Report only — run cleanup_accessions.py to apply."
        if n
        else f"No hash/rev-comp consolidations among {stats.get('sequences', 0)} sequences."
    )
    return {
        "job": "accession_cleanup_analyze",
        "mode": "report",
        "columns": _ACCESSION_CLEANUP_COLS,
        "rows": preview,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_taxonomy_ncbi_backfill():
    from src.database.cleanup.utils.fill_taxonomy_from_ncbi import (
        fill_taxonomy_from_ncbi,
    )

    summary, rows, _counts = fill_taxonomy_from_ncbi(dry_run=False, overwrite=False)
    return {
        "job": "taxonomy_ncbi_backfill",
        "mode": "apply",
        "columns": _TAX_NCBI_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_genome_coordinates_validate():
    """Check contig/chr IDs and coordinates against linked assembly accessions."""
    from src.database.cleanup.utils.validate_genome_coordinates import (
        analyze_genome_coordinates,
    )

    rows, stats = analyze_genome_coordinates(
        include_submissions=True,
        validate_ncbi=True,
        validate_sequences=True,
    )
    n = len(rows)
    by_type: dict[str, int] = {}
    for r in rows:
        t = r.get("issue_type") or "unknown"
        by_type[t] = by_type.get(t, 0) + 1
    type_bits = ", ".join(f"{k}={v}" for k, v in sorted(by_type.items())[:6])
    jgi_skip = stats.get("jgi_rows_skipped", 0)
    seq_ok = stats.get("sequence_matches", 0)
    summary = (
        f"{n} issue(s) across {stats.get('rows_checked', 0)} row(s)"
        + (f": {type_bits}" if type_bits else "")
        + (
            f"; NCBI metadata on {stats.get('ncbi_rows_checked', 0)} row(s) "
            f"({stats.get('ncbi_assemblies', 0)} assemblies)"
            if stats.get("ncbi_rows_checked")
            else ""
        )
        + (
            f"; {seq_ok} full sequence match(es), "
            f"{stats.get('sequence_mismatches', 0)} mismatch(es)"
            if seq_ok or stats.get("sequence_mismatches")
            else ""
        )
        + (f"; {jgi_skip} JGI row(s) skipped external lookup" if jgi_skip else "")
        if n or jgi_skip or seq_ok
        else f"No genome coordinate issues in {stats.get('rows_checked', 0)} row(s)."
    )
    return {
        "job": "genome_coordinates_validate",
        "mode": "report",
        "columns": _GENOME_COORD_COLS,
        "rows": rows,
        "proposed_changes": [],
        "summary": summary,
    }


def _job_genome_coordinates_fix():
    """Align ship sequences to NCBI contigs and stage coordinate corrections."""
    from src.database.cleanup.utils.validate_genome_coordinates import (
        analyze_genome_coordinate_fixes,
    )

    proposed, preview, stats = analyze_genome_coordinate_fixes(require_perfect=True)
    n = len(proposed)
    field_updates = len({(c["row_id"], c["col_id"]) for c in proposed})
    staged_ships = stats.get("fixes_staged", 0)
    summary = (
        f"{field_updates} field update(s) staged for {staged_ships} ship(s)"
        + (
            f"; {stats.get('already_correct', 0)} already match GenBank, "
            f"{stats.get('partial_match', 0)} partial alignment(s), "
            f"{stats.get('no_alignment', 0)} not located"
            if stats.get("rows_checked")
            else ""
        )
        if n
        else (
            f"No coordinate fixes to stage ({stats.get('rows_checked', 0)} checked, "
            f"{stats.get('already_correct', 0)} already correct, "
            f"{stats.get('no_alignment', 0)} not located)."
        )
    )
    return {
        "job": "genome_coordinates_fix",
        "mode": "stage",
        "columns": _GENOME_COORD_FIX_COLS,
        "rows": preview,
        "proposed_changes": proposed,
        "summary": summary,
    }


JOB_REGISTRY = {
    "md5_validation": {
        "label": "MD5 Duplicate Detection",
        "description": "Find ships with duplicate or missing MD5 / rev_comp_md5 hashes.",
        "fn": _job_md5_validation,
        "mode": "report",
    },
    "backfill_rev_comp_md5": {
        "label": "Backfill missing rev_comp_md5",
        "description": "Compute rev_comp_md5 for ships where it is NULL (does not overwrite existing values).",
        "fn": _job_backfill_rev_comp_md5,
        "mode": "stage",
    },
    "taxonomy_validate": {
        "label": "Taxonomy consistency check",
        "description": "Report duplicates, orphaned taxonomy rows, and family inconsistencies.",
        "fn": _job_taxonomy_validate,
        "mode": "report",
    },
    "taxonomy_clean": {
        "label": "Clean taxonomy fields",
        "description": "Stage whitespace, species, strain, and name fixes (TaxonomyConsolidator).",
        "fn": _job_taxonomy_clean,
        "mode": "stage",
    },
    "taxonomy_hierarchy_backfill": {
        "label": "Backfill taxonomy hierarchy",
        "description": "Fill empty hierarchy fields from other rows with same family/genus.",
        "fn": _job_taxonomy_hierarchy_backfill,
        "mode": "stage",
    },
    "taxonomy_link_ships": {
        "label": "Link missing ship taxonomy",
        "description": "Stage joined_ships.tax_id links from genome, ome, or gluck_thaler_2025 rules.",
        "fn": _job_taxonomy_link_ships,
        "mode": "stage",
    },
    "taxonomy_ncbi_backfill": {
        "label": "Fill taxonomy from NCBI",
        "description": (
            "Look up NCBI taxon from GCA/GCF assembly accessions or SAMN/SAMEA biosamples, "
            "create taxonomy rows if needed, and set genomes.taxonomy_id and joined_ships.tax_id "
            "(writes immediately; skips rows that already have taxonomy; requires NCBI_API_KEY in .env)."
        ),
        "fn": _job_taxonomy_ncbi_backfill,
        "mode": "apply",
    },
    "accession_status": {
        "label": "Accession status",
        "description": "List ships with sequences but no accession_id.",
        "fn": _job_accession_status,
        "mode": "report",
    },
    "accession_diagnose": {
        "label": "Diagnose accession sync",
        "description": "Report out-of-sync joined_ships, FK mismatches, duplicates, and orphans.",
        "fn": _job_accession_diagnose,
        "mode": "report",
    },
    "accession_sync": {
        "label": "Sync joined_ships accession_id",
        "description": "Stage joined_ships.accession_id updates from ships.accession_id.",
        "fn": _job_accession_sync,
        "mode": "stage",
    },
    "accession_assign_preview": {
        "label": "Preview missing accession assignment",
        "description": (
            f"Dry-run AccessionManager pipeline (max {_ACCESSION_ASSIGN_PREVIEW_LIMIT} ships). "
            "Apply via accession_manager.py --fix-missing --apply."
        ),
        "fn": _job_accession_assign_preview,
        "mode": "report",
    },
    "accession_cleanup_analyze": {
        "label": "Analyze accession consolidations",
        "description": "Report hash-duplicate and reverse-complement groups (fast scan; apply via cleanup_accessions CLI).",
        "fn": _job_accession_cleanup_analyze,
        "mode": "report",
    },
    "accession_fix_duplicates_preview": {
        "label": "Preview duplicate ship_id cleanup",
        "description": "Show joined_ships rows that would be cleared when starshipID != ships.header.",
        "fn": _job_accession_fix_duplicates_preview,
        "mode": "report",
    },
    "genome_coordinates_validate": {
        "label": "Validate genome coordinates",
        "description": (
            "Check contig/chr IDs and coordinates against genomes.assembly_accession. "
            "NCBI GCA/GCF: validates contig/range via Datasets, then compares the "
            "full ship sequence to the GenBank interval (set NCBI_API_KEY in .env). "
            "Includes pending submissions."
        ),
        "fn": _job_genome_coordinates_validate,
        "mode": "report",
    },
    "genome_coordinates_fix": {
        "label": "Fix genome coordinates from sequence",
        "description": (
            "Align ship sequences to NCBI contigs (minimap2, same thresholds as "
            "check_contained_match). Stages perfect-match updates to "
            "starship_features elementBegin/elementEnd/elementLength/strand only — "
            "never modifies ship sequences. Run Validate first to review mismatches."
        ),
        "fn": _job_genome_coordinates_fix,
        "mode": "stage",
    },
}


def _merge_pending_changes(existing, new_changes):
    """Merge job/manual pending entries; later entry wins per table/row/col."""
    pending = list(existing or [])
    for ch in new_changes:
        entry = {k: v for k, v in ch.items() if k != "note"}
        if "source" not in entry:
            entry["source"] = "job"
        idx = next(
            (
                i
                for i, p in enumerate(pending)
                if p.get("table") == entry["table"]
                and p.get("row_id") == entry["row_id"]
                and p.get("col_id") == entry["col_id"]
            ),
            -1,
        )
        if idx >= 0:
            pending[idx] = entry
        else:
            pending.append(entry)
    return pending


def _overlay_job_changes_on_rowdata(table_key, changes):
    """Re-fetch rowData and overlay staged job values with _dirty='job'."""
    rows = _refetch_rowdata(table_key)
    if rows is no_update:
        return no_update
    assert isinstance(rows, list)
    table_changes = [c for c in changes if c.get("table") == table_key]
    if not table_changes:
        return no_update
    by_row = {}
    for ch in table_changes:
        by_row.setdefault(ch["row_id"], []).append(ch)
    out = []
    for row in rows:
        rid = row.get("id")
        nr = {**row, "_dirty": row.get("_dirty", False)}
        if rid in by_row:
            for ch in by_row[rid]:
                nr[ch["col_id"]] = ch.get("new_value", "")
            nr["_dirty"] = "job"
        out.append(nr)
    return out


# ---------------------------------------------------------------------------
# Promote helper
# ---------------------------------------------------------------------------


def _promote_submission(sub_id: int):
    """
    Load a staging submission and insert it into the main starbase DB via
    perform_database_insertion (SubmissionProcessor).

    Returns (success: bool, accession: str|None, error: str|None).
    """
    from src.utils.web_submission_adapter import (
        perform_database_insertion,
        parse_submission_fasta,
    )

    try:
        with get_submissions_session() as session:
            row = session.execute(
                text("SELECT * FROM submissions WHERE id = :id"), {"id": sub_id}
            ).fetchone()

        if not row:
            return False, None, f"Submission {sub_id} not found"

        row = dict(row._mapping)

        if row.get("processing_status") != "processed":
            return (
                False,
                None,
                f"Submission {sub_id} must be processed before promotion "
                f"(current status: {row.get('processing_status') or 'pending'})",
            )

        if not row.get("seq_contents"):
            return False, None, "No sequence data stored in this submission"

        staging_accession = str(row.get("accession_tag") or "").strip()
        if not staging_accession:
            return (
                False,
                None,
                f"Submission {sub_id} has no accession_tag; run Process first",
            )

        seq_id, sequence = parse_submission_fasta(row["seq_contents"])

        strand = row.get("shipstrand") or "+"
        processed_data = {
            "sequence": sequence,
            "starshipID": seq_id,
            "evidence": row.get("evidence") or "",
            "source": "web_submission_promoted",
            "genus": row.get("genus") or "",
            "species": row.get("species") or "",
            "strain": row.get("strain") or None,
            "contig_id": row.get("hostchr") or "",
            "assembly_accession": row.get("assembly_accession") or None,
            "element_start": row.get("shipstart"),
            "element_end": row.get("shipend"),
            "element_strand": strand,
            "curator": row.get("uploader") or "",
            "curated_status": "curated",
            "notes": row.get("comment") or "",
            "trust_staging": True,
            "staging_accession_tag": staging_accession,
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
                text(
                    """
                    UPDATE submissions
                    SET needs_review = 0, processing_status = 'promoted'
                    WHERE id = :id
                    """
                ),
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


def _make_grid(df, grid_id, editable_cols, row_selection=False):
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

    dash_grid_options = {
        "pagination": True,
        "paginationPageSize": 50,
        "suppressPropertyNamesCheck": True,
        "rowHeight": 40,
        "headerHeight": 44,
        "rowClassRules": {
            "admin-unsaved-row": "params.data._dirty === 'manual' || params.data._dirty === true",
            "admin-job-row": "params.data._dirty === 'job'",
        },
    }
    if row_selection:
        dash_grid_options["rowSelection"] = "multiple"

    return dag.AgGrid(
        id=grid_id,
        columnDefs=col_defs,
        rowData=rows,
        defaultColDef={"resizable": True, "minWidth": 80},
        dashGridOptions=dash_grid_options,
        getRowId="params.data.id",
        className="ag-theme-alpine",
        style={"width": "100%", "height": "600px"},
    )


def _make_result_table(columns, rows, table_id):
    """Read-only AG Grid for job report / diff preview."""
    if not rows:
        return dmc.Text("No rows.", size="sm", c="dimmed")
    col_defs = [
        {
            "field": c,
            "headerName": c,
            "filter": True,
            "sortable": True,
            "resizable": True,
            "minWidth": 90,
        }
        for c in columns
    ]
    return dag.AgGrid(
        id=table_id,
        columnDefs=col_defs,
        rowData=[{k: ("" if v is None else v) for k, v in r.items()} for r in rows],
        defaultColDef={"resizable": True, "minWidth": 80},
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": 25,
            "suppressPropertyNamesCheck": True,
        },
        className="ag-theme-alpine",
        style={"width": "100%", "height": "360px"},
    )


def _build_jobs_tab():
    job_options = [
        {"label": meta["label"], "value": key} for key, meta in JOB_REGISTRY.items()
    ]
    default_job = job_options[0]["value"] if job_options else None
    default_desc = JOB_REGISTRY[default_job]["description"] if default_job else ""

    return html.Div(
        [
            dmc.Text(
                "Run consistency checks or stage cleanup changes. "
                "Staged changes merge into the pending-changes queue — review and Save from the toolbar.",
                size="sm",
                c="dimmed",
                mb="sm",
                mt="xs",
            ),
            dmc.Paper(
                dmc.Stack(
                    [
                        dmc.Select(
                            id="admin-job-select",
                            label="Job",
                            data=job_options,
                            value=default_job,
                            allowDeselect=False,
                        ),
                        dmc.Text(
                            id="admin-job-description",
                            size="sm",
                            c="dimmed",
                            children=default_desc,
                        ),
                        dmc.Group(
                            [
                                dmc.Button(
                                    "Run job",
                                    id="admin-run-job-btn",
                                    color="blue",
                                    size="sm",
                                ),
                                dmc.Button(
                                    "Apply to pending",
                                    id="admin-apply-job-btn",
                                    color="teal",
                                    size="sm",
                                    variant="light",
                                    disabled=True,
                                    style={"display": "none"},
                                ),
                            ],
                            gap="xs",
                        ),
                    ],
                    gap="sm",
                ),
                p="md",
                withBorder=True,
                radius="sm",
                mb="md",
            ),
            html.Div(id="admin-job-result-panel"),
        ]
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


def _process_modal():
    return dmc.Modal(
        id="admin-process-modal",
        title="Processing submissions",
        centered=True,
        opened=False,
        closeOnClickOutside=False,
        withCloseButton=False,
        children=dmc.Stack(
            [
                dmc.Loader(type="dots", size="lg"),
                dmc.Text(
                    id="admin-process-modal-text",
                    children="Running checks, classification, and accession assignment…",
                    size="sm",
                ),
            ],
            align="center",
            gap="md",
            py="md",
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
                                    "Add row",
                                    id="admin-add-row-btn",
                                    variant="light",
                                    color="blue",
                                    size="sm",
                                ),
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
            dbc.Tabs(  # pyright: ignore[reportCallIssue]
                id="admin-tabs",
                children=[
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            js_df,
                            "admin-joined-ships-grid",
                            EDITABLE_COLS["joined_ships"],
                        ),
                        label=f"Starships ({len(js_df):,})",
                        tab_id="joined_ships",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        html.Div(
                            [
                                dmc.Group(
                                    [
                                        dmc.Text(
                                            "Select row(s), then Process (classify + accession) or Promote (main DB). "
                                            "Re-processing skips rows that already have an accession and classification.",
                                            size="sm",
                                            c="dimmed",
                                        ),
                                        dmc.Group(
                                            [
                                                dmc.Button(
                                                    "Process submission(s)",
                                                    id="admin-process-btn",
                                                    color="blue",
                                                    size="sm",
                                                    disabled=True,
                                                ),
                                                dmc.Button(
                                                    "Promote to Main DB",
                                                    id="admin-promote-btn",
                                                    color="green",
                                                    size="sm",
                                                    disabled=True,
                                                ),
                                            ],
                                            gap="sm",
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
                                    row_selection=True,
                                ),
                            ]
                        ),
                        label=f"Submissions ({len(sub_df):,})",
                        tab_id="submissions",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            tax_df,
                            "admin-taxonomy-grid",
                            EDITABLE_COLS["taxonomy"],
                        ),
                        label=f"Taxonomy ({len(tax_df):,})",
                        tab_id="taxonomy",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            pap_df,
                            "admin-papers-grid",
                            EDITABLE_COLS["papers"],
                        ),
                        label=f"Papers ({len(pap_df):,})",
                        tab_id="papers",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            fam_df,
                            "admin-family-names-grid",
                            EDITABLE_COLS["family_names"],
                        ),
                        label=f"Families ({len(fam_df):,})",
                        tab_id="family_names",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            nav_df,
                            "admin-navis-names-grid",
                            EDITABLE_COLS["navis_names"],
                        ),
                        label=f"Navis ({len(nav_df):,})",
                        tab_id="navis_names",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            hap_df,
                            "admin-haplotype-names-grid",
                            EDITABLE_COLS["haplotype_names"],
                        ),
                        label=f"Haplotypes ({len(hap_df):,})",
                        tab_id="haplotype_names",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            acc_df,
                            "admin-accessions-grid",
                            EDITABLE_COLS["accessions"],
                        ),
                        label=f"Accessions ({len(acc_df):,})",
                        tab_id="accessions",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            sacc_df,
                            "admin-ship-accessions-grid",
                            EDITABLE_COLS["ship_accessions"],
                        ),
                        label=f"Ship Accessions ({len(sacc_df):,})",
                        tab_id="ship_accessions",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _make_grid(
                            gen_df,
                            "admin-genomes-grid",
                            EDITABLE_COLS["genomes"],
                        ),
                        label=f"Genomes ({len(gen_df):,})",
                        tab_id="genomes",
                    ),
                    dbc.Tab(  # pyright: ignore[reportCallIssue]
                        _build_jobs_tab(),
                        label="Jobs",
                        tab_id="jobs",
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
        dcc.Store(id="admin-submissions-rowdata", data=None),
        dcc.Store(id="admin-job-result", data=None),
        dcc.Store(id="admin-job-applied", data=False),
        dcc.Store(id="admin-busy", data=False),
        dcc.Store(id="admin-job-running", data=False),
        dcc.Store(id="admin-busy-reset", data=0),
        dcc.Store(id="admin-process-meta", data=None),
        dcc.Store(id="admin-authed", data=False),
        dcc.Store(id="admin-ui-state", data={}),
        dcc.Store(
            id="admin-apply-job-meta",
            data={
                "disabled": True,
                "style": {"display": "none"},
                "children": "Apply to pending",
            },
        ),
        html.Div(id="admin-content"),
        _version_bump_modal(),
        _promote_modal(),
        _process_modal(),
    ]
)

# ---------------------------------------------------------------------------
# Callbacks — token gate
# ---------------------------------------------------------------------------


@callback(
    [
        Output("admin-content", "children"),
        Output("admin-authed", "data"),
    ],
    Input("admin-url", "search"),
)
def render_admin_content(search):
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
    return _build_admin_layout(), True


# ---------------------------------------------------------------------------
# Clientside — mark busy on action button click (cleared by Python when done)
# ---------------------------------------------------------------------------

clientside_callback(
    """
    function() {
        return true;
    }
    """,
    Output("admin-busy", "data"),
    Input("admin-apply-job-btn", "n_clicks"),
    Input("admin-save-btn", "n_clicks"),
    Input("admin-discard-btn", "n_clicks"),
    Input("admin-process-btn", "n_clicks"),
    Input("admin-promote-confirm", "n_clicks"),
    Input("admin-version-confirm", "n_clicks"),
    prevent_initial_call=True,
)

clientside_callback(
    """
    function(n, selected) {
        if (!n) return [
            window.dash_clientside.no_update,
            window.dash_clientside.no_update,
        ];
        var count = (selected || []).length;
        var label = count === 1
            ? 'Processing 1 submission…'
            : ('Processing ' + count + ' submissions…');
        return [true, label];
    }
    """,
    Output("admin-process-modal", "opened"),
    Output("admin-process-modal-text", "children"),
    Input("admin-process-btn", "n_clicks"),
    State("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)

# Run job blocks synchronously on the server — set busy immediately in the browser
# so the button disables before the round-trip completes.
clientside_callback(
    """
    function(n) {
        if (!n) return [
            window.dash_clientside.no_update,
            window.dash_clientside.no_update,
            window.dash_clientside.no_update,
            window.dash_clientside.no_update,
        ];
        return [true, true, true, 'Running job...'];
    }
    """,
    Output("admin-job-running", "data"),
    Output("admin-run-job-btn", "disabled", allow_duplicate=True),
    Output("admin-run-job-btn", "loading", allow_duplicate=True),
    Output("admin-run-job-btn", "children", allow_duplicate=True),
    Input("admin-run-job-btn", "n_clicks"),
    prevent_initial_call=True,
)

# ---------------------------------------------------------------------------
# Clientside — compute UI state into static store (safe before auth layout loads)
# ---------------------------------------------------------------------------

clientside_callback(
    """
    function(busy, jobRunning, pending, applyMeta, selected) {
        var isBusy = !!busy || !!jobRunning;
        var n = (pending || []).length;
        var meta = applyMeta || {
            disabled: true,
            style: {display: 'none'},
            children: 'Apply to pending'
        };
        return {
            isBusy: isBusy,
            saveDis: isBusy || n === 0,
            discardDis: isBusy || n === 0,
            saveLbl: isBusy ? 'Saving...'
                : (n > 0 ? ('Save ' + n + ' change' + (n !== 1 ? 's' : '')) : 'Save changes'),
            hint: isBusy ? 'Working — please wait...'
                : (n > 0 ? (n + ' unsaved change' + (n !== 1 ? 's' : '') + ' \\u2014 click Save or Discard.')
                         : 'No pending changes. Yellow cells are editable \\u2014 edit freely, then save. Blue rows = job-staged.'),
            applyDis: isBusy || !!meta.disabled,
            applyStyle: meta.style || {display: 'none'},
            applyLbl: isBusy ? 'Applying...' : (meta.children || 'Apply to pending'),
            promoteDis: isBusy || !(selected && selected.length),
            processDis: isBusy || !(selected && selected.length),
            jobSelectDis: isBusy,
            runJobDis: isBusy,
            runJobLbl: jobRunning ? 'Running job...' : 'Run job',
        };
    }
    """,
    Output("admin-ui-state", "data"),
    Input("admin-busy", "data"),
    Input("admin-job-running", "data"),
    Input("admin-pending-changes", "data"),
    Input("admin-apply-job-meta", "data"),
    Input("admin-selected-submissions", "data"),
)


# Single duplicate writer clears busy flags after server actions complete.
@callback(
    [
        Output("admin-busy", "data", allow_duplicate=True),
        Output("admin-job-running", "data", allow_duplicate=True),
    ],
    Input("admin-busy-reset", "data"),
    prevent_initial_call=True,
)
def _clear_admin_busy(_token):
    return False, False


@callback(
    Output("admin-process-modal", "opened", allow_duplicate=True),
    Input("admin-busy-reset", "data"),
    prevent_initial_call=True,
)
def _close_process_modal(_token):
    return False


# ---------------------------------------------------------------------------
# Clientside — apply UI state to auth-gated components (only exist after login)
# ---------------------------------------------------------------------------

clientside_callback(
    """
    function(uiState, authed, content) {
        var nu = window.dash_clientside.no_update;
        if (!authed || !content || !uiState || uiState.saveDis === undefined) {
            return [nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu, nu];
        }
        var b = uiState.isBusy;
        return [
            uiState.saveDis, uiState.discardDis, uiState.saveLbl, uiState.hint,
            b, b,
            uiState.runJobDis, b, uiState.runJobLbl,
            uiState.applyDis, uiState.applyStyle, uiState.applyLbl, b,
            uiState.processDis, b,
            uiState.promoteDis, b,
            uiState.jobSelectDis,
        ];
    }
    """,
    Output("admin-save-btn", "disabled"),
    Output("admin-discard-btn", "disabled"),
    Output("admin-save-btn", "children"),
    Output("admin-pending-label", "children"),
    Output("admin-save-btn", "loading"),
    Output("admin-discard-btn", "loading"),
    Output("admin-run-job-btn", "disabled"),
    Output("admin-run-job-btn", "loading"),
    Output("admin-run-job-btn", "children"),
    Output("admin-apply-job-btn", "disabled"),
    Output("admin-apply-job-btn", "style"),
    Output("admin-apply-job-btn", "children"),
    Output("admin-apply-job-btn", "loading"),
    Output("admin-process-btn", "disabled"),
    Output("admin-process-btn", "loading"),
    Output("admin-promote-btn", "disabled"),
    Output("admin-promote-btn", "loading"),
    Output("admin-job-select", "disabled"),
    Input("admin-ui-state", "data"),
    Input("admin-authed", "data"),
    Input("admin-content", "children"),
    prevent_initial_call=True,
)

# ---------------------------------------------------------------------------
# Clientside — modal confirm buttons live in static layout; always safe to update
# ---------------------------------------------------------------------------

clientside_callback(
    """
    function(uiState) {
        var nu = window.dash_clientside.no_update;
        if (!uiState || uiState.isBusy === undefined) return [nu, nu, nu, nu];
        var b = uiState.isBusy;
        return [b, b, b, b];
    }
    """,
    Output("admin-promote-confirm", "disabled"),
    Output("admin-promote-confirm", "loading"),
    Output("admin-version-confirm", "disabled"),
    Output("admin-version-confirm", "loading"),
    Input("admin-ui-state", "data"),
)

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
    if (idx >= 0) {{ p[idx].new_value = nv; p[idx].source = 'manual'; }} else {{ p.push({{table:tk,row_id:rid,col_id:cid,old_value:ov,new_value:nv,source:'manual'}}); }}
    var rd = (rowData||[]).map(function(r) {{ return r.id===rid ? Object.assign({{}},r,{{_dirty:'manual'}}) : r; }});
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
# Clientside callbacks — "Add row" button inserts a blank editable row into
# whichever grid is on the active tab, and stages a __new_row__ marker so
# save_all_pending knows to INSERT rather than UPDATE for that row_id.
# ---------------------------------------------------------------------------

_ADD_ROW_JS = """
function(n, activeTab, rowData) {{
    if (!n || activeTab !== "{tk}") return window.dash_clientside.no_update;
    var tempId = -n;
    var editableCols = {editable_cols_json};
    var newRow = {{id: tempId, _dirty: 'new'}};
    for (var i = 0; i < editableCols.length; i++) {{ newRow[editableCols[i]] = ''; }}
    var rd = (rowData || []).slice();
    rd.unshift(newRow);
    return rd;
}}
"""

for _tk, _gid in GRID_IDS.items():
    clientside_callback(
        _ADD_ROW_JS.format(
            tk=_tk, editable_cols_json=list(EDITABLE_COLS.get(_tk, set()))
        ),
        Output(_gid, "rowData", allow_duplicate=True),
        Input("admin-add-row-btn", "n_clicks"),
        State("admin-tabs", "active_tab"),
        State(_gid, "rowData"),
        prevent_initial_call=True,
    )

# tempId must match the -n_clicks value used above so the marker's row_id
# lines up with the blank row's id for whichever grid actually added one.
clientside_callback(
    """
    function(n, activeTab, pending) {
        if (!n || !activeTab) return window.dash_clientside.no_update;
        var tempId = -n;
        var p = JSON.parse(JSON.stringify(pending || []));
        p.push({
            table: activeTab, row_id: tempId, col_id: "__new_row__",
            old_value: null, new_value: {}, source: 'manual'
        });
        return p;
    }
    """,
    Output("admin-pending-changes", "data", allow_duplicate=True),
    Input("admin-add-row-btn", "n_clicks"),
    State("admin-tabs", "active_tab"),
    State("admin-pending-changes", "data"),
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

clientside_callback(
    """
    function(rowData) {
        if (rowData === null || rowData === undefined) {
            return window.dash_clientside.no_update;
        }
        return rowData;
    }
    """,
    Output("admin-submissions-grid", "rowData", allow_duplicate=True),
    Input("admin-submissions-rowdata", "data"),
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
        Output("admin-busy-reset", "data", allow_duplicate=True),
    ],
    Input("admin-save-btn", "n_clicks"),
    State("admin-pending-changes", "data"),
    prevent_initial_call=True,
)
def save_all_pending(n_clicks, pending):
    n_out = len(GRID_IDS) + 6  # pending + grids + modal×3 + notif + busy-reset
    if not n_clicks or not pending:
        return [no_update] * (n_out - 1) + [time.time()]

    logger.info("save_all_pending: %d change(s): %s", len(pending), pending)
    errors, successes = [], []

    new_row_markers = [ch for ch in pending if ch.get("col_id") == "__new_row__"]
    new_row_keys = {(m["table"], m["row_id"]) for m in new_row_markers}
    updates = [
        ch
        for ch in pending
        if (ch["table"], ch["row_id"]) not in new_row_keys
    ]

    for table, row_id in new_row_keys:
        col_values = {
            ch["col_id"]: ch.get("new_value")
            for ch in pending
            if ch["table"] == table
            and ch["row_id"] == row_id
            and ch["col_id"] != "__new_row__"
        }
        ok, err, new_id = _do_insert(table, col_values)
        marker = {"table": table, "row_id": row_id, "col_id": "__new_row__"}
        (successes if ok else errors).append((marker, err))

    for ch in updates:
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
        return [[], *fresh, False, no_update, no_update, notif, time.time()]

    info, desc = _make_version_modal_content([c for c, _ in successes])
    return [[], *fresh, True, info, desc, no_update, time.time()]


@callback(
    [
        Output("admin-pending-changes", "data", allow_duplicate=True),
        *[Output(gid, "rowData", allow_duplicate=True) for gid in GRID_IDS.values()],
        Output("admin-busy-reset", "data", allow_duplicate=True),
    ],
    Input("admin-discard-btn", "n_clicks"),
    prevent_initial_call=True,
)
def discard_pending(n_clicks):
    if not n_clicks:
        return [no_update] * (len(GRID_IDS) + 1) + [time.time()]
    return [[], *[_refetch_rowdata(k) for k in GRID_IDS], time.time()]


# ---------------------------------------------------------------------------
# Callbacks — version modal confirm / cancel
# ---------------------------------------------------------------------------


@callback(
    [
        Output("admin-version-modal", "opened", allow_duplicate=True),
        Output("notifications-container", "children", allow_duplicate=True),
        Output("admin-busy-reset", "data", allow_duplicate=True),
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
        return no_update, no_update, no_update

    if bump_type == "skip":
        notification = dmc.Notification(
            id=f"admin-ok-{uuid.uuid4().hex[:6]}",
            title="Saved",
            message="Change saved. Version not bumped.",
            color="green",
            action="show",
            autoClose=4000,
        )
        return False, notification, time.time()

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
    return False, notification, time.time()


@callback(
    Output("admin-version-modal", "opened", allow_duplicate=True),
    Input("admin-version-cancel", "n_clicks"),
    prevent_initial_call=True,
)
def cancel_version_bump(n_clicks):
    return False


# ---------------------------------------------------------------------------
# Callbacks — process submission
# ---------------------------------------------------------------------------


@callback(
    Output("notifications-container", "children", allow_duplicate=True),
    Output("admin-submissions-rowdata", "data", allow_duplicate=True),
    Output("admin-busy-reset", "data", allow_duplicate=True),
    Input("admin-process-btn", "n_clicks"),
    State("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)
def run_process(n_clicks, selected_rows):
    if not n_clicks or not selected_rows:
        raise PreventUpdate

    from src.utils.web_submission_adapter import (
        process_staging_submissions_ordered,
        summarize_staging_process_results,
    )

    sub_ids = sorted({row.get("id") for row in selected_rows if row.get("id")})
    results = process_staging_submissions_ordered(sub_ids)
    grid_rows = _fetch_submissions().fillna("").to_dict("records")

    failures = [r for r in results if not r.get("success")]
    detail = summarize_staging_process_results(results)
    if failures:
        notif = dmc.Notification(
            id=f"admin-err-{uuid.uuid4().hex[:6]}",
            title=f"Processing finished with {len(failures)} error(s)",
            message=detail,
            color="red",
            action="show",
            autoClose=15000,
        )
        return notif, grid_rows, time.time()

    skipped = sum(1 for r in results if r.get("already_processed"))
    classified = sum(1 for r in results if r.get("classified"))
    notif = dmc.Notification(
        id=f"admin-ok-{uuid.uuid4().hex[:6]}",
        title="Processing complete",
        message=detail
        or (
            f"{len(results)} submission(s): "
            f"{classified} classified, {skipped} skipped (already done)"
        ),
        color="green",
        action="show",
        autoClose=12000,
    )
    return notif, grid_rows, time.time()


# ---------------------------------------------------------------------------
# Callbacks — promote submission
# ---------------------------------------------------------------------------


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
                                    row.get("processing_status") or "pending",
                                    color="blue"
                                    if row.get("processing_status") == "processed"
                                    else "gray",
                                    variant="light",
                                    size="xs",
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
        Output("admin-busy-reset", "data", allow_duplicate=True),
        Output("admin-submissions-rowdata", "data", allow_duplicate=True),
    ],
    Input("admin-promote-confirm", "n_clicks"),
    State("admin-selected-submissions", "data"),
    prevent_initial_call=True,
)
def run_promotion(n_clicks, selected_rows):
    if not n_clicks or not selected_rows:
        return (no_update,) * 7

    results = []
    for row in selected_rows:
        sub_id = row.get("id")
        success, accession, error = _promote_submission(sub_id)
        results.append((sub_id, success, accession, error))

    grid_rows = _fetch_submissions().fillna("").to_dict("records")
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
        return False, False, no_update, no_update, notif, time.time(), grid_rows

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
    return False, True, info, desc, no_update, time.time(), grid_rows


@callback(
    Output("admin-promote-modal", "opened", allow_duplicate=True),
    Input("admin-promote-cancel", "n_clicks"),
    prevent_initial_call=True,
)
def cancel_promote(n_clicks):
    return False


# ---------------------------------------------------------------------------
# Callbacks — admin jobs
# ---------------------------------------------------------------------------


@callback(
    Output("admin-job-description", "children"),
    Input("admin-job-select", "value"),
)
def update_job_description(job_key):
    if not job_key or job_key not in JOB_REGISTRY:
        return ""
    return JOB_REGISTRY[job_key]["description"]


@callback(
    [
        Output("admin-job-result", "data"),
        Output("admin-job-applied", "data"),
        Output("notifications-container", "children", allow_duplicate=True),
        Output("admin-busy-reset", "data", allow_duplicate=True),
    ],
    Input("admin-run-job-btn", "n_clicks"),
    State("admin-job-select", "value"),
    prevent_initial_call=True,
)
def run_job(n_clicks, job_key):
    if not n_clicks or not job_key:
        return no_update, no_update, no_update, time.time()
    meta = JOB_REGISTRY.get(job_key)
    if not meta:
        return no_update, no_update, no_update, time.time()
    try:
        result = meta["fn"]()
        logger.info("admin job %s: %s", job_key, result.get("summary"))
        return result, False, no_update, time.time()
    except Exception as exc:
        logger.error("admin job %s failed: %s", job_key, exc)
        notif = dmc.Notification(
            id=f"admin-job-err-{uuid.uuid4().hex[:6]}",
            title="Job failed",
            message=str(exc),
            color="red",
            action="show",
            autoClose=10000,
        )
        return no_update, no_update, notif, time.time()


@callback(
    [
        Output("admin-job-result-panel", "children"),
        Output("admin-apply-job-meta", "data"),
    ],
    [
        Input("admin-job-result", "data"),
        Input("admin-job-applied", "data"),
        Input("admin-authed", "data"),
    ],
    prevent_initial_call=True,
)
def render_job_result(result, applied, authed):
    if not authed or not result:
        return no_update, no_update

    mode = result.get("mode", "report")
    summary = result.get("summary", "")
    rows = result.get("rows") or []
    columns = result.get("columns") or []
    proposed = result.get("proposed_changes") or []
    job_label = JOB_REGISTRY.get(result.get("job", ""), {}).get("label", "Job")

    parts = [
        dmc.Group(
            [
                dmc.Text(job_label, fw=600, size="sm"),
                dmc.Badge(
                    "Report only"
                    if mode == "report"
                    else ("Applied" if mode == "apply" else "Stage changes"),
                    color="gray"
                    if mode == "report"
                    else ("green" if mode == "apply" else "teal"),
                    variant="light",
                    size="sm",
                ),
            ],
            gap="xs",
        ),
        dmc.Text(summary, size="sm"),
    ]
    if mode == "apply":
        parts.append(
            dmc.Alert(
                "Changes were written to the database.",
                color="green",
                variant="light",
            )
        )
    elif mode == "stage":
        parts.append(
            dmc.Alert(
                "Review proposed changes below, then Apply to pending. "
                "Nothing is written until you Save from the toolbar.",
                color="blue",
                variant="light",
            )
        )
    if rows:
        parts.append(_make_result_table(columns, rows, "admin-job-result-grid"))
    else:
        parts.append(dmc.Text("No rows to display.", size="sm", c="dimmed"))

    show_apply = mode == "stage" and bool(proposed) and not applied
    apply_meta = {
        "disabled": not show_apply,
        "style": {} if show_apply else {"display": "none"},
        "children": f"Apply {len(proposed)} change(s) to pending",
    }

    return dmc.Stack(parts, gap="sm"), apply_meta


@callback(
    [
        Output("admin-pending-changes", "data", allow_duplicate=True),
        Output("admin-job-applied", "data", allow_duplicate=True),
        *[Output(gid, "rowData", allow_duplicate=True) for gid in GRID_IDS.values()],
        Output("notifications-container", "children", allow_duplicate=True),
        Output("admin-busy-reset", "data", allow_duplicate=True),
    ],
    Input("admin-apply-job-btn", "n_clicks"),
    [
        State("admin-job-result", "data"),
        State("admin-pending-changes", "data"),
    ],
    prevent_initial_call=True,
)
def apply_job_changes(n_clicks, result, pending):
    n_out = len(GRID_IDS) + 4
    if not n_clicks or not result:
        return [no_update] * (n_out - 1) + [time.time()]

    proposed = result.get("proposed_changes") or []
    if not proposed:
        return [no_update] * (n_out - 1) + [time.time()]

    merged = _merge_pending_changes(pending, proposed)
    changed_tables = {c["table"] for c in proposed}
    grid_updates = [
        _overlay_job_changes_on_rowdata(k, proposed)
        if k in changed_tables
        else no_update
        for k in GRID_IDS
    ]

    n = len(proposed)
    notif = dmc.Notification(
        id=f"admin-job-ok-{uuid.uuid4().hex[:6]}",
        title="Changes staged",
        message=f"{n} job change(s) added to pending queue. Review and Save.",
        color="blue",
        action="show",
        autoClose=6000,
    )
    return [merged, True, *grid_updates, notif, time.time()]

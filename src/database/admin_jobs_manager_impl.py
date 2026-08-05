"""Local implementations of admin consistency/cleanup jobs.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/admin.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.

Only the self-contained, mostly-dry-run jobs live here (12 of 15) --
taxonomy_ncbi_backfill (writes immediately, hits NCBI) and the two
genome_coordinates jobs (NCBI Datasets + minimap2, slower) stay in
src/pages/admin.py for now, deferred to a later slice.
"""

import pandas as pd

from src.database.sql_engine import get_starbase_session

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


_JOBS = {
    "taxonomy_validate": _job_taxonomy_validate,
    "taxonomy_clean": _job_taxonomy_clean,
    "taxonomy_hierarchy_backfill": _job_taxonomy_hierarchy_backfill,
    "taxonomy_link_ships": _job_taxonomy_link_ships,
    "md5_validation": _job_md5_validation,
    "backfill_rev_comp_md5": _job_backfill_rev_comp_md5,
    "accession_status": _job_accession_status,
    "accession_assign_preview": _job_accession_assign_preview,
    "accession_sync": _job_accession_sync,
    "accession_diagnose": _job_accession_diagnose,
    "accession_fix_duplicates_preview": _job_accession_fix_duplicates_preview,
    "accession_cleanup_analyze": _job_accession_cleanup_analyze,
}


def run_job(job_key: str) -> dict:
    if job_key not in _JOBS:
        raise ValueError(f"Unknown admin job_key: {job_key}")
    return _JOBS[job_key]()

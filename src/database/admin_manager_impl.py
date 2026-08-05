"""Local SQLAlchemy implementation of admin grid read/write operations.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/admin.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.
"""

import pandas as pd
from sqlalchemy import text

from src.config.logging import get_logger
from src.database.admin_table_config import (
    BOOLEAN_COLS,
    EDITABLE_COLS,
    TABLE_CONFIG,
    sql_col_ref,
)
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

# One DB (starbase.sqlite) holds every admin-editable table now, submissions
# included -- see alembic/versions/m6n7o8p9q0r1_*.
_QUERIES = {
    "joined_ships": """
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
    "submissions": """
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
    "taxonomy": """
        SELECT id, name, taxID, genus, species, strain, kingdom, phylum
        FROM taxonomy
        ORDER BY id
        """,
    "papers": """
        SELECT id, Title, Author, PublicationYear, DOI, Url,
               shortCitation, starshipMentioned, typePaper
        FROM papers
        ORDER BY id
        """,
    "family_names": """
        SELECT id, familyName, longFamilyID, oldFamilyID, otherFamilyID,
               clade, newFamilyID, type_element_reference, notes
        FROM family_names
        ORDER BY id
        """,
    "navis_names": """
        SELECT n.id, n.navis_name, n.previous_navis_name, n.activity,
               f.familyName
        FROM navis_names n
        LEFT JOIN family_names f ON n.ship_family_id = f.id
        ORDER BY n.id
        """,
    "haplotype_names": """
        SELECT h.id, h.haplotype_name, h.previous_haplotype_name, h.activity,
               n.navis_name, f.familyName
        FROM haplotype_names h
        LEFT JOIN navis_names n ON h.navis_id = n.id
        LEFT JOIN family_names f ON h.ship_family_id = f.id
        ORDER BY h.id
        """,
    "accessions": """
        SELECT a.id, a.accession_tag, a.version_tag
        FROM accessions a
        ORDER BY a.id
        """,
    "ship_accessions": """
        SELECT sa.id, sa.ship_accession_tag, sa.ship_accession_display,
               sa.ship_version_tag, sa.ship_id
        FROM ship_accessions sa
        ORDER BY sa.id
        """,
    "genomes": """
        SELECT g.id, g.ome, g.version, g.genomeSource, g.citation,
               g.biosample, g.acquisition_date, g.assembly_accession, g.taxonomy_id,
               t.name AS taxonomy_name
        FROM genomes g
        LEFT JOIN taxonomy t ON g.taxonomy_id = t.id
        ORDER BY g.id
        """,
}


def fetch_admin_table(table_key: str) -> pd.DataFrame:
    """Fetch rowData for one admin grid, keyed by table_key."""
    if table_key not in _QUERIES:
        raise ValueError(f"Unknown admin table_key: {table_key}")
    with get_starbase_session() as session:
        return pd.read_sql_query(_QUERIES[table_key], session.bind)


def do_insert(table_key: str, col_values: dict):
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

    col_refs = ", ".join(sql_col_ref(c) for c in cols)
    placeholders = ", ".join(f":{c}" for c in cols)
    sql = text(
        f"INSERT INTO {config['sql_table']} ({col_refs}) VALUES ({placeholders})"
    )

    try:
        with get_starbase_session() as session:
            result = session.execute(sql, values)
            session.commit()
            new_id = result.lastrowid
        logger.info("admin INSERT %s: %r (id=%s)", config["sql_table"], values, new_id)
        return True, None, new_id
    except Exception as exc:
        logger.error("admin INSERT error (%s): %s", table_key, exc)
        return False, str(exc), None


def do_update(table_key: str, row_id, col_id: str, new_value):
    """Run a whitelisted UPDATE. Returns (success: bool, error: str|None)."""
    allowed = EDITABLE_COLS.get(table_key, set())
    if col_id not in allowed:
        return False, f"Column '{col_id}' is not editable."

    config = TABLE_CONFIG[table_key]

    if col_id in BOOLEAN_COLS:
        new_value = 1 if str(new_value).lower() in ("true", "1", "yes") else 0

    sql = text(
        f"UPDATE {config['sql_table']} SET {sql_col_ref(col_id)} = :val WHERE {config['pk']} = :pk"
    )

    try:
        with get_starbase_session() as session:
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

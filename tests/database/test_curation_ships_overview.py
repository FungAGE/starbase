"""Read-only tests for the curation page's Ships overview query."""

from sqlalchemy import text

from src.database.curation_manager_impl import fetch_ships_overview
from src.database.sql_engine import get_starbase_session

_EXPECTED_KEYS = {
    "id",
    "starshipID",
    "curated_status",
    "evidence",
    "source",
    "ship_accession_id",
    "accession_tag",
    "ship_accession_tag",
    "familyName",
    "navis_name",
    "haplotype_name",
    "taxonomy_name",
    "genus",
    "species",
    "sequence_length",
}


def test_fetch_ships_overview_row_count_matches_source_table():
    with get_starbase_session() as session:
        expected_count = session.execute(
            text("SELECT COUNT(*) FROM joined_ships WHERE is_deleted = 0")
        ).scalar()

    rows = fetch_ships_overview()
    assert len(rows) == expected_count


def test_fetch_ships_overview_row_shape():
    rows = fetch_ships_overview()
    assert rows, "expected at least one row in the scratch DB"
    assert _EXPECTED_KEYS.issubset(rows[0].keys())


def test_fetch_ships_overview_is_read_only():
    with get_starbase_session() as session:
        before = session.execute(
            text("SELECT COUNT(*) FROM joined_ships")
        ).scalar()

    fetch_ships_overview()

    with get_starbase_session() as session:
        after = session.execute(
            text("SELECT COUNT(*) FROM joined_ships")
        ).scalar()

    assert before == after

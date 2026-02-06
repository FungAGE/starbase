"""add accession_display and ship_accession_display generated columns

Revision ID: f1a2b3c4d5e6
Revises: v0_1_0_pre
Create Date: 2026-02-06

Adds computed columns to accessions and ship_accessions so display strings
(accession.version) are stored in the database instead of computed in Python.
Requires SQLite 3.31+ for generated columns; downgrade requires SQLite 3.35+ for DROP COLUMN.
Idempotent: skips ADD COLUMN if column already exists (e.g. after a partial run).
"""
from typing import Sequence, Union

from alembic import op
from sqlalchemy import text
from sqlalchemy.exc import OperationalError

# revision identifiers, used by Alembic.
revision: str = "f1a2b3c4d5e6"
down_revision: Union[str, None] = "v0_1_0_pre"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None

# SQL expression matching _make_accession_display: show "accession.version" when both
# are set and accession has no dot; else accession only; else ''.
_ACCESSION_DISPLAY_EXPR = (
    "CASE WHEN (accession_tag IS NOT NULL AND trim(coalesce(accession_tag,''))<>'' "
    "AND version_tag IS NOT NULL AND trim(cast(version_tag AS text))<>'' "
    "AND instr(accession_tag,'.')=0) THEN accession_tag||'.'||cast(version_tag AS text) "
    "WHEN (accession_tag IS NOT NULL AND trim(coalesce(accession_tag,''))<>'') THEN trim(accession_tag) ELSE '' END"
)
# ship_accessions uses ship_version_tag (not version_tag) as the version column name
_SHIP_ACCESSION_DISPLAY_EXPR = (
    "CASE WHEN (ship_accession_tag IS NOT NULL AND trim(coalesce(ship_accession_tag,''))<>'' "
    "AND ship_version_tag IS NOT NULL AND trim(cast(ship_version_tag AS text))<>'' "
    "AND instr(ship_accession_tag,'.')=0) THEN ship_accession_tag||'.'||cast(ship_version_tag AS text) "
    "WHEN (ship_accession_tag IS NOT NULL AND trim(coalesce(ship_accession_tag,''))<>'') THEN trim(ship_accession_tag) ELSE '' END"
)


def _has_column(connection, table: str, column: str) -> bool:
    result = connection.execute(text(f"PRAGMA table_info({table})"))
    rows = result.fetchall()
    return any(len(row) > 1 and row[1] == column for row in rows)


def _add_column_if_missing(connection, table: str, column: str, alter_sql: str) -> None:
    if _has_column(connection, table, column):
        return
    try:
        op.execute(alter_sql)
    except OperationalError as e:
        if "duplicate column name" not in str(e).lower():
            raise
        # Column was added by a previous partial run (e.g. SQLite non-transactional DDL)


def upgrade() -> None:
    conn = op.get_bind()
    _add_column_if_missing(
        conn,
        "accessions",
        "accession_display",
        f"ALTER TABLE accessions ADD COLUMN accession_display TEXT "
        f"GENERATED ALWAYS AS ({_ACCESSION_DISPLAY_EXPR})",
    )
    _add_column_if_missing(
        conn,
        "ship_accessions",
        "ship_accession_display",
        f"ALTER TABLE ship_accessions ADD COLUMN ship_accession_display TEXT "
        f"GENERATED ALWAYS AS ({_SHIP_ACCESSION_DISPLAY_EXPR})",
    )


def downgrade() -> None:
    # Requires SQLite 3.35+ for DROP COLUMN; only drop if present
    conn = op.get_bind()
    if _has_column(conn, "accessions", "accession_display"):
        op.execute("ALTER TABLE accessions DROP COLUMN accession_display")
    if _has_column(conn, "ship_accessions", "ship_accession_display"):
        op.execute("ALTER TABLE ship_accessions DROP COLUMN ship_accession_display")

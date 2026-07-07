"""drop unique constraint on database_versions.semantic_version

Revision ID: h1i2j3k4l5m6
Revises: g7h8i9j0k1l2
Create Date: 2026-06-14

"""

from typing import Sequence, Union

from alembic import op

revision: str = "h1i2j3k4l5m6"
down_revision: Union[str, None] = "g7h8i9j0k1l2"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # SQLite does not support DROP CONSTRAINT, so recreate the table without
    # the UNIQUE constraint on semantic_version.
    op.execute("""
        CREATE TABLE database_versions_new (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            semantic_version TEXT NOT NULL,
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            created_by TEXT DEFAULT 'system'
        )
    """)
    op.execute(
        "INSERT INTO database_versions_new SELECT id, semantic_version, description, created_at, created_by FROM database_versions"
    )
    op.drop_table("database_versions")
    op.execute("ALTER TABLE database_versions_new RENAME TO database_versions")


def downgrade() -> None:
    op.execute("""
        CREATE TABLE database_versions_new (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            semantic_version TEXT NOT NULL UNIQUE,
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            created_by TEXT DEFAULT 'system'
        )
    """)
    # Keep only the latest row per version string when restoring the unique constraint.
    op.execute("""
        INSERT INTO database_versions_new
        SELECT id, semantic_version, description, created_at, created_by
        FROM database_versions
        GROUP BY semantic_version
        HAVING created_at = MAX(created_at)
    """)
    op.drop_table("database_versions")
    op.execute("ALTER TABLE database_versions_new RENAME TO database_versions")

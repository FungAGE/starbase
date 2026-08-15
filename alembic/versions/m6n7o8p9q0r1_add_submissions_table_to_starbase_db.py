"""add submissions table to starbase db (DB unification)

Submission has always lived on the same SQLAlchemy Base as everything
else in this file (see src/database/models/schema.py); it only ever ran
in a physically separate SQLite file because get_submissions_session()
pointed at a different engine. src/config/settings.py now points
DB_PATHS["submissions"] at the same file as DB_PATHS["starbase"], so this
creates the table here for the first time via Alembic.

Column list mirrors the live submissions.sqlite schema exactly (verified
against PRAGMA table_info) and src/database/models/schema.py's Submission
model -- no indexes/FKs on the original table, so none added here.

This migration only creates the table (idempotent via IF NOT EXISTS data
migration lives outside Alembic since it reads from a second SQLite file
-- see src/database/migrations/merge_submissions_into_starbase.py.

Revision ID: m6n7o8p9q0r1
Revises: l5m6n7o8p9q0
Create Date: 2026-08-05

"""

from typing import Sequence, Union

from alembic import op

revision: str = "m6n7o8p9q0r1"
down_revision: Union[str, None] = "l5m6n7o8p9q0"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.execute("""
        CREATE TABLE IF NOT EXISTS submissions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            seq_contents TEXT NOT NULL,
            seq_filename VARCHAR(255) NOT NULL,
            seq_date VARCHAR(50),
            anno_contents TEXT,
            anno_filename VARCHAR(255),
            anno_date VARCHAR(50),
            uploader VARCHAR(255),
            evidence VARCHAR(100),
            genus VARCHAR(255),
            species VARCHAR(255),
            strain VARCHAR(255),
            hostchr VARCHAR(255),
            assembly_accession VARCHAR(50),
            shipstart INTEGER,
            shipend INTEGER,
            shipstrand VARCHAR(10),
            comment TEXT,
            ship_accession_tag VARCHAR(50),
            accession_tag VARCHAR(50),
            needs_review BOOLEAN,
            classification_source VARCHAR(50),
            classification_family VARCHAR(100),
            classification_navis VARCHAR(100),
            classification_haplotype VARCHAR(100),
            closest_match VARCHAR(50),
            classification_confidence VARCHAR(20),
            submission_group_id VARCHAR(36),
            processing_status VARCHAR(20) DEFAULT 'pending'
        )
    """)


def downgrade() -> None:
    op.execute("DROP TABLE IF EXISTS submissions")

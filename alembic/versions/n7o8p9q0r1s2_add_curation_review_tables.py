"""add curation review tables (annotations, gene_features, annotation_history, result tables)

Phase 3: ported from MAS4starships' starship.Annotation / starship.Feature /
result_viewer models -- viewer/schema only, no data migration (starting
fresh per plan). Result-generation pipeline is a deferred slice; the 4
result tables exist so the viewer has somewhere to read from.

gene_features is intentionally not named starship_features -- that table
already exists and means something different (the starship element's own
begin/end within its host genome, not a gene call within a ship).

Revision ID: n7o8p9q0r1s2
Revises: m6n7o8p9q0r1
Create Date: 2026-08-05

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "n7o8p9q0r1s2"
down_revision: Union[str, None] = "m6n7o8p9q0r1"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "annotations",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("sequence", sa.Text(), nullable=False, unique=True),
        sa.Column(
            "annotation", sa.String(255), nullable=True, server_default="No Annotation"
        ),
        sa.Column("public_notes", sa.Text(), nullable=True, server_default=""),
        sa.Column("private_notes", sa.Text(), nullable=True, server_default=""),
        sa.Column("flag", sa.Integer(), nullable=False, server_default="7"),
        sa.Column("assigned_to", sa.String(255), nullable=True),
        sa.Column(
            "created_at",
            sa.DateTime(),
            nullable=False,
            server_default=sa.func.current_timestamp(),
        ),
        sa.Column("updated_at", sa.DateTime(), nullable=True),
    )
    op.create_index("idx_annotations_flag", "annotations", ["flag"])
    op.create_index("idx_annotations_assigned_to", "annotations", ["assigned_to"])

    op.create_table(
        "gene_features",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("joined_ship_id", sa.Integer(), nullable=False),
        sa.Column("annotation_id", sa.Integer(), nullable=True),
        sa.Column("start", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("stop", sa.Integer(), nullable=False, server_default="0"),
        sa.Column("type", sa.String(50), nullable=True),
        sa.Column("strand", sa.String(1), nullable=True),
        sa.ForeignKeyConstraint(
            ["joined_ship_id"], ["joined_ships.id"], ondelete="CASCADE"
        ),
        sa.ForeignKeyConstraint(
            ["annotation_id"], ["annotations.id"], ondelete="SET NULL"
        ),
    )
    op.create_index(
        "idx_gene_features_joined_ship_id", "gene_features", ["joined_ship_id"]
    )
    op.create_index(
        "idx_gene_features_annotation_id", "gene_features", ["annotation_id"]
    )

    op.create_table(
        "annotation_history",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("annotation_id", sa.Integer(), nullable=False),
        sa.Column("changed_by", sa.String(255), nullable=True),
        sa.Column(
            "changed_at",
            sa.DateTime(),
            nullable=False,
            server_default=sa.func.current_timestamp(),
        ),
        sa.Column("old_flag", sa.Integer(), nullable=True),
        sa.Column("new_flag", sa.Integer(), nullable=True),
        sa.Column("old_annotation", sa.Text(), nullable=True),
        sa.Column("new_annotation", sa.Text(), nullable=True),
        sa.Column("old_public_notes", sa.Text(), nullable=True),
        sa.Column("new_public_notes", sa.Text(), nullable=True),
        sa.Column("old_private_notes", sa.Text(), nullable=True),
        sa.Column("new_private_notes", sa.Text(), nullable=True),
        sa.ForeignKeyConstraint(
            ["annotation_id"], ["annotations.id"], ondelete="CASCADE"
        ),
    )
    op.create_index(
        "idx_annotation_history_annotation_id", "annotation_history", ["annotation_id"]
    )

    for table_name, db_comment in [
        ("blastp_results", "swissprot, protein, nr"),
        ("rpsblast_results", "cdd"),
        ("hhsearch_results", "pdb"),
        ("interpro_results", "interpro"),
    ]:
        op.create_table(
            table_name,
            sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
            sa.Column("annotation_id", sa.Integer(), nullable=False),
            sa.Column("database", sa.String(50), nullable=False),
            sa.Column("result", sa.Text(), nullable=True),
            sa.Column("run_date", sa.DateTime(), nullable=True),
            sa.Column("status", sa.Integer(), nullable=False, server_default="0"),
            sa.ForeignKeyConstraint(
                ["annotation_id"], ["annotations.id"], ondelete="CASCADE"
            ),
            sa.UniqueConstraint(
                "annotation_id", "database", name=f"uq_{table_name}_annotation_database"
            ),
        )
        op.create_index(
            f"idx_{table_name}_annotation_id", table_name, ["annotation_id"]
        )


def downgrade() -> None:
    for table_name in [
        "interpro_results",
        "hhsearch_results",
        "rpsblast_results",
        "blastp_results",
    ]:
        op.drop_index(f"idx_{table_name}_annotation_id", table_name)
        op.drop_table(table_name)

    op.drop_index("idx_annotation_history_annotation_id", "annotation_history")
    op.drop_table("annotation_history")

    op.drop_index("idx_gene_features_annotation_id", "gene_features")
    op.drop_index("idx_gene_features_joined_ship_id", "gene_features")
    op.drop_table("gene_features")

    op.drop_index("idx_annotations_assigned_to", "annotations")
    op.drop_index("idx_annotations_flag", "annotations")
    op.drop_table("annotations")

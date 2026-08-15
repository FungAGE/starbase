"""add starfish pipeline tables (starfish_runs, starfish_run_genomes, starfish_elements)

Phase 4 slice 4a: ported from MAS4starships' starship.StarfishRun /
StarfishRunGenome / StarfishElement -- schema + run-management only, no
pipeline execution yet (that's 4b). No data migration, starting fresh.

Renamed MAS4starships' "pid" field to pid_threshold (percent-identity
threshold, a nextflow --pid CLI flag) to avoid confusion with process_pid
(the actual OS process id of the nextflow subprocess, needed to fix
MAS4starships' unreliable cancel -- see report).

Revision ID: o8p9q0r1s2t3
Revises: n7o8p9q0r1s2
Create Date: 2026-08-05

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


revision: str = "o8p9q0r1s2t3"
down_revision: Union[str, None] = "n7o8p9q0r1s2"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "starfish_runs",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("run_name", sa.String(255), nullable=False, unique=True),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("status", sa.String(20), nullable=False, server_default="pending"),
        sa.Column(
            "created_at",
            sa.DateTime(),
            nullable=False,
            server_default=sa.func.current_timestamp(),
        ),
        sa.Column("started_at", sa.DateTime(), nullable=True),
        sa.Column("completed_at", sa.DateTime(), nullable=True),
        sa.Column("created_by", sa.String(255), nullable=True),
        sa.Column("model", sa.String(50), nullable=False, server_default="tyr"),
        sa.Column("threads", sa.Integer(), nullable=False, server_default="20"),
        sa.Column("missing", sa.Integer(), nullable=False, server_default="1"),
        sa.Column("maxcopy", sa.Integer(), nullable=False, server_default="5"),
        sa.Column("pid_threshold", sa.Integer(), nullable=False, server_default="90"),
        sa.Column("hsp", sa.Integer(), nullable=False, server_default="1000"),
        sa.Column("flank", sa.Integer(), nullable=False, server_default="6"),
        sa.Column(
            "neighbourhood", sa.Integer(), nullable=False, server_default="10000"
        ),
        sa.Column("samplesheet_path", sa.String(500), nullable=True),
        sa.Column("output_dir", sa.String(500), nullable=True),
        sa.Column("log_file", sa.String(500), nullable=True),
        sa.Column("celery_task_id", sa.String(255), nullable=True),
        sa.Column("process_pid", sa.Integer(), nullable=True),
        sa.Column("error_message", sa.Text(), nullable=True),
        sa.Column("num_genomes", sa.Integer(), nullable=True),
        sa.Column("num_elements_found", sa.Integer(), nullable=True),
    )
    op.create_index("idx_starfish_runs_status", "starfish_runs", ["status"])

    op.create_table(
        "starfish_run_genomes",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("run_id", sa.Integer(), nullable=False),
        sa.Column("genome_id", sa.String(255), nullable=False),
        sa.Column("tax_id", sa.String(50), nullable=True),
        sa.Column("fna_path", sa.String(500), nullable=False),
        sa.Column("gff3_path", sa.String(500), nullable=False),
        sa.Column("emapper_path", sa.String(500), nullable=True),
        sa.Column("cds_path", sa.String(500), nullable=True),
        sa.Column("faa_path", sa.String(500), nullable=True),
        sa.Column("num_elements", sa.Integer(), nullable=True),
        sa.Column("status", sa.String(20), nullable=False, server_default="pending"),
        sa.Column("error_message", sa.Text(), nullable=True),
        sa.ForeignKeyConstraint(["run_id"], ["starfish_runs.id"], ondelete="CASCADE"),
        sa.UniqueConstraint("run_id", "genome_id"),
    )
    op.create_index(
        "idx_starfish_run_genomes_run_id", "starfish_run_genomes", ["run_id"]
    )

    op.create_table(
        "starfish_elements",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("element_id", sa.String(255), nullable=False, unique=True),
        sa.Column("run_id", sa.Integer(), nullable=False),
        sa.Column("genome_id", sa.Integer(), nullable=False),
        sa.Column("contig_id", sa.String(255), nullable=True),
        sa.Column("start", sa.Integer(), nullable=False),
        sa.Column("end", sa.Integer(), nullable=False),
        sa.Column("strand", sa.String(1), nullable=True),
        sa.Column("sequence", sa.Text(), nullable=True),
        sa.Column("family", sa.String(255), nullable=True),
        sa.Column("navis", sa.String(255), nullable=True),
        sa.Column("haplotype", sa.String(255), nullable=True),
        sa.Column("quality_score", sa.Integer(), nullable=True),
        sa.Column("confidence", sa.String(50), nullable=True),
        sa.Column(
            "created_at",
            sa.DateTime(),
            nullable=False,
            server_default=sa.func.current_timestamp(),
        ),
        sa.Column("notes", sa.Text(), nullable=True),
        sa.Column("imported_submission_id", sa.Integer(), nullable=True),
        sa.ForeignKeyConstraint(["run_id"], ["starfish_runs.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(
            ["genome_id"], ["starfish_run_genomes.id"], ondelete="CASCADE"
        ),
    )
    op.create_index("idx_starfish_elements_run_id", "starfish_elements", ["run_id"])
    op.create_index(
        "idx_starfish_elements_genome_id", "starfish_elements", ["genome_id"]
    )


def downgrade() -> None:
    op.drop_index("idx_starfish_elements_genome_id", "starfish_elements")
    op.drop_index("idx_starfish_elements_run_id", "starfish_elements")
    op.drop_table("starfish_elements")

    op.drop_index("idx_starfish_run_genomes_run_id", "starfish_run_genomes")
    op.drop_table("starfish_run_genomes")

    op.drop_index("idx_starfish_runs_status", "starfish_runs")
    op.drop_table("starfish_runs")

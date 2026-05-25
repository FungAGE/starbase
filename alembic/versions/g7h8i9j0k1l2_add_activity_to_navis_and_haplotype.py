"""add activity column to navis_names and haplotype_names

Revision ID: g7h8i9j0k1l2
Revises: 55c7e854319a
Create Date: 2026-03-06

Adds activity column to navis_names and haplotype_names tables.
When activity = 0, the navis/haplotype info is not shown in wiki modals.
Default is 1 (active) for backward compatibility.
"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa

# revision identifiers, used by Alembic.
revision: str = "g7h8i9j0k1l2"
down_revision: Union[str, None] = "55c7e854319a"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # activity: 0 = inactive (don't show in modals), 1 or NULL = active (show)
    # NULL for existing rows = backward compatible (treat as active)
    op.add_column(
        "navis_names",
        sa.Column("activity", sa.Integer(), nullable=True),
    )
    op.add_column(
        "haplotype_names",
        sa.Column("activity", sa.Integer(), nullable=True),
    )


def downgrade() -> None:
    op.drop_column("haplotype_names", "activity")
    op.drop_column("navis_names", "activity")

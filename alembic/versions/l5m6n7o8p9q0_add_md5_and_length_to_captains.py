"""Add md5 and sequence_length to captains, matching ships table conventions

Revision ID: l5m6n7o8p9q0
Revises: k4l5m6n7o8p9
Create Date: 2026-07-28 00:00:00.000000

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = "l5m6n7o8p9q0"
down_revision: Union[str, None] = "k4l5m6n7o8p9"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Add md5 and sequence_length columns to captains (mirrors ships), plus an
    index on md5 to make sequence-duplicate lookups cheap. Data backfill for
    existing rows is handled separately by a cleanup script, not here."""

    op.add_column("captains", sa.Column("md5", sa.Text(), nullable=True))
    op.add_column("captains", sa.Column("sequence_length", sa.Integer(), nullable=True))
    op.create_index("idx_captains_md5", "captains", ["md5"])


def downgrade() -> None:
    op.drop_index("idx_captains_md5", table_name="captains")
    op.drop_column("captains", "sequence_length")
    op.drop_column("captains", "md5")

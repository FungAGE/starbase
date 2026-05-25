"""add_rev_comp_md5_to_ships

Revision ID: 23eff9174136
Revises: d40154e33df5
Create Date: 2025-08-29 16:24:06.144717

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = "23eff9174136"
down_revision: Union[str, None] = "d40154e33df5"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # Add rev_comp_md5 column to ships table
    op.add_column("ships", sa.Column("rev_comp_md5", sa.String(), nullable=True))


def downgrade() -> None:
    # Remove rev_comp_md5 column from ships table
    op.drop_column("ships", "rev_comp_md5")

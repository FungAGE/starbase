"""merge_heads_for_semantic_versioning

Revision ID: 55c7e854319a
Revises: e789abc12345, a6df69041421
Create Date: 2025-11-18 14:06:51.366427

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = '55c7e854319a'
down_revision: Union[str, None] = ('e789abc12345', 'a6df69041421')
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass

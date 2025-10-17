"""Test migration creation

Revision ID: d40154e33df5
Revises: b29b0a6f346b
Create Date: 2025-08-29 12:03:11.366601

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = 'd40154e33df5'
down_revision: Union[str, None] = 'b29b0a6f346b'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass

"""Enforce sequence_length and elementLength invariants

Revision ID: f8e9a1b2c3d4
Revises: 55c7e854319a
Create Date: 2025-02-18

Enforces:
- ships.sequence_length = length of ships.sequence (character length)
- starship_features.elementLength = elementEnd - elementBegin + 1
"""

from typing import Sequence, Union

from alembic import op


revision: str = "f8e9a1b2c3d4"
down_revision: Union[str, None] = "55c7e854319a"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # ships.sequence_length must equal the length of the string in ships.sequence
    op.execute(
        """
        UPDATE ships
        SET sequence_length = LENGTH(sequence)
        WHERE sequence IS NOT NULL
        """
    )
    # starship_features.elementLength must equal ABS(elementEnd - elementBegin) + 1
    # ABS() handles negative-strand features where elementBegin > elementEnd
    op.execute(
        """
        UPDATE starship_features
        SET elementLength = (ABS(elementEnd - elementBegin) + 1)
        WHERE elementBegin IS NOT NULL AND elementEnd IS NOT NULL
        """
    )


def downgrade() -> None:
    # No reversible way to restore previous values
    pass

"""add_ship_quality_tags_table

Revision ID: a6df69041421
Revises: 23eff9174136
Create Date: 2025-10-08 17:42:03.203152

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = "a6df69041421"
down_revision: Union[str, None] = "23eff9174136"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # Create ship_quality_tags table
    op.create_table(
        "ship_quality_tags",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("joined_ship_id", sa.Integer(), nullable=False),
        sa.Column("tag_type", sa.String(50), nullable=False),
        sa.Column("tag_value", sa.String(100), nullable=True),
        sa.Column(
            "created_at",
            sa.DateTime(),
            nullable=False,
            server_default=sa.func.current_timestamp(),
        ),
        sa.Column("created_by", sa.String(50), nullable=True, default="auto"),
        sa.ForeignKeyConstraint(
            ["joined_ship_id"], ["joined_ships.id"], ondelete="CASCADE"
        ),
        sa.UniqueConstraint(
            "joined_ship_id", "tag_type", name="uq_ship_quality_tags_ship_tag"
        ),
    )

    # Create indexes for efficient querying
    op.create_index(
        "idx_ship_quality_tags_joined_ship_id", "ship_quality_tags", ["joined_ship_id"]
    )
    op.create_index("idx_ship_quality_tags_tag_type", "ship_quality_tags", ["tag_type"])


def downgrade() -> None:
    # Drop indexes first
    op.drop_index("idx_ship_quality_tags_tag_type", "ship_quality_tags")
    op.drop_index("idx_ship_quality_tags_joined_ship_id", "ship_quality_tags")

    # Drop table
    op.drop_table("ship_quality_tags")

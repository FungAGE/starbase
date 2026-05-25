"""Add accession_id to joined_ships for central table design

Revision ID: e789abc12345
Revises: d40154e33df5
Create Date: 2025-09-23 14:30:00.000000

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = "e789abc12345"
down_revision: Union[str, None] = "d40154e33df5"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Add accession_id column to joined_ships and migrate data"""

    # Add the new accession_id column with foreign key constraint
    op.add_column(
        "joined_ships",
        sa.Column(
            "accession_id", sa.Integer(), sa.ForeignKey("accessions.id"), nullable=True
        ),
    )

    # Note: Foreign key constraint is included in the column definition above
    # This is the proper way to add a column with FK constraint in Alembic

    # Add index for performance
    op.create_index("idx_joined_ships_accession_id", "joined_ships", ["accession_id"])

    # Make starshipID not null (duplicates allowed for cases where multiple entries reference same ship)
    op.alter_column(
        "joined_ships", "starshipID", existing_type=sa.String(), nullable=False
    )
    # Note: No unique constraint on starshipID to allow duplicates when needed

    # Migrate existing ship_id references to accession_id
    # This requires a data migration that should be done carefully
    # The migration script migrate_to_central_joined_ships.py should be run separately


def downgrade() -> None:
    """Remove accession_id column and related constraints"""

    # Remove index and foreign key constraint
    op.drop_index("idx_joined_ships_accession_id", table_name="joined_ships")
    # Note: FK constraint is dropped automatically when column is dropped

    # Remove the column
    op.drop_column("joined_ships", "accession_id")

    # Restore starshipID to nullable
    op.alter_column(
        "joined_ships", "starshipID", existing_type=sa.String(), nullable=True
    )

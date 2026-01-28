"""add_ship_accessions_table

Revision ID: 8adfce5239bb
Revises: a6df69041421
Create Date: 2026-01-15 16:57:16.000000

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = '8adfce5239bb'
down_revision: Union[str, None] = 'a6df69041421'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # Create ship_accessions table
    op.create_table(
        'ship_accessions',
        sa.Column('id', sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column('ship_accession_tag', sa.String(), nullable=False),
        sa.Column('version_tag', sa.Integer(), nullable=False, server_default='1'),
        sa.Column('ship_id', sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(['ship_id'], ['ships.id'], ondelete='CASCADE'),
        sa.UniqueConstraint('ship_accession_tag', name='uq_ship_accessions_tag'),
        sa.UniqueConstraint('ship_id', name='uq_ship_accessions_ship_id')
    )
    
    # Create indexes for efficient querying
    op.create_index('idx_ship_accessions_ship_id', 'ship_accessions', ['ship_id'])
    op.create_index('idx_ship_accessions_tag', 'ship_accessions', ['ship_accession_tag'])


def downgrade() -> None:
    # Drop indexes first
    op.drop_index('idx_ship_accessions_tag', 'ship_accessions')
    op.drop_index('idx_ship_accessions_ship_id', 'ship_accessions')
    
    # Drop table
    op.drop_table('ship_accessions')

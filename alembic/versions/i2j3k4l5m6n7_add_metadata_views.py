"""add ships_with_metadata, captains_with_metadata, ship_table_view

Merges the 4 divergent heads and adds read-only SQL views that mirror the
JOINs previously duplicated as inline CTEs/queries in sql_manager.py:
- ships_with_metadata: mirrors fetch_ships()'s CTE (full join incl. captains/features)
- captains_with_metadata: mirrors fetch_captains()'s CTE
- ship_table_view: mirrors fetch_ship_table()'s narrower join (no features/genomes/
  navis/haplotype/captains, to avoid introducing row fan-out from starship_features)

Revision ID: i2j3k4l5m6n7
Revises: 8adfce5239bb, f1a2b3c4d5e6, f8e9a1b2c3d4, h1i2j3k4l5m6
Create Date: 2026-07-23

"""

from typing import Sequence, Union

from alembic import op

revision: str = "i2j3k4l5m6n7"
down_revision: Union[str, Sequence[str], None] = (
    "8adfce5239bb",
    "f1a2b3c4d5e6",
    "f8e9a1b2c3d4",
    "h1i2j3k4l5m6",
)
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.execute("DROP VIEW IF EXISTS ships_with_metadata")
    op.execute("""
        CREATE VIEW ships_with_metadata AS
        SELECT DISTINCT
            j.id as joined_ship_id,
            j.ship_id,
            j.source,
            j.evidence,
            sa.ship_accession_tag,
            sa.ship_version_tag,
            sa.ship_accession_display,
            j.curated_status,
            j.starshipID,
            sf.elementBegin, sf.elementEnd, sf.elementLength, sf.contigID, sf.upDR, sf.downDR,
            t.name, t.family, t."order", t.taxID, t.strain,
            f.familyName, f.type_element_reference,
            n.navis_name, n.activity as navis_activity,
            h.haplotype_name, h.activity as haplotype_activity,
            g.ome, g.version as genome_version, g.genomeSource, g.citation, g.assembly_accession,
            c.captainID,
            a.accession_tag, a.version_tag, a.accession_display,
            s.sequence, s.md5, s.rev_comp_md5, s.type_ship
        FROM joined_ships j
        LEFT JOIN ships s ON s.id = j.ship_id
        LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        LEFT JOIN starship_features sf ON j.ship_id = sf.ship_id
        LEFT JOIN captains c ON j.captain_id = c.id
        LEFT JOIN accessions a ON j.accession_id = a.id
        WHERE j.ship_id IS NOT NULL
    """)

    op.execute("DROP VIEW IF EXISTS captains_with_metadata")
    op.execute("""
        CREATE VIEW captains_with_metadata AS
        SELECT DISTINCT
            sa.ship_accession_tag,
            sa.ship_version_tag,
            sa.ship_accession_display,
            j.curated_status,
            j.starshipID,
            c.id as captain_row_id,
            c.captainID,
            c.sequence as captain_sequence,
            n.navis_name,
            h.haplotype_name,
            a.accession_tag,
            a.version_tag,
            a.accession_display
        FROM joined_ships j
        LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        LEFT JOIN captains c ON j.captain_id = c.id
        LEFT JOIN starship_features sf ON j.ship_id = sf.ship_id
        LEFT JOIN accessions a ON j.accession_id = a.id
        WHERE j.ship_id IS NOT NULL
    """)

    op.execute("DROP VIEW IF EXISTS ship_table_view")
    op.execute("""
        CREATE VIEW ship_table_view AS
        SELECT DISTINCT
            js.ship_id,
            js.source,
            js.curated_status,
            sa.ship_accession_tag,
            sa.ship_version_tag,
            sa.ship_accession_display,
            f.familyName,
            t.name,
            a.accession_tag, a.version_tag, a.accession_display
        FROM joined_ships js
        LEFT JOIN ship_accessions sa ON sa.ship_id = js.ship_id
        LEFT JOIN taxonomy t ON js.tax_id = t.id
        LEFT JOIN family_names f ON js.ship_family_id = f.id
        LEFT JOIN accessions a ON js.accession_id = a.id
        WHERE js.ship_id IS NOT NULL
    """)


def downgrade() -> None:
    op.execute("DROP VIEW IF EXISTS ship_table_view")
    op.execute("DROP VIEW IF EXISTS captains_with_metadata")
    op.execute("DROP VIEW IF EXISTS ships_with_metadata")

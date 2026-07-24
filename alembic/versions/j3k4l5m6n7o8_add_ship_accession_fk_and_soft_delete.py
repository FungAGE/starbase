"""add ship_accession_id fk and is_deleted soft-delete flag to joined_ships

joined_ships already links directly to accessions (SSA, group-level) via
accession_id. The SSB (ship-level) link only existed indirectly, via
joined_ships.ship_id -> ships.id -> ship_accessions.ship_id. This adds a
direct ship_accession_id FK on joined_ships mirroring accession_id, and
backfills it from the existing ship_id linkage.

Also adds is_deleted (soft delete flag) to joined_ships, and updates the
3 metadata views (see i2j3k4l5m6n7) to exclude soft-deleted rows so every
existing caller gets safe behavior with no code changes.

Revision ID: j3k4l5m6n7o8
Revises: i2j3k4l5m6n7
Create Date: 2026-07-24

"""

from typing import Sequence, Union

from alembic import op

revision: str = "j3k4l5m6n7o8"
down_revision: Union[str, None] = "i2j3k4l5m6n7"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.execute(
        "ALTER TABLE joined_ships ADD COLUMN ship_accession_id INTEGER REFERENCES ship_accessions(id)"
    )
    op.execute("""
        UPDATE joined_ships
        SET ship_accession_id = (
            SELECT sa.id FROM ship_accessions sa WHERE sa.ship_id = joined_ships.ship_id
        )
        WHERE ship_id IS NOT NULL
    """)

    op.execute(
        "ALTER TABLE joined_ships ADD COLUMN is_deleted BOOLEAN NOT NULL DEFAULT 0"
    )

    _create_views(deleted_filter=True)


def downgrade() -> None:
    _create_views(deleted_filter=False)

    op.execute("ALTER TABLE joined_ships DROP COLUMN is_deleted")
    op.execute("ALTER TABLE joined_ships DROP COLUMN ship_accession_id")


def _create_views(deleted_filter: bool) -> None:
    deleted_clause = "AND j.is_deleted = 0" if deleted_filter else ""

    op.execute("DROP VIEW IF EXISTS ships_with_metadata")
    op.execute(f"""
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
        WHERE j.ship_id IS NOT NULL {deleted_clause}
    """)

    op.execute("DROP VIEW IF EXISTS captains_with_metadata")
    op.execute(f"""
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
        WHERE j.ship_id IS NOT NULL {deleted_clause}
    """)

    op.execute("DROP VIEW IF EXISTS ship_table_view")
    op.execute(f"""
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
        WHERE js.ship_id IS NOT NULL {deleted_clause.replace("j.", "js.")}
    """)

"""add genus/species/class taxonomy columns to ships_with_metadata

wiki.py's taxa search is being fixed to search exactly: name, class, order,
family, genus, species. order/family/name were already present in
ships_with_metadata (and therefore fetch_meta_data, which backs the wiki page's
meta-data store), but genus, species, and class were never selected from
taxonomy -- so search on those 3 levels was always a silent no-op.

Revision ID: k4l5m6n7o8p9
Revises: j3k4l5m6n7o8
Create Date: 2026-07-24

"""

from typing import Sequence, Union

from alembic import op

revision: str = "k4l5m6n7o8p9"
down_revision: Union[str, None] = "j3k4l5m6n7o8"
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
            t.name, t.genus, t.species, t.family, t."order", t.class, t.taxID, t.strain,
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
        WHERE j.ship_id IS NOT NULL AND j.is_deleted = 0
    """)


def downgrade() -> None:
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
        WHERE j.ship_id IS NOT NULL AND j.is_deleted = 0
    """)

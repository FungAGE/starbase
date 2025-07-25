from sqlalchemy import text
from src.config.database import StarbaseSession

from src.config.logging import get_logger

logger = get_logger(__name__)


def create_database_indexes():
    """Create indexes to optimize database queries"""
    session = StarbaseSession()

    indexes = [
        # joined_ships indexes
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_curated ON joined_ships(curated_status)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_id ON joined_ships(ship_id)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_tax_id ON joined_ships(tax_id)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_family_id ON joined_ships(ship_family_id)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_navis_id ON joined_ships(ship_navis_id)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_haplotype_id ON joined_ships(ship_haplotype_id)",
        "CREATE INDEX IF NOT EXISTS idx_joined_ships_genome_id ON joined_ships(genome_id)",
        # accessions indexes
        "CREATE INDEX IF NOT EXISTS idx_accessions_tag ON accessions(accession_tag)",
        "CREATE INDEX IF NOT EXISTS idx_accessions_id ON accessions(id)",
        # taxonomy indexes
        "CREATE INDEX IF NOT EXISTS idx_taxonomy_name ON taxonomy(name)",
        "CREATE INDEX IF NOT EXISTS idx_taxonomy_taxid ON taxonomy(taxID)",
        # family_names indexes
        "CREATE INDEX IF NOT EXISTS idx_family_names_reference ON family_names(type_element_reference)",
        "CREATE INDEX IF NOT EXISTS idx_family_names_family ON family_names(familyName)",
        # ships indexes
        "CREATE INDEX IF NOT EXISTS idx_ships_accession_id ON ships(accession_id)",
        # gff indexes
        "CREATE INDEX IF NOT EXISTS idx_gff_ship_id ON gff(ship_id)",
    ]

    try:
        for index_sql in indexes:
            logger.debug(f"Creating index: {index_sql}")
            session.execute(text(index_sql))
        session.commit()
        logger.debug("Successfully created all database indexes")
    except Exception as e:
        logger.error(f"Error creating indexes: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()

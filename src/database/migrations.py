from sqlalchemy import text
from src.config.database import StarbaseSession

from src.config.logging import get_logger

logger = get_logger(__name__)


def create_database_indexes():
    """Create indexes to optimize database queries"""
    session = StarbaseSession()

    # Check if indexes already exist by querying sqlite_master
    existing_indexes_query = """
    SELECT name FROM sqlite_master 
    WHERE type='index' AND name LIKE 'idx_%'
    """
    existing_indexes = {row[0] for row in session.execute(text(existing_indexes_query))}
    
    logger.info(f"Found {len(existing_indexes)} existing indexes")

    indexes = [
        # joined_ships indexes
        ("idx_joined_ships_curated", "CREATE INDEX IF NOT EXISTS idx_joined_ships_curated ON joined_ships(curated_status)"),
        ("idx_joined_ships_accession_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_accession_id ON joined_ships(accession_id)"),
        ("idx_joined_ships_tax_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_tax_id ON joined_ships(tax_id)"),
        ("idx_joined_ships_ship_family_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_family_id ON joined_ships(ship_family_id)"),
        ("idx_joined_ships_ship_navis_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_navis_id ON joined_ships(ship_navis_id)"),
        ("idx_joined_ships_ship_haplotype_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_ship_haplotype_id ON joined_ships(ship_haplotype_id)"),
        ("idx_joined_ships_genome_id", "CREATE INDEX IF NOT EXISTS idx_joined_ships_genome_id ON joined_ships(genome_id)"),
        # accessions indexes
        ("idx_accessions_tag", "CREATE INDEX IF NOT EXISTS idx_accessions_tag ON accessions(accession_tag)"),
        ("idx_accessions_id", "CREATE INDEX IF NOT EXISTS idx_accessions_id ON accessions(id)"),
        # taxonomy indexes
        ("idx_taxonomy_name", "CREATE INDEX IF NOT EXISTS idx_taxonomy_name ON taxonomy(name)"),
        ("idx_taxonomy_taxid", "CREATE INDEX IF NOT EXISTS idx_taxonomy_taxid ON taxonomy(taxID)"),
        # family_names indexes
        ("idx_family_names_reference", "CREATE INDEX IF NOT EXISTS idx_family_names_reference ON family_names(type_element_reference)"),
        ("idx_family_names_family", "CREATE INDEX IF NOT EXISTS idx_family_names_family ON family_names(familyName)"),
        # ships indexes
        ("idx_ships_accession_id", "CREATE INDEX IF NOT EXISTS idx_ships_accession_id ON ships(accession_id)"),
        # gff indexes
        ("idx_gff_ship_id", "CREATE INDEX IF NOT EXISTS idx_gff_ship_id ON gff(ship_id)"),
    ]

    created_count = 0
    try:
        for index_name, index_sql in indexes:
            if index_name not in existing_indexes:
                logger.debug(f"Creating index: {index_name}")
                session.execute(text(index_sql))
                created_count += 1
            else:
                logger.debug(f"Index {index_name} already exists, skipping")
        
        if created_count > 0:
            session.commit()
            logger.info(f"Successfully created {created_count} new database indexes")
        else:
            logger.info("All database indexes already exist, no new indexes created")
    except Exception as e:
        logger.error(f"Error creating indexes: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()

import sys
import os

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
)

from sqlalchemy import create_engine
import pandas as pd
from src.database.models.schema import Base  # Import all models
from src.config.settings import DB_PATHS
from datetime import datetime
import logging
from src.utils.seq_utils import generate_checksum
from inspect import signature

logger = logging.getLogger(__name__)


def setup_engines():
    """Create connections to both old and new databases"""
    old_engine = create_engine(f"sqlite:///{DB_PATHS['starbase']}")
    new_engine = create_engine(f"sqlite:///{DB_PATHS['starbase']}_new")
    return old_engine, new_engine


def create_new_schema(engine):
    """Create all tables with new schema"""
    Base.metadata.create_all(engine)


def add_ship_sequences(old_engine, new_engine, num_rows=10):
    """Create ship sequence table with foreign key to accessions"""
    # First, get the accession IDs from the new database
    query_new = "SELECT id, accession_tag FROM accessions"
    accessions_df = pd.read_sql(query_new, new_engine)

    # Then get the sequences from the old database
    query_old = f"SELECT a.accession_tag, s.sequence FROM accessions a JOIN ships s ON a.id = s.accession LIMIT {num_rows}"
    sequences_df = pd.read_sql(query_old, old_engine)

    # Merge the dataframes to get the new accession IDs
    df = pd.merge(sequences_df, accessions_df, on="accession_tag", how="inner")

    # Add computed columns
    df["md5"] = df["sequence"].apply(generate_checksum)
    df["sequence_length"] = df["sequence"].apply(len)

    # Select only the columns we want to insert
    df = df[["id", "sequence", "md5", "sequence_length"]].rename(
        columns={"id": "accession_id"}
    )

    return df


def create_accession_table(old_engine, num_rows=10):
    """Create accession table"""
    query = f"""
    SELECT 
        ship_name,
        accession_tag,
        accession_tag as accession
    FROM accessions 
    LIMIT {num_rows}
    """
    df = pd.read_sql(query, old_engine)

    # Ensure required fields are not null
    df["ship_name"] = df["ship_name"].fillna("Unknown")

    # Add timestamps and soft delete flag
    df["created_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["is_deleted"] = False

    return df


def create_starship_features_table(old_engine, num_rows=10):
    """Create starship features table"""
    # get first 10 rows from existing starship features table
    query = f"SELECT * FROM starship_features LIMIT {num_rows}"
    df = pd.read_sql(query, old_engine)

    return df


# TODO: get accession tags another way, instead of from input
def create_captain_table(new_engine):
    """Create captain table"""

    # get accession tags from new database
    query = "SELECT accession_tag FROM accessions"
    accession_tags_df = pd.read_sql(query, new_engine)
    accession_tags = accession_tags_df["accession_tag"].tolist()

    # use existing method to get captain table
    from src.database.sql_manager import fetch_captains

    df = fetch_captains(
        accession_tags=accession_tags,
        curated=False,
        dereplicate=False,
        with_sequence=True,
    )

    # add checksum to dataframe
    df["md5"] = df["sequence"].apply(generate_checksum)

    # count the length of the sequence
    df["sequence_length"] = df["sequence"].apply(len)

    # write to new captain table
    df.to_sql("captain", new_engine, if_exists="append", index=False)


def create_classification_table(old_engine, new_engine, num_rows=10):
    """Create classification table"""
    # get first 10 rows from existing classification table
    query = f"SELECT * FROM classification LIMIT {num_rows}"
    df = pd.read_sql(query, old_engine)

    return df


def migrate_table(
    old_engine=None, new_engine=None, table_name=None, transform_func=None
):
    """Migrate a single table with optional data transformation"""
    try:
        # Handle transformation
        if transform_func:
            # Check if function expects one or two arguments
            params = signature(transform_func).parameters
            if len(params) == 2:  # Function expects both engines
                df = transform_func(old_engine, new_engine)
            else:  # Function expects only old_engine
                df = transform_func(old_engine)

        # Write to new table
        if new_engine is None:
            new_engine = old_engine  # Use old_engine if new_engine is None

        with new_engine.begin() as conn:
            df.to_sql(table_name, conn, if_exists="append", index=False)

        return len(df)
    except Exception as e:
        logger.error(f"Error migrating table {table_name}: {e}")
        raise


def transform_starship_features(df):
    """Example transformation function for starship_features table"""
    # Convert string numbers to integers
    for col in ["elementBegin", "elementEnd", "elementLength"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Validate and clean strand values
    df["strand"] = df["strand"].apply(lambda x: x if x in ["+", "-"] else "+")

    # Set default boundary type
    df["boundaryType"] = df["boundaryType"].fillna("unknown")

    # Add timestamps
    df["created_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    return df


def create_taxonomy_table(old_engine):
    return None


def create_family_names_table(old_engine):
    return None


def create_genomes_table(old_engine):
    return None


def create_gff_table(old_engine):
    return None


def create_haplotype_names_table(old_engine):
    return None


def create_joined_ships_table(old_engine):
    return None


def create_navis_names_table(old_engine):
    return None


def create_papers_table(old_engine):
    return None


def check_accessions(old_engine, new_engine):
    """Check if accessions have been assigned correctly"""
    query = """
    SELECT a.id, a.accession_tag, s.sequence_length, s.md5
    FROM accession a
    INNER JOIN ship_sequence s ON a.id = s.ship_id
    """
    df = pd.read_sql(query, new_engine)
    return df


def check_classification(old_engine, new_engine):
    """Check if classification has been assigned correctly"""
    query = """
    SELECT * FROM classification
    """
    df = pd.read_sql(query, new_engine)
    return df


def validate_relationships(new_engine):
    """Validate that relationships between tables are maintained"""
    # Check ships-accessions relationship
    query = """
    SELECT 
        a.accession_tag,
        s.sequence_length,
        s.md5
    FROM accessions a
    LEFT JOIN ships s ON a.id = s.accession_id
    WHERE s.accession_id IS NULL
    """
    df = pd.read_sql(query, new_engine)
    if not df.empty:
        logger.warning(
            f"Found {len(df)} accessions without corresponding ship sequences"
        )
        logger.warning(df)


def validate_foreign_keys(new_engine, table_name):
    """Validate foreign key relationships"""
    # Check if all ships have valid accession_ids
    query = f"SELECT s.id, s.accession_id, a.id as valid_accession FROM ships s LEFT JOIN {table_name} a ON s.accession_id = a.id WHERE a.id IS NULL"
    invalid_relations = pd.read_sql(query, new_engine)

    if not invalid_relations.empty:
        logger.error(
            f"Found {len(invalid_relations)} invalid accession_ids in {table_name}"
        )
        logger.error(invalid_relations)
        raise ValueError(f"Foreign key validation failed for table: {table_name}")


def main():
    """
    steps to create database from existing database
    1. create new database with updated schema
    2. create ship sequence table
    3. create accession table
    4. create captain sequence table
    5. create starship features table
    6. create classification results table
    """

    # remove new database, if it exists
    if os.path.exists(f"{DB_PATHS['starbase']}_new"):
        os.remove(f"{DB_PATHS['starbase']}_new")

    logger.info("Setting up databases...")
    old_engine, new_engine = setup_engines()

    logger.info("Creating new database schema...")
    create_new_schema(new_engine)

    # Migrate tables in order of dependencies
    # First migrate accessions (parent table)
    logger.info("Migrating accessions table...")
    rows_migrated = migrate_table(
        old_engine=old_engine,
        new_engine=None,
        table_name="accessions",
        transform_func=create_accession_table,
    )
    logger.info(f"Migrated {rows_migrated} rows from accessions")

    # Then migrate dependent tables
    dependent_tables = {
        "ships": add_ship_sequences,
        "captains": create_captain_table,
        # "starship_features": transform_starship_features,
        # "classification": create_classification_table,
        # "taxonomy": create_taxonomy_table,
        # "family_names": create_family_names_table,
        # "genomes": create_genomes_table,
        # "gff": create_gff_table,
        # "haplotype_names": create_haplotype_names_table,
        # "joined_ships": create_joined_ships_table,
        # "haplotype_names": create_haplotype_names_table,
        # "navis_names": create_navis_names_table,
        # "papers": create_papers_table,
    }

    for table_name, transform_func in dependent_tables.items():
        logger.info(f"Migrating table: {table_name}")
        rows_migrated = migrate_table(
            old_engine=old_engine,
            new_engine=new_engine,
            table_name=table_name,
            transform_func=transform_func,
        )
        logger.info(f"Migrated {rows_migrated} rows from {table_name}")

        # Validate foreign keys
        logger.info(f"Validating foreign key relationships for {table_name}...")
        validate_foreign_keys(new_engine, table_name)

    logger.info("Migration completed successfully")


if __name__ == "__main__":
    main()

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
        accession_tag
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


def fetch_captains_for_migration(
    old_engine,
    accession_tags=None,
    curated=False,
    dereplicate=True,
    with_sequence=False,
):
    """
    Migration-specific version of fetch_captains that works with the old database schema.
    This function directly queries the old database without using SQLAlchemy session.

    Args:
        old_engine: Database engine for the old database
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all captains.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing captain data
    """

    # First, let's try a simpler query to see what tables and columns exist
    if with_sequence:
        query = """
        SELECT DISTINCT 
            a.id, 
            a.accession_tag,
            j.curated_status,
            j.starshipID,
            j.captainID,
            j.captain_id,
            c.sequence
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN captains c ON j.captain_id = c.id
        WHERE 1=1
        """
    else:
        query = """
        SELECT DISTINCT 
            a.id, 
            a.accession_tag,
            j.curated_status,
            j.starshipID,
            j.captainID,
            j.captain_id
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        WHERE 1=1
        """

    if accession_tags:
        query += " AND a.accession_tag IN ({})".format(
            ",".join(f"'{tag}'" for tag in accession_tags)
        )
    if curated:
        query += " AND j.curated_status = 'curated'"

    if with_sequence:
        query += " AND (c.sequence IS NOT NULL)"

    try:
        df = pd.read_sql_query(query, old_engine)

        if dereplicate:
            df = df.drop_duplicates(subset="accession_tag")

        if df.empty:
            logger.warning("Fetched captains DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching captains data for migration: {str(e)}")
        # If the above query fails, try a simpler fallback
        try:
            fallback_query = """
            SELECT DISTINCT 
                a.id, 
                a.accession_tag,
                j.curated_status,
                j.starshipID,
                j.captainID,
                j.captain_id
            FROM joined_ships j
            INNER JOIN accessions a ON j.ship_id = a.id
            WHERE 1=1
            """

            if accession_tags:
                fallback_query += " AND a.accession_tag IN ({})".format(
                    ",".join(f"'{tag}'" for tag in accession_tags)
                )
            if curated:
                fallback_query += " AND j.curated_status = 'curated'"

            df = pd.read_sql_query(fallback_query, old_engine)

            # Add empty sequence column if requested
            if with_sequence:
                df["sequence"] = None

            if dereplicate:
                df = df.drop_duplicates(subset="accession_tag")

            return df
        except Exception as fallback_error:
            logger.error(f"Fallback query also failed: {str(fallback_error)}")
            raise e


# TODO: get accession tags another way, instead of from input
def create_captain_table(old_engine, new_engine):
    """Create captain table"""

    # get accession tags from new database
    query = "SELECT id, accession_tag FROM accessions"
    accessions_df = pd.read_sql(query, new_engine)
    accession_tags = accessions_df["accession_tag"].tolist()

    df = fetch_captains_for_migration(
        old_engine=old_engine,
        accession_tags=accession_tags,
        curated=False,
        dereplicate=False,
        with_sequence=True,
    )

    # Map accession_tag to ship_id (foreign key)
    # Rename the id column in accessions_df to avoid conflicts
    accessions_df = accessions_df.rename(columns={"id": "accession_id"})
    df = pd.merge(df, accessions_df, on="accession_tag", how="inner")

    # Select and rename columns to match the schema
    # Schema columns: id, captain_name, sequence, accession_id, reviewed, evidence
    df_mapped = pd.DataFrame(
        {
            "captain_name": df.get(
                "captainID", df.get("captain_id", "unknown")
            ),  # Use available captain ID
            "sequence": df["sequence"],
            "accession_id": df["accession_id"],  # Now this should work
            "reviewed": df.get(
                "curated_status", "not curated"
            ),  # Map curated_status to reviewed
            "evidence": df.get(
                "evidence", "unknown"
            ),  # Use evidence if available, otherwise default
        }
    )

    return df_mapped


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
            # Check function signature to determine how to call it
            params = signature(transform_func).parameters
            param_names = list(params.keys())

            # Check if function expects new_engine as a parameter
            if "new_engine" in param_names:
                df = transform_func(old_engine, new_engine)
            else:
                # Function only expects old_engine (and possibly num_rows)
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
        new_engine=new_engine,
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

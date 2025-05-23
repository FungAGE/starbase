from sqlalchemy import create_engine
import pandas as pd
from src.database.models.schema import Base  # Import all models
from src.config.settings import DB_PATHS
from datetime import datetime
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def setup_engines():
    """Create connections to both old and new databases"""
    old_engine = create_engine(f"sqlite:///{DB_PATHS['starbase']}")
    new_engine = create_engine(f"sqlite:///{DB_PATHS['starbase']}_new")
    return old_engine, new_engine


def create_new_schema(engine):
    """Create all tables with new schema"""
    Base.metadata.create_all(engine)


def migrate_table(old_engine, new_engine, table_name, transform_func=None):
    """Migrate a single table with optional data transformation"""
    try:
        # Read data from old table
        query = f"SELECT * FROM {table_name}"
        df = pd.read_sql(query, old_engine)

        if df.empty:
            logger.warning(f"No data found in table {table_name}")
            return 0

        if transform_func:
            df = transform_func(df)

        # Write to new table
        with new_engine.begin() as conn:
            df.to_sql(table_name, conn, if_exists="append", index=False)

        return len(df)
    except Exception as e:
        logger.error(f"Error migrating table {table_name}: {e}")
        raise


def transform_accessions(df):
    """Example transformation function for accessions table"""
    # Add new required columns
    df["created_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["is_deleted"] = False

    # Ensure required fields are not null
    df["ship_name"] = df["ship_name"].fillna("Unknown")

    # Truncate strings that are too long
    df["ship_name"] = df["ship_name"].str.slice(0, 255)
    df["accession"] = df["accession"].str.slice(0, 50)

    return df


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


def main():
    old_engine, new_engine = setup_engines()

    # Create new database with updated schema
    logger.info("Creating new database schema...")
    create_new_schema(new_engine)

    # Dictionary of tables and their transformation functions
    migrations = {
        "accessions": transform_accessions,
        "ships": None,  # No transformation needed
        "captains": None,
        "starship_features": transform_starship_features,
        # Add other tables...
    }

    # Migrate each table
    for table_name, transform_func in migrations.items():
        logger.info(f"Migrating table: {table_name}")
        rows_migrated = migrate_table(
            old_engine, new_engine, table_name, transform_func
        )
        logger.info(f"Migrated {rows_migrated} rows from {table_name}")

    logger.info("Migration completed successfully")


if __name__ == "__main__":
    main()

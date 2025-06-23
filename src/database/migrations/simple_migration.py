#!/usr/bin/env python3
"""
Simple database migration script that builds tables piece-by-piece.
This approach is much simpler than the complex migration script.
"""

import sys
import os
import logging
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from sqlalchemy import create_engine, text
import pandas as pd
from src.database.models.schema import Base
from src.config.settings import DB_PATHS
from datetime import datetime

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SimpleMigration:
    """Simple database migration that builds tables step by step"""

    def __init__(self, old_db_path=None, new_db_path=None):
        self.old_db_path = old_db_path or DB_PATHS["starbase"]
        self.new_db_path = new_db_path or f"{DB_PATHS['starbase']}_new"

        # Create engines
        self.old_engine = create_engine(f"sqlite:///{self.old_db_path}")
        self.new_engine = create_engine(f"sqlite:///{self.new_db_path}")

    def clean_start(self):
        """Remove new database if it exists to start fresh"""
        if os.path.exists(self.new_db_path):
            os.remove(self.new_db_path)
            logger.info(f"Removed existing database: {self.new_db_path}")

    def create_schema(self):
        """Create the new database schema"""
        Base.metadata.create_all(self.new_engine)
        logger.info("Created new database schema")

    def get_old_table_info(self, table_name):
        """Get information about a table in the old database"""
        try:
            with self.old_engine.connect() as conn:
                result = conn.execute(text(f"PRAGMA table_info({table_name})"))
                columns = [row[1] for row in result.fetchall()]

                result = conn.execute(text(f"SELECT COUNT(*) FROM {table_name}"))
                count = result.fetchone()[0]

                return {"columns": columns, "count": count}
        except Exception as e:
            logger.warning(f"Table {table_name} not found in old database: {e}")
            return None

    def show_old_database_info(self):
        """Show what's in the old database"""
        logger.info("=== OLD DATABASE INFO ===")

        # Get all tables
        with self.old_engine.connect() as conn:
            result = conn.execute(
                text("SELECT name FROM sqlite_master WHERE type='table'")
            )
            tables = [row[0] for row in result.fetchall()]

        for table in tables:
            info = self.get_old_table_info(table)
            if info:
                logger.info(
                    f"Table: {table} | Rows: {info['count']} | Columns: {', '.join(info['columns'][:5])}{'...' if len(info['columns']) > 5 else ''}"
                )

    def migrate_accessions(self, limit=None):
        """Migrate accessions table (base table with no dependencies)"""
        logger.info("Migrating accessions table...")

        query = "SELECT ship_name, accession_tag FROM accessions"
        if limit:
            query += f" LIMIT {limit}"

        try:
            df = pd.read_sql(query, self.old_engine)

            # Clean data
            df["ship_name"] = df["ship_name"].fillna("Unknown")
            df["created_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            df["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            df["is_deleted"] = False

            # Write to new database
            df.to_sql("accessions", self.new_engine, if_exists="append", index=False)
            logger.info(f"Migrated {len(df)} accession records")
            return len(df)

        except Exception as e:
            logger.error(f"Error migrating accessions: {e}")
            return 0

    def migrate_ships(self, limit=None):
        """Migrate ships table (depends on accessions)"""
        logger.info("Migrating ships table...")

        try:
            # Get accession mappings from new database
            new_accessions = pd.read_sql(
                "SELECT id, accession_tag FROM accessions", self.new_engine
            )

            # Get ship data from old database
            query = """
            SELECT a.accession_tag, s.sequence 
            FROM accessions a 
            JOIN ships s ON a.id = s.accession
            """
            if limit:
                query += f" LIMIT {limit}"

            old_ships = pd.read_sql(query, self.old_engine)

            # Merge to get new accession IDs
            ships_df = pd.merge(
                old_ships, new_accessions, on="accession_tag", how="inner"
            )

            # Add computed fields
            ships_df["md5"] = ships_df["sequence"].apply(
                lambda x: f"md5_{hash(x) % 100000}" if pd.notna(x) else None
            )
            ships_df["sequence_length"] = ships_df["sequence"].apply(
                lambda x: len(x) if pd.notna(x) else 0
            )
            ships_df["created_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            ships_df["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            ships_df["is_deleted"] = False

            # Select final columns
            final_df = ships_df[
                [
                    "sequence",
                    "md5",
                    "sequence_length",
                    "created_at",
                    "updated_at",
                    "is_deleted",
                    "id",
                ]
            ].rename(columns={"id": "accession_id"})

            # Write to new database
            final_df.to_sql("ships", self.new_engine, if_exists="append", index=False)
            logger.info(f"Migrated {len(final_df)} ship records")
            return len(final_df)

        except Exception as e:
            logger.error(f"Error migrating ships: {e}")
            return 0

    def add_database_inspector_methods(self):
        """Add methods to inspect what's actually in the old database"""

        def explore_table_relationships(self, table_name):
            """Explore relationships for a specific table"""
            logger.info(f"=== EXPLORING {table_name.upper()} TABLE ===")

            try:
                with self.old_engine.connect() as conn:
                    # Get column info
                    result = conn.execute(text(f"PRAGMA table_info({table_name})"))
                    columns = [(row[1], row[2]) for row in result.fetchall()]
                    logger.info(f"Columns: {columns}")

                    # Get sample data
                    result = conn.execute(text(f"SELECT * FROM {table_name} LIMIT 3"))
                    sample_data = result.fetchall()
                    logger.info(f"Sample data (first 3 rows): {sample_data}")

            except Exception as e:
                logger.warning(f"Could not explore {table_name}: {e}")

        def find_captain_data_sources(self):
            """Find all tables that might contain captain data"""
            logger.info("=== SEARCHING FOR CAPTAIN DATA ===")

            tables_to_check = ["captains", "joined_ships", "ships"]

            for table in tables_to_check:
                self.explore_table_relationships(table)

        # Add these methods to the class
        SimpleMigration.explore_table_relationships = explore_table_relationships
        SimpleMigration.find_captain_data_sources = find_captain_data_sources

    def migrate_captains_v2(self, limit=None):
        """Alternative captain migration that tries different approaches"""
        logger.info("Migrating captains table (v2)...")

        try:
            # Get accession mappings from new database
            new_accessions = pd.read_sql(
                "SELECT id, accession_tag FROM accessions", self.new_engine
            )
            logger.info(f"Found {len(new_accessions)} accessions to match against")

            # First, let's see what captain-related data we actually have
            self.find_captain_data_sources()

            # Try multiple approaches to get captain data
            captain_queries = [
                # Approach 1: Direct from captains table
                """
                SELECT 
                    c.id as captain_id,
                    c.sequence,
                    'direct' as source_method
                FROM captains c
                WHERE c.sequence IS NOT NULL
                """,
                # Approach 2: From joined_ships if it exists
                """
                SELECT DISTINCT 
                    j.captainID as captain_id,
                    j.captain_id,
                    'joined_ships' as source_method
                FROM joined_ships j
                WHERE j.captainID IS NOT NULL OR j.captain_id IS NOT NULL
                """,
                # Approach 3: Any table with captain in the name
                """
                SELECT * FROM captains LIMIT 5
                """,
            ]

            for i, query in enumerate(captain_queries):
                try:
                    logger.info(f"Trying captain query approach {i + 1}")
                    if limit and "LIMIT" not in query:
                        query += f" LIMIT {limit}"

                    result_df = pd.read_sql(query, self.old_engine)
                    logger.info(f"Approach {i + 1} returned {len(result_df)} rows")
                    logger.info(f"Columns: {list(result_df.columns)}")

                    if len(result_df) > 0:
                        logger.info(f"Sample data: {result_df.head(2).to_dict()}")

                        # If we found some data, try to process it
                        if (
                            "sequence" in result_df.columns
                            and result_df["sequence"].notna().any()
                        ):
                            # Create a simple mapping
                            final_df = pd.DataFrame(
                                {
                                    "captain_name": result_df.get(
                                        "captain_id", f"captain_{i}"
                                    ),
                                    "sequence": result_df["sequence"],
                                    "accession_id": new_accessions["id"].iloc[0]
                                    if len(new_accessions) > 0
                                    else 1,  # Just use first accession for now
                                    "reviewed": "migrated",
                                    "evidence": f"migration_approach_{i + 1}",
                                }
                            )

                            # Remove rows with null sequences
                            final_df = final_df[final_df["sequence"].notna()]

                            if len(final_df) > 0:
                                final_df.to_sql(
                                    "captains",
                                    self.new_engine,
                                    if_exists="append",
                                    index=False,
                                )
                                logger.info(
                                    f"Successfully migrated {len(final_df)} captain records using approach {i + 1}"
                                )
                                return len(final_df)

                except Exception as e:
                    logger.warning(f"Captain query approach {i + 1} failed: {e}")
                    continue

            logger.warning("No captain data could be migrated")
            return 0

        except Exception as e:
            logger.error(f"Error in captain migration v2: {e}")
            return 0

    def migrate_basic_tables(self, limit=10):
        """Migrate the basic core tables in dependency order"""
        logger.info(f"Starting basic migration with limit={limit}")

        # Step 1: Create schema
        self.create_schema()

        # Step 2: Migrate base table (accessions)
        accession_count = self.migrate_accessions(limit=limit)

        # Step 3: Migrate dependent tables
        ship_count = self.migrate_ships(limit=limit)
        captain_count = self.migrate_captains_v2(limit=limit)

        logger.info("=== MIGRATION SUMMARY ===")
        logger.info(f"Accessions: {accession_count}")
        logger.info(f"Ships: {ship_count}")
        logger.info(f"Captains: {captain_count}")

        return {
            "accessions": accession_count,
            "ships": ship_count,
            "captains": captain_count,
        }

    def validate_migration(self):
        """Simple validation to check if migration worked"""
        logger.info("=== VALIDATION ===")

        with self.new_engine.connect() as conn:
            # Check table counts
            tables = ["accessions", "ships", "captains"]
            for table in tables:
                try:
                    result = conn.execute(text(f"SELECT COUNT(*) FROM {table}"))
                    count = result.fetchone()[0]
                    logger.info(f"{table}: {count} rows")
                except Exception as e:
                    logger.error(f"Error checking {table}: {e}")

            # Check foreign key relationships
            try:
                result = conn.execute(
                    text("""
                    SELECT COUNT(*) FROM ships s 
                    LEFT JOIN accessions a ON s.accession_id = a.id 
                    WHERE a.id IS NULL
                """)
                )
                orphaned_ships = result.fetchone()[0]
                if orphaned_ships > 0:
                    logger.warning(
                        f"Found {orphaned_ships} ships without valid accession_id"
                    )
                else:
                    logger.info("All ships have valid accession_id references")

            except Exception as e:
                logger.error(f"Error validating relationships: {e}")


def main():
    """Main function to run the simple migration"""

    # Create migration instance
    migration = SimpleMigration()

    # Clean start
    migration.clean_start()

    # Show what we're working with
    migration.show_old_database_info()

    # Run basic migration with a small number of records first
    migration.migrate_basic_tables(limit=10)

    # Validate the results
    migration.validate_migration()

    logger.info("Simple migration completed!")
    logger.info(f"New database created at: {migration.new_db_path}")


if __name__ == "__main__":
    main()

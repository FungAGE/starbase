#!/usr/bin/env python3
"""
Enhanced debug version of the simple migration script.
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


class SimpleMigrationDebug:
    """Enhanced debug version of simple database migration"""

    def __init__(self, old_db_path=None, new_db_path=None):
        self.old_db_path = old_db_path or DB_PATHS["starbase"]
        self.new_db_path = new_db_path or f"{DB_PATHS['starbase']}_debug"

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

    def explore_table_sample(self, table_name, limit=3):
        """Show sample data from a table"""
        try:
            with self.old_engine.connect() as conn:
                result = conn.execute(text(f"SELECT * FROM {table_name} LIMIT {limit}"))
                rows = result.fetchall()
                columns = result.keys()

                logger.info(f"=== SAMPLE DATA FROM {table_name.upper()} ===")
                for i, row in enumerate(rows):
                    row_dict = dict(zip(columns, row))
                    logger.info(f"Row {i + 1}: {row_dict}")

        except Exception as e:
            logger.warning(f"Could not get sample data from {table_name}: {e}")

    def debug_captain_data(self):
        """Debug what captain data is available in the old database"""
        logger.info("=== DEBUGGING CAPTAIN DATA ===")

        # List of tables that might contain captain data
        potential_tables = ["captains", "joined_ships", "ships"]

        for table in potential_tables:
            info = self.get_old_table_info(table)
            if info and info["count"] > 0:
                logger.info(f"\n--- {table.upper()} TABLE ---")
                logger.info(f"Columns: {info['columns']}")
                self.explore_table_sample(table, limit=2)

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

    def migrate_captains_debug(self, limit=None):
        """Debug version of captain migration that tries multiple approaches"""
        logger.info("Migrating captains table (debug version)...")

        # First show what we're working with
        self.debug_captain_data()

        try:
            # Get accession mappings from new database
            new_accessions = pd.read_sql(
                "SELECT id, accession_tag FROM accessions", self.new_engine
            )
            logger.info(f"Have {len(new_accessions)} accessions to work with")

            # Try multiple approaches
            approaches = [
                {
                    "name": "Direct captains table",
                    "query": "SELECT id, captain_name, sequence FROM captains WHERE sequence IS NOT NULL",
                },
                {
                    "name": "Captains via joined_ships",
                    "query": """
                    SELECT DISTINCT 
                        a.accession_tag,
                        COALESCE(j.captain_id, j.captainID) as captain_name,
                        c.sequence
                    FROM joined_ships j
                    INNER JOIN accessions a ON j.ship_id = a.id
                    LEFT JOIN captains c ON j.captain_id = c.id
                    WHERE c.sequence IS NOT NULL
                    """,
                },
                {
                    "name": "All non-null captain sequences",
                    "query": "SELECT * FROM captains WHERE sequence IS NOT NULL AND LENGTH(sequence) > 50",
                },
            ]

            for i, approach in enumerate(approaches):
                try:
                    logger.info(f"\n--- APPROACH {i + 1}: {approach['name']} ---")

                    query = approach["query"]
                    if limit:
                        query += f" LIMIT {limit}"

                    df = pd.read_sql(query, self.old_engine)
                    logger.info(f"Query returned {len(df)} rows")

                    if len(df) > 0:
                        logger.info(f"Columns: {list(df.columns)}")
                        logger.info(
                            f"Sample row: {df.iloc[0].to_dict() if len(df) > 0 else 'None'}"
                        )

                        # Try to create captain records
                        if "sequence" in df.columns:
                            captain_records = []

                            for idx, row in df.iterrows():
                                # Determine accession_id
                                if "accession_tag" in df.columns and pd.notna(
                                    row["accession_tag"]
                                ):
                                    # Find matching accession
                                    matching_acc = new_accessions[
                                        new_accessions["accession_tag"]
                                        == row["accession_tag"]
                                    ]
                                    if len(matching_acc) > 0:
                                        accession_id = matching_acc.iloc[0]["id"]
                                    else:
                                        accession_id = (
                                            new_accessions.iloc[0]["id"]
                                            if len(new_accessions) > 0
                                            else 1
                                        )
                                else:
                                    accession_id = (
                                        new_accessions.iloc[0]["id"]
                                        if len(new_accessions) > 0
                                        else 1
                                    )

                                # Determine captain name
                                captain_name = row.get(
                                    "captain_name", row.get("id", f"captain_{idx}")
                                )

                                captain_records.append(
                                    {
                                        "captain_name": str(captain_name),
                                        "sequence": row["sequence"],
                                        "accession_id": accession_id,
                                        "reviewed": "migrated",
                                        "evidence": approach["name"],
                                    }
                                )

                            if captain_records:
                                captain_df = pd.DataFrame(captain_records)
                                captain_df.to_sql(
                                    "captains",
                                    self.new_engine,
                                    if_exists="append",
                                    index=False,
                                )
                                logger.info(
                                    f"✓ Successfully migrated {len(captain_df)} captain records using: {approach['name']}"
                                )
                                return len(captain_df)

                except Exception as e:
                    logger.warning(f"Approach {i + 1} failed: {e}")
                    continue

            logger.warning("No captain approaches succeeded")
            return 0

        except Exception as e:
            logger.error(f"Fatal error in captain migration: {e}")
            return 0

    def run_debug_migration(self, limit=5):
        """Run a small debug migration"""
        logger.info(f"=== RUNNING DEBUG MIGRATION (limit={limit}) ===")

        # Clean start
        self.clean_start()

        # Show database info
        self.show_old_database_info()

        # Create schema
        self.create_schema()

        # Migrate tables
        results = {}
        results["accessions"] = self.migrate_accessions(limit=limit)
        results["ships"] = self.migrate_ships(limit=limit)
        results["captains"] = self.migrate_captains_debug(limit=limit)

        # Show results
        logger.info("=== DEBUG MIGRATION RESULTS ===")
        for table, count in results.items():
            logger.info(f"{table}: {count} records migrated")

        return results


def main():
    """Run the debug migration"""
    migration = SimpleMigrationDebug()
    migration.run_debug_migration(limit=5)

    logger.info("Debug migration completed!")
    logger.info(f"Debug database created at: {migration.new_db_path}")


if __name__ == "__main__":
    main()

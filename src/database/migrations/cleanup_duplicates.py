#!/usr/bin/env python3
"""
Cleanup Duplicate Records in joined_ships table

This script removes incomplete duplicate records and fixes schema constraints.
It's a Python version of the cleanup_duplicates.sql script.

Usage:
    python cleanup_duplicates.py [database_path]

Example:
    python cleanup_duplicates.py src/database/db/starbase.sqlite
"""

import sqlite3
import sys
import os
from pathlib import Path
from typing import Tuple, Dict, Any
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class DuplicateCleanupManager:
    """Manages cleanup of duplicate records in the joined_ships table."""

    def __init__(self, db_path: str):
        """
        Initialize the cleanup manager.

        Args:
            db_path: Path to the SQLite database file
        """
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database file not found: {db_path}")

        self.connection = None
        self.cursor = None

    def __enter__(self):
        """Context manager entry."""
        self.connection = sqlite3.connect(str(self.db_path))
        self.connection.execute("PRAGMA foreign_keys = ON")
        self.cursor = self.connection.cursor()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.connection:
            if exc_type is None:
                self.connection.commit()
                logger.info("Changes committed successfully")
            else:
                self.connection.rollback()
                logger.error(f"Error occurred, rolling back: {exc_val}")
            self.connection.close()

    def analyze_current_state(self) -> Dict[str, int]:
        """Analyze the current duplication problem."""
        logger.info("ANALYSIS: Current State")

        analyses = {
            "total_records": self._get_count("SELECT COUNT(*) FROM joined_ships"),
            "unique_starshipids": self._get_count(
                "SELECT COUNT(DISTINCT starshipID) FROM joined_ships"
            ),
            "records_with_ship_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NOT NULL"
            ),
            "records_without_ship_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NULL"
            ),
        }

        for key, value in analyses.items():
            logger.info(f"{key.replace('_', ' ').title()}: {value}")

        return analyses

    def create_backup(self) -> int:
        """Create a backup table with current data."""
        logger.info("Creating backup table...")

        # Drop existing backup if it exists
        self.cursor.execute("DROP TABLE IF EXISTS joined_ships_before_cleanup")

        # Create new backup
        self.cursor.execute(
            "CREATE TABLE joined_ships_before_cleanup AS SELECT * FROM joined_ships"
        )

        backup_count = self._get_count(
            "SELECT COUNT(*) FROM joined_ships_before_cleanup"
        )
        logger.info(f"Backup created with {backup_count} records")

        return backup_count

    def show_sample_incomplete_records(self, limit: int = 5) -> None:
        """Show sample of records to be removed (NULL ship_id)."""
        logger.info("Sample of records to be removed (NULL ship_id):")

        query = """
        SELECT starshipID, evidence, source, ship_family_id, ship_id, tax_id 
        FROM joined_ships 
        WHERE ship_id IS NULL 
        LIMIT ?
        """

        results = self.cursor.execute(query, (limit,)).fetchall()
        for row in results:
            logger.info(f"  {row[0]}|{row[1]}|{row[2]}|{row[3]}|{row[4]}|{row[5]}")

    def analyze_duplicates(self) -> Tuple[int, int]:
        """Identify which starshipIDs have both complete and incomplete records."""
        logger.info("Analyzing duplicate patterns...")

        # Create temporary analysis table
        self.cursor.execute("DROP TABLE IF EXISTS duplicate_analysis")
        self.cursor.execute("""
        CREATE TEMPORARY TABLE duplicate_analysis AS
        SELECT 
            starshipID,
            COUNT(*) as total_count,
            COUNT(CASE WHEN ship_id IS NOT NULL THEN 1 END) as complete_records,
            COUNT(CASE WHEN ship_id IS NULL THEN 1 END) as incomplete_records,
            MAX(CASE WHEN ship_id IS NOT NULL THEN 1 ELSE 0 END) as has_complete_record
        FROM joined_ships 
        GROUP BY starshipID
        """)

        # Count starshipIDs with both complete and incomplete records
        both_types = self._get_count("""
        SELECT COUNT(*) 
        FROM duplicate_analysis 
        WHERE complete_records > 0 AND incomplete_records > 0
        """)

        # Count starshipIDs with only incomplete records
        only_incomplete = self._get_count("""
        SELECT COUNT(*) 
        FROM duplicate_analysis 
        WHERE complete_records = 0 AND incomplete_records > 0
        """)

        logger.info(
            f"StarshipIDs with both complete and incomplete records: {both_types}"
        )
        logger.info(f"StarshipIDs with only incomplete records: {only_incomplete}")

        return both_types, only_incomplete

    def show_example_duplicates(self, starship_id: str = "altals1_s00058") -> None:
        """Show example of duplicates for a specific starshipID."""
        logger.info(f"Example of duplicates for starshipID {starship_id}:")

        query = """
        SELECT starshipID, evidence, source, ship_id, ship_family_id, curated_status
        FROM joined_ships 
        WHERE starshipID = ?
        ORDER BY ship_id DESC NULLS LAST
        """

        results = self.cursor.execute(query, (starship_id,)).fetchall()
        for row in results:
            logger.info(f"  {row[0]}|{row[1]}|{row[2]}|{row[3]}|{row[4]}|{row[5]}")

    def remove_incomplete_duplicates(self) -> int:
        """Remove incomplete records where complete records exist."""
        logger.info("Removing incomplete duplicates where complete records exist...")

        query = """
        DELETE FROM joined_ships 
        WHERE ship_id IS NULL 
          AND starshipID IN (
            SELECT starshipID 
            FROM joined_ships 
            WHERE ship_id IS NOT NULL
          )
        """

        self.cursor.execute(query)
        removed_count = self.cursor.rowcount
        remaining_count = self._get_count("SELECT COUNT(*) FROM joined_ships")

        logger.info(f"Removed {removed_count} incomplete duplicates")
        logger.info(f"Records remaining: {remaining_count}")

        return removed_count

    def keep_one_incomplete_per_starship(self) -> int:
        """For starshipIDs that only have incomplete records, keep only one record."""
        logger.info("Keeping only one incomplete record per starshipID...")

        # Create temporary table with one record per starshipID (for incomplete-only records)
        self.cursor.execute("DROP TABLE IF EXISTS keep_one_incomplete")
        self.cursor.execute("""
        CREATE TEMPORARY TABLE keep_one_incomplete AS
        SELECT MIN(rowid) as keep_rowid, starshipID
        FROM joined_ships 
        WHERE ship_id IS NULL
        GROUP BY starshipID
        """)

        # Remove all but one incomplete record for each starshipID
        query = """
        DELETE FROM joined_ships 
        WHERE ship_id IS NULL 
          AND rowid NOT IN (SELECT keep_rowid FROM keep_one_incomplete)
        """

        self.cursor.execute(query)
        removed_count = self.cursor.rowcount
        remaining_count = self._get_count("SELECT COUNT(*) FROM joined_ships")

        logger.info(f"Removed {removed_count} extra incomplete records")
        logger.info(f"Records remaining: {remaining_count}")

        return removed_count

    def final_analysis(self) -> Dict[str, int]:
        """Perform final analysis after cleanup."""
        logger.info("FINAL ANALYSIS: After Cleanup")

        analyses = {
            "total_records": self._get_count("SELECT COUNT(*) FROM joined_ships"),
            "unique_starshipids": self._get_count(
                "SELECT COUNT(DISTINCT starshipID) FROM joined_ships"
            ),
            "records_with_ship_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NOT NULL"
            ),
            "records_without_ship_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NULL"
            ),
        }

        for key, value in analyses.items():
            logger.info(f"{key.replace('_', ' ').title()}: {value}")

        return analyses

    def check_remaining_duplicates(self, limit: int = 5) -> None:
        """Check for any remaining duplicates with ship_id."""
        logger.info("Checking for duplicate ship_ids:")

        query = """
        SELECT ship_id, COUNT(*) as count 
        FROM joined_ships 
        WHERE ship_id IS NOT NULL 
        GROUP BY ship_id 
        HAVING COUNT(*) > 1
        LIMIT ?
        """

        results = self.cursor.execute(query, (limit,)).fetchall()
        if results:
            for row in results:
                logger.info(f"  ship_id {row[0]}: {row[1]} records")
        else:
            logger.info("  No duplicate ship_ids found")

    def show_cleaned_sample(self) -> None:
        """Show sample of cleaned data."""
        logger.info("Sample of cleaned data:")

        query = """
        SELECT starshipID, evidence, source, ship_id, ship_family_id, curated_status
        FROM joined_ships 
        WHERE starshipID IN ('altals1_s00058', 'altalt11_s00073')
        ORDER BY starshipID, ship_id DESC NULLS LAST
        """

        results = self.cursor.execute(query).fetchall()
        for row in results:
            logger.info(f"  {row[0]}|{row[1]}|{row[2]}|{row[3]}|{row[4]}|{row[5]}")

    def verify_data_integrity(self, backup_count: int) -> int:
        """Verify data integrity by checking records removed."""
        current_count = self._get_count("SELECT COUNT(*) FROM joined_ships")
        removed_count = backup_count - current_count

        logger.info(f"Data integrity check - records removed: {removed_count}")
        return removed_count

    def cleanup_backup(self) -> None:
        """Remove the backup table (optional)."""
        logger.info("Cleaning up backup table...")
        self.cursor.execute("DROP TABLE IF EXISTS joined_ships_before_cleanup")
        logger.info("Backup table removed")

    def run_full_cleanup(self, keep_backup: bool = True) -> Dict[str, Any]:
        """
        Run the complete cleanup process.

        Args:
            keep_backup: Whether to keep the backup table after cleanup

        Returns:
            Dictionary with cleanup results
        """
        results = {}

        # Step 1: Analyze current state
        results["initial_analysis"] = self.analyze_current_state()

        # Step 2: Create backup
        backup_count = self.create_backup()
        results["backup_count"] = backup_count

        # Step 3: Show sample incomplete records
        self.show_sample_incomplete_records()

        # Step 4 & 5: Analyze duplicates and show examples
        both_types, only_incomplete = self.analyze_duplicates()
        results["duplicate_analysis"] = {
            "both_types": both_types,
            "only_incomplete": only_incomplete,
        }
        self.show_example_duplicates()

        # Step 6: Remove incomplete duplicates
        removed_incomplete = self.remove_incomplete_duplicates()
        results["removed_incomplete"] = removed_incomplete

        # Step 7: Keep one incomplete per starshipID
        removed_extras = self.keep_one_incomplete_per_starship()
        results["removed_extras"] = removed_extras

        # Step 8: Final analysis
        results["final_analysis"] = self.final_analysis()

        # Step 9: Check remaining duplicates
        self.check_remaining_duplicates()

        # Step 10: Show cleaned sample
        self.show_cleaned_sample()

        # Step 11: Verify integrity
        total_removed = self.verify_data_integrity(backup_count)
        results["total_removed"] = total_removed

        # Step 12: Optional cleanup
        if not keep_backup:
            self.cleanup_backup()

        results["success"] = True
        return results

    def _get_count(self, query: str) -> int:
        """Helper method to get a count from a query."""
        return self.cursor.execute(query).fetchone()[0]


def main():
    """Main function to run the cleanup script."""
    if len(sys.argv) > 1:
        db_path = sys.argv[1]
    else:
        # Default path
        db_path = "src/database/db/starbase.sqlite"

    if not os.path.exists(db_path):
        logger.error(f"Database file not found: {db_path}")
        sys.exit(1)

    try:
        with DuplicateCleanupManager(db_path) as cleanup_manager:
            results = cleanup_manager.run_full_cleanup(keep_backup=True)

            logger.info("=" * 50)
            logger.info("CLEANUP COMPLETED SUCCESSFULLY!")
            logger.info("=" * 50)
            logger.info(
                f"Initial records: {results['initial_analysis']['total_records']}"
            )
            logger.info(f"Final records: {results['final_analysis']['total_records']}")
            logger.info(f"Total removed: {results['total_removed']}")
            logger.info("Backup table preserved: joined_ships_before_cleanup")

    except Exception as e:
        logger.error(f"Cleanup failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

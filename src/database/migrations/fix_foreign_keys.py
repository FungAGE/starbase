#!/usr/bin/env python3
"""
Foreign Key Migration Script

This script adds proper foreign key constraints to the joined_ships table
while maintaining data integrity. It's a Python version of the fix_foreign_keys.sql script.

Usage:
    python fix_foreign_keys.py [database_path]

Example:
    python fix_foreign_keys.py src/database/db/starbase.sqlite
"""

import sqlite3
import sys
import os
from pathlib import Path
from typing import Dict, Any, List, Tuple
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class ForeignKeyMigrationManager:
    """Manages foreign key migration for the joined_ships table."""

    def __init__(self, db_path: str):
        """
        Initialize the migration manager.

        Args:
            db_path: Path to the SQLite database file
        """
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database file not found: {db_path}")

        self.connection = None
        self.cursor = None

        # Foreign key constraints to add (excluding ships table due to PK issues)
        self.foreign_keys = [
            ("ship_family_id", "family_names", "id"),
            ("tax_id", "taxonomy", "id"),
            ("genome_id", "genomes", "id"),
            ("captain_id", "captains", "id"),
            ("ship_navis_id", "navis_names", "id"),
            ("ship_haplotype_id", "haplotype_names", "id"),
        ]

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
                logger.info("Migration committed successfully")
            else:
                self.connection.rollback()
                logger.error(f"Error occurred, rolling back: {exc_val}")
            self.connection.close()

    def create_backup(self) -> int:
        """Create backup of current data."""
        logger.info("Creating backup of current joined_ships data...")

        # Drop existing backup if it exists
        self.cursor.execute("DROP TABLE IF EXISTS joined_ships_backup")

        # Create new backup
        self.cursor.execute(
            "CREATE TABLE joined_ships_backup AS SELECT * FROM joined_ships"
        )

        backup_count = self._get_count("SELECT COUNT(*) FROM joined_ships_backup")
        logger.info(f"Backup created with {backup_count} records")

        return backup_count

    def check_current_state(self) -> int:
        """Check current state of joined_ships table."""
        record_count = self._get_count("SELECT COUNT(*) FROM joined_ships")
        logger.info(f"Current joined_ships records: {record_count}")
        return record_count

    def create_new_table_with_fk(self) -> None:
        """Create new table with proper foreign key constraints."""
        logger.info("Creating new table with foreign key constraints...")

        # Drop existing table if it exists
        self.cursor.execute("DROP TABLE IF EXISTS joined_ships_with_fk")

        # Create new table with foreign key constraints
        create_table_sql = """
        CREATE TABLE joined_ships_with_fk (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            starshipID TEXT NOT NULL,
            evidence TEXT,
            source TEXT,
            curated_status TEXT,
            
            -- Foreign key columns
            ship_family_id INTEGER,
            tax_id INTEGER,
            ship_id INTEGER,  -- Keep as regular column (ships table lacks proper PK)
            genome_id INTEGER,
            captain_id INTEGER,
            ship_navis_id INTEGER,
            ship_haplotype_id INTEGER,
            
            -- Timestamps
            created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
            
            -- Foreign key constraints
            CONSTRAINT fk_joined_ships_family 
                FOREIGN KEY (ship_family_id) REFERENCES family_names(id),
            CONSTRAINT fk_joined_ships_taxonomy 
                FOREIGN KEY (tax_id) REFERENCES taxonomy(id),
            CONSTRAINT fk_joined_ships_genome 
                FOREIGN KEY (genome_id) REFERENCES genomes(id),
            CONSTRAINT fk_joined_ships_captain 
                FOREIGN KEY (captain_id) REFERENCES captains(id),
            CONSTRAINT fk_joined_ships_navis 
                FOREIGN KEY (ship_navis_id) REFERENCES navis_names(id),
            CONSTRAINT fk_joined_ships_haplotype 
                FOREIGN KEY (ship_haplotype_id) REFERENCES haplotype_names(id)
        )
        """

        self.cursor.execute(create_table_sql)
        logger.info("New table created successfully")

    def create_indexes(self) -> None:
        """Create performance indexes on the new table."""
        logger.info("Creating performance indexes...")

        indexes = [
            "CREATE INDEX idx_joined_ships_final_starshipid ON joined_ships_with_fk(starshipID)",
            "CREATE INDEX idx_joined_ships_final_family ON joined_ships_with_fk(ship_family_id)",
            "CREATE INDEX idx_joined_ships_final_taxonomy ON joined_ships_with_fk(tax_id)",
            "CREATE INDEX idx_joined_ships_final_ship ON joined_ships_with_fk(ship_id)",
            "CREATE INDEX idx_joined_ships_final_genome ON joined_ships_with_fk(genome_id)",
            "CREATE INDEX idx_joined_ships_final_captain ON joined_ships_with_fk(captain_id)",
            "CREATE INDEX idx_joined_ships_final_navis ON joined_ships_with_fk(ship_navis_id)",
            "CREATE INDEX idx_joined_ships_final_haplotype ON joined_ships_with_fk(ship_haplotype_id)",
        ]

        for index_sql in indexes:
            self.cursor.execute(index_sql)

        logger.info(f"Created {len(indexes)} performance indexes")

    def check_foreign_key_violations(self) -> Dict[str, int]:
        """Check for foreign key violations before migration."""
        logger.info("Checking for foreign key violations before migration...")

        violations = {}

        # Check each foreign key relationship
        violation_checks = [
            ("family_id", "ship_family_id", "family_names"),
            ("taxonomy_id", "tax_id", "taxonomy"),
            ("genome_id", "genome_id", "genomes"),
            ("captain_id", "captain_id", "captains"),
            ("navis_id", "ship_navis_id", "navis_names"),
            ("haplotype_id", "ship_haplotype_id", "haplotype_names"),
        ]

        for check_name, column_name, ref_table in violation_checks:
            if column_name == "genome_id":
                # Special handling for genome_id conversion
                query = f"""
                SELECT COUNT(*) 
                FROM joined_ships 
                WHERE {column_name} IS NOT NULL 
                  AND CAST({column_name} AS INTEGER) NOT IN (SELECT id FROM {ref_table})
                """
            else:
                query = f"""
                SELECT COUNT(*) 
                FROM joined_ships 
                WHERE {column_name} IS NOT NULL 
                  AND {column_name} NOT IN (SELECT id FROM {ref_table})
                """

            violation_count = self._get_count(query)
            violations[check_name] = violation_count
            logger.info(f"Invalid {check_name} references: {violation_count}")

        return violations

    def migrate_data(self) -> int:
        """Migrate data from old table to new table."""
        logger.info("Migrating data to new table...")

        migrate_sql = """
        INSERT INTO joined_ships_with_fk (
            starshipID, evidence, source, curated_status,
            ship_family_id, tax_id, ship_id, captain_id,
            ship_navis_id, ship_haplotype_id, genome_id
        )
        SELECT 
            starshipID, evidence, source, curated_status,
            ship_family_id, tax_id, ship_id, captain_id,
            ship_navis_id, ship_haplotype_id,
            -- Handle genome_id conversion
            CASE 
                WHEN genome_id IS NULL THEN NULL
                WHEN genome_id = '' THEN NULL
                WHEN genome_id GLOB '[0-9]*' THEN CAST(genome_id AS INTEGER)
                ELSE NULL 
            END as genome_id
        FROM joined_ships
        WHERE starshipID IS NOT NULL
        """

        self.cursor.execute(migrate_sql)
        migrated_count = self.cursor.rowcount
        logger.info(f"Migrated {migrated_count} records to new table")

        return migrated_count

    def verify_migration(self, original_count: int) -> Tuple[int, int]:
        """Verify migration success."""
        logger.info("Verifying migration...")

        original_records = self._get_count("SELECT COUNT(*) FROM joined_ships")
        migrated_records = self._get_count("SELECT COUNT(*) FROM joined_ships_with_fk")

        logger.info(f"Original records: {original_records}")
        logger.info(f"Migrated records: {migrated_records}")

        if original_records != migrated_records:
            logger.warning(
                f"Record count mismatch! Original: {original_records}, Migrated: {migrated_records}"
            )

        return original_records, migrated_records

    def replace_old_table(self) -> None:
        """Replace the old table with the new one."""
        logger.info("Replacing old table with new foreign key table...")

        # Rename old table to backup
        self.cursor.execute("ALTER TABLE joined_ships RENAME TO joined_ships_pre_fk")

        # Rename new table to original name
        self.cursor.execute("ALTER TABLE joined_ships_with_fk RENAME TO joined_ships")

        logger.info("Table replacement completed")

    def verify_foreign_keys(self) -> List[Tuple]:
        """Verify that foreign key constraints are in place."""
        logger.info("Verifying foreign key constraints...")

        foreign_keys = self.cursor.execute(
            "PRAGMA foreign_key_list(joined_ships)"
        ).fetchall()

        logger.info("Foreign key constraints:")
        for fk in foreign_keys:
            logger.info(f"  {fk[3]} -> {fk[2]}({fk[4]})")

        return foreign_keys

    def final_record_counts(self) -> Dict[str, int]:
        """Get final record counts for verification."""
        logger.info("Final record counts:")

        counts = {
            "total_records": self._get_count("SELECT COUNT(*) FROM joined_ships"),
            "records_with_genome_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE genome_id IS NOT NULL"
            ),
            "records_with_ship_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NOT NULL"
            ),
            "records_with_family_id": self._get_count(
                "SELECT COUNT(*) FROM joined_ships WHERE ship_family_id IS NOT NULL"
            ),
        }

        for key, value in counts.items():
            logger.info(f"  {key.replace('_', ' ').title()}: {value}")

        return counts

    def test_relationships(self, limit: int = 5) -> List[Tuple]:
        """Test foreign key relationships with sample data."""
        logger.info("Testing sample foreign key relationships:")

        test_query = """
        SELECT 
            js.starshipID,
            js.ship_family_id,
            fn.familyName,
            js.ship_navis_id,
            nn.navis_name,
            js.genome_id,
            g.ome
        FROM joined_ships js
        LEFT JOIN family_names fn ON js.ship_family_id = fn.id
        LEFT JOIN navis_names nn ON js.ship_navis_id = nn.id
        LEFT JOIN genomes g ON js.genome_id = g.id
        WHERE js.ship_family_id IS NOT NULL 
        LIMIT ?
        """

        results = self.cursor.execute(test_query, (limit,)).fetchall()

        for row in results:
            logger.info(
                f"  {row[0]} -> family: {row[2]}, navis: {row[4]}, genome: {row[6]}"
            )

        return results

    def cleanup_backup_tables(self) -> None:
        """Clean up backup tables (optional)."""
        logger.info("Cleaning up backup tables...")

        backup_tables = ["joined_ships_backup", "joined_ships_pre_fk"]

        for table in backup_tables:
            self.cursor.execute(f"DROP TABLE IF EXISTS {table}")

        logger.info("Backup tables cleaned up")

    def run_full_migration(self, keep_backups: bool = True) -> Dict[str, Any]:
        """
        Run the complete foreign key migration process.

        Args:
            keep_backups: Whether to keep backup tables after migration

        Returns:
            Dictionary with migration results
        """
        results = {}

        # Step 1: Create backup
        backup_count = self.create_backup()
        results["backup_count"] = backup_count

        # Step 2: Check current state
        original_count = self.check_current_state()
        results["original_count"] = original_count

        # Step 3: Create new table with foreign keys
        self.create_new_table_with_fk()

        # Step 4: Create indexes
        self.create_indexes()

        # Step 5: Check for violations
        violations = self.check_foreign_key_violations()
        results["violations"] = violations

        # Check if there are any violations
        if any(count > 0 for count in violations.values()):
            raise ValueError(f"Foreign key violations found: {violations}")

        # Step 6: Migrate data
        migrated_count = self.migrate_data()
        results["migrated_count"] = migrated_count

        # Step 7: Verify migration
        original_verify, migrated_verify = self.verify_migration(original_count)
        results["verification"] = {
            "original": original_verify,
            "migrated": migrated_verify,
        }

        # Step 8: Replace old table
        self.replace_old_table()

        # Step 9: Verify foreign keys
        foreign_keys = self.verify_foreign_keys()
        results["foreign_keys"] = foreign_keys

        # Step 10: Final counts
        final_counts = self.final_record_counts()
        results["final_counts"] = final_counts

        # Step 11: Test relationships
        sample_relationships = self.test_relationships()
        results["sample_relationships"] = sample_relationships

        # Step 12: Optional cleanup
        if not keep_backups:
            self.cleanup_backup_tables()

        results["success"] = True
        results["note"] = (
            "ship_id column kept as regular INTEGER (ships table lacks proper PRIMARY KEY)"
        )

        return results

    def _get_count(self, query: str) -> int:
        """Helper method to get a count from a query."""
        return self.cursor.execute(query).fetchone()[0]


def main():
    """Main function to run the foreign key migration script."""
    if len(sys.argv) > 1:
        db_path = sys.argv[1]
    else:
        # Default path
        db_path = "src/database/db/starbase.sqlite"

    if not os.path.exists(db_path):
        logger.error(f"Database file not found: {db_path}")
        sys.exit(1)

    try:
        with ForeignKeyMigrationManager(db_path) as migration_manager:
            results = migration_manager.run_full_migration(keep_backups=True)

            logger.info("=" * 60)
            logger.info("FOREIGN KEY MIGRATION COMPLETED SUCCESSFULLY!")
            logger.info("=" * 60)
            logger.info(f"Original records: {results['original_count']}")
            logger.info(f"Migrated records: {results['migrated_count']}")
            logger.info(f"Foreign key constraints: {len(results['foreign_keys'])}")
            logger.info(
                f"Final record count: {results['final_counts']['total_records']}"
            )
            logger.info(f"Note: {results['note']}")
            logger.info("Backup tables preserved for safety")

    except Exception as e:
        logger.error(f"Migration failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

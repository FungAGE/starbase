#!/usr/bin/env python3
"""
Database normalization script for starbase.sqlite
This script normalizes the joined_ships table to use foreign key references
instead of storing navis and haplotype values directly.
"""

import sqlite3
import sys
from pathlib import Path
from typing import Tuple


class DatabaseNormalizer:
    """Handles database normalization operations for starbase.sqlite"""

    def __init__(self, db_path: str):
        """Initialize with database path"""
        self.db_path = Path(db_path)
        if not self.db_path.exists():
            raise FileNotFoundError(f"Database file not found: {db_path}")

        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row  # Enable column access by name

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.conn:
            self.conn.close()

    def populate_navis_names(self) -> int:
        """Populate navis_names table with unique values from joined_ships"""
        print("Step 1: Populating navis_names table...")

        cursor = self.conn.cursor()

        # Insert unique navis values
        cursor.execute("""
            INSERT INTO navis_names (navis_name)
            SELECT DISTINCT starship_navis 
            FROM joined_ships 
            WHERE starship_navis IS NOT NULL
            ORDER BY starship_navis
        """)

        count = cursor.rowcount
        print(f"  → Inserted {count} unique navis records")
        return count

    def populate_haplotype_names(self) -> int:
        """Populate haplotype_names table with unique values from joined_ships"""
        print("Step 2: Populating haplotype_names table...")

        cursor = self.conn.cursor()

        # Insert unique haplotype values
        cursor.execute("""
            INSERT INTO haplotype_names (haplotype_name)
            SELECT DISTINCT starship_haplotype 
            FROM joined_ships 
            WHERE starship_haplotype IS NOT NULL
            ORDER BY starship_haplotype
        """)

        count = cursor.rowcount
        print(f"  → Inserted {count} unique haplotype records")
        return count

    def add_foreign_key_columns(self) -> None:
        """Add new foreign key columns to joined_ships table"""
        print("Step 3: Adding foreign key columns to joined_ships table...")

        cursor = self.conn.cursor()

        # Check if columns already exist
        cursor.execute("PRAGMA table_info(joined_ships)")
        columns = [row[1] for row in cursor.fetchall()]

        if "navis_id" not in columns:
            cursor.execute("ALTER TABLE joined_ships ADD COLUMN navis_id INTEGER")
            print("  → Added navis_id column")
        else:
            print("  → navis_id column already exists")

        if "haplotype_id" not in columns:
            cursor.execute("ALTER TABLE joined_ships ADD COLUMN haplotype_id INTEGER")
            print("  → Added haplotype_id column")
        else:
            print("  → haplotype_id column already exists")

    def update_foreign_key_references(self) -> Tuple[int, int]:
        """Update the new foreign key columns with references to lookup tables"""
        print("Step 4: Updating foreign key references...")

        cursor = self.conn.cursor()

        # Update navis_id references
        cursor.execute("""
            UPDATE joined_ships 
            SET navis_id = (
                SELECT n.id 
                FROM navis_names n 
                WHERE n.navis_name = joined_ships.starship_navis
            )
            WHERE starship_navis IS NOT NULL
        """)
        navis_updates = cursor.rowcount
        print(f"  → Updated {navis_updates} navis_id references")

        # Update haplotype_id references
        cursor.execute("""
            UPDATE joined_ships 
            SET haplotype_id = (
                SELECT h.id 
                FROM haplotype_names h 
                WHERE h.haplotype_name = joined_ships.starship_haplotype
            )
            WHERE starship_haplotype IS NOT NULL
        """)
        haplotype_updates = cursor.rowcount
        print(f"  → Updated {haplotype_updates} haplotype_id references")

        return navis_updates, haplotype_updates

    def create_indexes(self) -> None:
        """Create indexes for performance optimization"""
        print("Step 5: Creating performance indexes...")

        cursor = self.conn.cursor()

        # Create indexes (ignore if they already exist)
        try:
            cursor.execute(
                "CREATE INDEX idx_joined_ships_navis_id ON joined_ships(navis_id)"
            )
            print("  → Created navis_id index")
        except sqlite3.OperationalError as e:
            if "already exists" in str(e):
                print("  → navis_id index already exists")
            else:
                raise

        try:
            cursor.execute(
                "CREATE INDEX idx_joined_ships_haplotype_id ON joined_ships(haplotype_id)"
            )
            print("  → Created haplotype_id index")
        except sqlite3.OperationalError as e:
            if "already exists" in str(e):
                print("  → haplotype_id index already exists")
            else:
                raise

    def verify_data_integrity(self) -> dict:
        """Verify the data integrity after normalization"""
        print("Step 6: Verifying data integrity...")

        cursor = self.conn.cursor()

        # Get counts for verification
        cursor.execute("SELECT COUNT(*) FROM navis_names")
        navis_count = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM haplotype_names")
        haplotype_count = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM joined_ships WHERE navis_id IS NOT NULL")
        navis_refs = cursor.fetchone()[0]

        cursor.execute(
            "SELECT COUNT(*) FROM joined_ships WHERE haplotype_id IS NOT NULL"
        )
        haplotype_refs = cursor.fetchone()[0]

        results = {
            "navis_records": navis_count,
            "haplotype_records": haplotype_count,
            "navis_references": navis_refs,
            "haplotype_references": haplotype_refs,
        }

        print(f"  → Navis records populated: {navis_count}")
        print(f"  → Haplotype records populated: {haplotype_count}")
        print(f"  → Joined_ships with navis_id: {navis_refs}")
        print(f"  → Joined_ships with haplotype_id: {haplotype_refs}")

        return results

    def show_sample_mappings(self, limit: int = 5) -> None:
        """Show sample mappings to verify correctness"""
        print(f"\nSample mappings (limit {limit}):")

        cursor = self.conn.cursor()

        # Sample navis mappings
        print("\nNavis mappings:")
        cursor.execute(
            """
            SELECT js.starshipID, js.starship_navis, n.navis_name, js.navis_id
            FROM joined_ships js
            LEFT JOIN navis_names n ON js.navis_id = n.id
            WHERE js.starship_navis IS NOT NULL
            LIMIT ?
        """,
            (limit,),
        )

        for row in cursor.fetchall():
            print(f"  {row[0]}: '{row[1]}' → '{row[2]}' (id: {row[3]})")

        # Sample haplotype mappings
        print("\nHaplotype mappings:")
        cursor.execute(
            """
            SELECT js.starshipID, js.starship_haplotype, h.haplotype_name, js.haplotype_id
            FROM joined_ships js
            LEFT JOIN haplotype_names h ON js.haplotype_id = h.id
            WHERE js.starship_haplotype IS NOT NULL
            LIMIT ?
        """,
            (limit,),
        )

        for row in cursor.fetchall():
            print(f"  {row[0]}: '{row[1]}' → '{row[2]}' (id: {row[3]})")

    def normalize_database(self, verify: bool = True) -> dict:
        """
        Perform complete database normalization

        Args:
            verify: Whether to run verification steps

        Returns:
            Dictionary with operation results
        """
        print("Starting database normalization...\n")

        try:
            # Begin transaction
            self.conn.execute("BEGIN TRANSACTION")

            # Perform normalization steps
            navis_inserted = self.populate_navis_names()
            haplotype_inserted = self.populate_haplotype_names()
            self.add_foreign_key_columns()
            navis_updates, haplotype_updates = self.update_foreign_key_references()
            self.create_indexes()

            if verify:
                verification_results = self.verify_data_integrity()
                self.show_sample_mappings()
            else:
                verification_results = {}

            # Commit transaction
            self.conn.commit()
            print("\n✅ Database normalization completed successfully!")

            return {
                "navis_inserted": navis_inserted,
                "haplotype_inserted": haplotype_inserted,
                "navis_updates": navis_updates,
                "haplotype_updates": haplotype_updates,
                "verification": verification_results,
            }

        except Exception as e:
            # Rollback on error
            self.conn.rollback()
            print(f"\n❌ Error during normalization: {e}")
            raise


def main():
    """Main function to run the normalization"""
    # Default database path
    default_db_path = "src/database/db/starbase.sqlite"

    # Get database path from command line argument or use default
    db_path = sys.argv[1] if len(sys.argv) > 1 else default_db_path

    try:
        with DatabaseNormalizer(db_path) as normalizer:
            results = normalizer.normalize_database()

            print("\n📊 Summary:")
            print(f"  • Navis records inserted: {results['navis_inserted']}")
            print(f"  • Haplotype records inserted: {results['haplotype_inserted']}")
            print(f"  • Navis references updated: {results['navis_updates']}")
            print(f"  • Haplotype references updated: {results['haplotype_updates']}")

    except Exception as e:
        print(f"❌ Failed to normalize database: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

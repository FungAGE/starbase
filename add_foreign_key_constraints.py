#!/usr/bin/env python3
"""
Add proper foreign key constraints to joined_ships table, specifically for accession_id.

This script recreates the joined_ships table with proper foreign key constraints.
Since SQLite doesn't support adding FK constraints to existing columns, we need to:
1. Create a new table with the correct constraints
2. Copy data from the old table
3. Drop the old table and rename the new one
"""

import sqlite3
import sys
import os
from datetime import datetime
import argparse

# Add the project root to sys.path so we can import our modules
PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.append(PROJECT_ROOT)

from src.config.logging import get_logger

logger = get_logger(__name__)


def add_foreign_key_constraints(db_path: str, dry_run: bool = True):
    """
    Add proper foreign key constraints to joined_ships table.
    
    Args:
        db_path: Path to SQLite database
        dry_run: If True, only show what would be changed
    """
    logger.info(f"Adding foreign key constraints to joined_ships table (dry_run={dry_run})")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Disable foreign keys during migration
        cursor.execute("PRAGMA foreign_keys = OFF")
        
        # Check current foreign key constraints and columns
        logger.info("Checking current foreign key constraints and columns...")
        cursor.execute("PRAGMA foreign_key_list(joined_ships)")
        current_fks = cursor.fetchall()
        accession_fk_exists = any(fk[2] == 'accessions' and fk[3] == 'accession_id' for fk in current_fks)

        cursor.execute("PRAGMA table_info(joined_ships)")
        cols = [c[1] for c in cursor.fetchall()]
        ship_id_exists = 'ship_id' in cols

        if accession_fk_exists and ship_id_exists:
            logger.info("Foreign key for accession_id exists and ship_id column is present")
            if dry_run:
                logger.info("No changes needed")
            return {"status": "no_changes_needed", "message": "FK present and ship_id exists"}

        if not accession_fk_exists:
            logger.info("accession_id foreign key constraint missing - will recreate table")
        if not ship_id_exists:
            logger.info("ship_id column missing - will recreate table to add ship_id with FK")
        
        if not dry_run:
            # Step 1: Create new table with proper foreign key constraints
            logger.info("Creating new joined_ships table with foreign key constraints...")
            
            # Drop the new table if it exists from a previous failed attempt
            cursor.execute("DROP TABLE IF EXISTS joined_ships_new")

            # Ensure referenced target column is eligible (PRIMARY KEY or UNIQUE)
            try:
                cursor.execute("PRAGMA table_info(ships)")
                ship_cols = cursor.fetchall()
                id_is_pk = any(col[1] == 'id' and col[5] == 1 for col in ship_cols)
                if not id_is_pk:
                    # Create a unique index on ships.id if not present
                    cursor.execute(
                        "CREATE UNIQUE INDEX IF NOT EXISTS ux_ships_id ON ships(id)"
                    )
                    logger.info("Ensured UNIQUE index on ships.id for FK target")
            except Exception as e:
                logger.warning(f"Could not verify/ensure ships.id uniqueness: {e}")
            
            create_new_table_sql = """
            CREATE TABLE joined_ships_new (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                starshipID TEXT NOT NULL,
                evidence TEXT,
                source TEXT,
                curated_status TEXT,
                
                -- Direct link to accession (when sequence data is available)
                -- Multiple ships can reference the same accession (duplicates allowed)
                accession_id INTEGER,
                
                -- Links to classification and annotation data
                ship_id INTEGER,
                ship_family_id INTEGER,
                tax_id INTEGER,
                genome_id INTEGER,
                captain_id INTEGER,
                ship_navis_id INTEGER,
                ship_haplotype_id INTEGER,
                created_at DATETIME,
                updated_at DATETIME,
                
                -- Foreign key constraints
                CONSTRAINT fk_joined_ships_accession 
                    FOREIGN KEY (accession_id) REFERENCES accessions(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_ship 
                    FOREIGN KEY (ship_id) REFERENCES ships(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_family 
                    FOREIGN KEY (ship_family_id) REFERENCES family_names(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_taxonomy 
                    FOREIGN KEY (tax_id) REFERENCES taxonomy(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_genome 
                    FOREIGN KEY (genome_id) REFERENCES genomes(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_captain 
                    FOREIGN KEY (captain_id) REFERENCES captains(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_navis 
                    FOREIGN KEY (ship_navis_id) REFERENCES navis_names(id)
                    ON DELETE SET NULL ON UPDATE CASCADE,
                CONSTRAINT fk_joined_ships_haplotype 
                    FOREIGN KEY (ship_haplotype_id) REFERENCES haplotype_names(id)
                    ON DELETE SET NULL ON UPDATE CASCADE
            )
            """
            
            cursor.execute(create_new_table_sql)
            logger.info("Created joined_ships_new table with foreign key constraints")
            
            # Step 2: Copy data from old table to new table
            logger.info("Copying data from old table to new table...")
            
            # Handle cases where old table may or may not have ship_id column
            cursor.execute("PRAGMA table_info(joined_ships)")
            old_cols = [c[1] for c in cursor.fetchall()]
            has_old_ship_id = 'ship_id' in old_cols

            if has_old_ship_id:
                # Cleanse invalid ship_id values during copy to avoid FK violations
                copy_data_sql = """
                INSERT INTO joined_ships_new (
                    id, starshipID, evidence, source, curated_status,
                    accession_id, ship_id, ship_family_id, tax_id, genome_id, captain_id,
                    ship_navis_id, ship_haplotype_id, created_at, updated_at
                )
                SELECT 
                    js.id, js.starshipID, js.evidence, js.source, js.curated_status,
                    js.accession_id,
                    CASE WHEN s.id IS NULL THEN NULL ELSE js.ship_id END as ship_id,
                    js.ship_family_id, js.tax_id, js.genome_id, js.captain_id,
                    js.ship_navis_id, js.ship_haplotype_id, js.created_at, js.updated_at
                FROM joined_ships js
                LEFT JOIN ships s ON js.ship_id = s.id
                """
            else:
                copy_data_sql = """
                INSERT INTO joined_ships_new (
                    id, starshipID, evidence, source, curated_status,
                    accession_id, ship_id, ship_family_id, tax_id, genome_id, captain_id,
                    ship_navis_id, ship_haplotype_id, created_at, updated_at
                )
                SELECT 
                    id, starshipID, evidence, source, curated_status,
                    accession_id, NULL as ship_id, ship_family_id, tax_id, genome_id, captain_id,
                    ship_navis_id, ship_haplotype_id, created_at, updated_at
                FROM joined_ships
                """
            
            cursor.execute(copy_data_sql)
            rows_copied = cursor.rowcount
            logger.info(f"Copied {rows_copied} rows to new table")
            
            # Step 3: Drop old table and rename new table
            logger.info("Replacing old table with new table...")
            cursor.execute("DROP TABLE joined_ships")
            cursor.execute("ALTER TABLE joined_ships_new RENAME TO joined_ships")
            
            # Step 4: Recreate indexes
            logger.info("Recreating indexes...")
            
            indexes = [
                "CREATE INDEX idx_joined_ships_final_starshipid ON joined_ships(starshipID)",
                "CREATE INDEX idx_joined_ships_final_family ON joined_ships(ship_family_id)",
                "CREATE INDEX idx_joined_ships_final_taxonomy ON joined_ships(tax_id)",
                "CREATE INDEX idx_joined_ships_final_genome ON joined_ships(genome_id)",
                "CREATE INDEX idx_joined_ships_final_captain ON joined_ships(captain_id)",
                "CREATE INDEX idx_joined_ships_final_navis ON joined_ships(ship_navis_id)",
                "CREATE INDEX idx_joined_ships_final_haplotype ON joined_ships(ship_haplotype_id)",
                "CREATE INDEX idx_joined_ships_curated ON joined_ships(curated_status)",
                "CREATE INDEX idx_joined_ships_accession_id ON joined_ships(accession_id)"
            ]
            
            for index_sql in indexes:
                try:
                    cursor.execute(index_sql)
                except sqlite3.Error as e:
                    logger.warning(f"Could not create index: {e}")
            
            logger.info("Recreated indexes")
            
            # Step 5: Re-enable foreign keys and check integrity for joined_ships only
            cursor.execute("PRAGMA foreign_keys = ON")
            
            # Check foreign key integrity specifically for joined_ships table
            cursor.execute("PRAGMA foreign_key_check(joined_ships)")
            fk_violations = cursor.fetchall()
            if fk_violations:
                logger.error(f"Found {len(fk_violations)} foreign key violations in joined_ships after migration:")
                for violation in fk_violations[:10]:  # Show first 10
                    logger.error(f"  {violation}")
                raise Exception("Foreign key violations found in joined_ships - migration failed")
            else:
                logger.info("No foreign key violations found in joined_ships - migration successful")
            
            # Verify the new foreign key constraints
            cursor.execute("PRAGMA foreign_key_list(joined_ships)")
            new_fks = cursor.fetchall()
            logger.info(f"New foreign key constraints: {len(new_fks)}")
            for fk in new_fks:
                logger.info(f"  {fk[3]} -> {fk[2]}.{fk[4]}")
            
            conn.commit()
            logger.info("Foreign key constraints added successfully!")
            
            return {
                "status": "success",
                "rows_migrated": rows_copied,
                "foreign_keys_added": len(new_fks),
                "message": "Table recreated with proper foreign key constraints"
            }
        
        else:
            # Dry run
            logger.info("Would recreate joined_ships table with proper foreign key constraints:")
            logger.info("  - accession_id -> accessions(id)")
            logger.info("  - ship_family_id -> family_names(id)")
            logger.info("  - tax_id -> taxonomy(id)")
            logger.info("  - genome_id -> genomes(id)")
            logger.info("  - captain_id -> captains(id)")
            logger.info("  - ship_navis_id -> navis_names(id)")
            logger.info("  - ship_haplotype_id -> haplotype_names(id)")
            logger.info("Would also recreate all indexes")
            
            return {
                "status": "dry_run",
                "message": "Would recreate table with foreign key constraints"
            }
        
    except Exception as e:
        logger.error(f"Error adding foreign key constraints: {str(e)}")
        if not dry_run:
            conn.rollback()
        raise
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Add foreign key constraints to joined_ships table")
    parser.add_argument("--db-path", required=True, help="Path to SQLite database file")
    parser.add_argument("--dry-run", action="store_true", default=True, 
                       help="Only show what would be changed (default: True)")
    parser.add_argument("--execute", action="store_true", 
                       help="Actually perform the changes (overrides --dry-run)")
    
    args = parser.parse_args()
    
    # Determine if this is a dry run
    dry_run = not args.execute
    
    if not dry_run:
        response = input("Are you sure you want to recreate the joined_ships table? This is a significant operation. (yes/no): ")
        if response.lower() != 'yes':
            print("Operation cancelled.")
            sys.exit(0)
    
    try:
        result = add_foreign_key_constraints(
            db_path=args.db_path,
            dry_run=dry_run
        )
        
        print(f"\nOperation completed {'(dry run)' if dry_run else 'successfully'}!")
        print(f"Status: {result['status']}")
        print(f"Message: {result['message']}")
        
    except Exception as e:
        logger.error(f"Operation failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()

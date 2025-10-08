#!/usr/bin/env python3
"""
Fix the Ships.id primary key issue by recreating the ships table with proper primary key constraint.

The ships table currently has its id column without a proper PRIMARY KEY constraint,
which is causing foreign key reference issues.
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


def fix_ships_primary_key(db_path: str, dry_run: bool = True):
    """
    Fix the ships table primary key by recreating it with proper constraints.
    
    Args:
        db_path: Path to SQLite database
        dry_run: If True, only show what would be changed
    """
    logger.info(f"Fixing ships table primary key (dry_run={dry_run})")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Check current ships table structure
        logger.info("Checking current ships table structure...")
        cursor.execute("PRAGMA table_info(ships)")
        current_cols = cursor.fetchall()
        
        # Check if id is already a primary key
        id_is_pk = any(col[1] == 'id' and col[5] == 1 for col in current_cols)
        if id_is_pk:
            logger.info("ships.id is already a primary key - no changes needed")
            return {"status": "no_changes_needed", "message": "Primary key already exists"}
        
        logger.info("ships.id is not a primary key - will recreate table")
        
        # Get current data count
        cursor.execute("SELECT COUNT(*) FROM ships")
        row_count = cursor.fetchone()[0]
        logger.info(f"Current ships table has {row_count} rows")
        
        if not dry_run:
            # Disable foreign keys during migration
            cursor.execute("PRAGMA foreign_keys = OFF")
            
            # Step 1: Create new ships table with proper primary key
            logger.info("Creating new ships table with proper primary key...")
            
            # Drop the new table if it exists from a previous failed attempt
            cursor.execute("DROP TABLE IF EXISTS ships_new")
            
            create_new_table_sql = """
            CREATE TABLE ships_new (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sequence TEXT,
                md5 TEXT,
                sequence_length INTEGER,
                header VARCHAR,
                accession_id INTEGER,
                rev_comp_md5 VARCHAR,
                
                -- Foreign key constraint
                CONSTRAINT fk_ships_accession 
                    FOREIGN KEY (accession_id) REFERENCES accessions(id)
                    ON DELETE SET NULL ON UPDATE CASCADE
            )
            """
            
            cursor.execute(create_new_table_sql)
            logger.info("Created ships_new table with proper primary key")
            
            # Step 2: Copy data from old table to new table (removing duplicates)
            logger.info("Copying data from old table to new table (removing duplicates)...")
            
            # First, check for and log duplicates
            cursor.execute("SELECT id, COUNT(*) as count FROM ships GROUP BY id HAVING COUNT(*) > 1")
            duplicates = cursor.fetchall()
            if duplicates:
                logger.warning(f"Found {len(duplicates)} duplicate IDs in ships table:")
                for dup_id, count in duplicates:
                    logger.warning(f"  ID {dup_id}: {count} duplicates")
            
            # First copy records with existing IDs
            logger.info("Copying records with existing IDs...")
            copy_existing_sql = """
            INSERT INTO ships_new (
                id, sequence, md5, sequence_length, header, accession_id, rev_comp_md5
            )
            SELECT 
                id, sequence, md5, sequence_length, header, accession_id, rev_comp_md5
            FROM ships
            WHERE id IS NOT NULL
            ORDER BY id
            """
            cursor.execute(copy_existing_sql)
            existing_count = cursor.rowcount
            logger.info(f"Copied {existing_count} records with existing IDs")
            
            # Then copy records with NULL IDs, assigning new sequential IDs
            logger.info("Copying records with NULL IDs and assigning new IDs...")
            copy_null_sql = """
            INSERT INTO ships_new (
                id, sequence, md5, sequence_length, header, accession_id, rev_comp_md5
            )
            SELECT 
                (SELECT COALESCE(MAX(id), 0) + ROW_NUMBER() OVER (ORDER BY rowid)) as id,
                sequence, md5, sequence_length, header, accession_id, rev_comp_md5
            FROM ships
            WHERE id IS NULL
            ORDER BY rowid
            """
            cursor.execute(copy_null_sql)
            null_count = cursor.rowcount
            logger.info(f"Copied {null_count} records with NULL IDs (assigned new IDs)")
            
            rows_copied = existing_count + null_count
            logger.info(f"Total copied {rows_copied} rows to new table")
            
            # Step 3: Drop old table and rename new table
            logger.info("Replacing old table with new table...")
            cursor.execute("DROP TABLE ships")
            cursor.execute("ALTER TABLE ships_new RENAME TO ships")
            
            # Step 4: Recreate indexes
            logger.info("Recreating indexes...")
            
            indexes = [
                "CREATE INDEX idx_ships_accession_id ON ships(accession_id)"
            ]
            
            for index_sql in indexes:
                try:
                    cursor.execute(index_sql)
                except sqlite3.Error as e:
                    logger.warning(f"Could not create index: {e}")
            
            logger.info("Recreated indexes")
            
            # Step 5: Re-enable foreign keys and check integrity
            cursor.execute("PRAGMA foreign_keys = ON")
            
            # Check foreign key integrity
            cursor.execute("PRAGMA foreign_key_check(ships)")
            fk_violations = cursor.fetchall()
            if fk_violations:
                logger.error(f"Found {len(fk_violations)} foreign key violations in ships after migration:")
                for violation in fk_violations[:10]:  # Show first 10
                    logger.error(f"  {violation}")
                raise Exception("Foreign key violations found in ships - migration failed")
            else:
                logger.info("No foreign key violations found in ships - migration successful")
            
            # Verify the new primary key
            cursor.execute("PRAGMA table_info(ships)")
            new_cols = cursor.fetchall()
            id_pk_status = any(col[1] == 'id' and col[5] == 1 for col in new_cols)
            if id_pk_status:
                logger.info("ships.id is now properly configured as primary key")
            else:
                raise Exception("Failed to set ships.id as primary key")
            
            conn.commit()
            logger.info("Ships table primary key fixed successfully!")
            
            return {
                "status": "success",
                "rows_migrated": rows_copied,
                "message": "Ships table recreated with proper primary key"
            }
        
        else:
            # Dry run
            logger.info("Would recreate ships table with proper primary key:")
            logger.info("  - id INTEGER PRIMARY KEY AUTOINCREMENT")
            logger.info("  - accession_id -> accessions(id) foreign key")
            logger.info("  - Recreate indexes")
            
            return {
                "status": "dry_run",
                "message": "Would recreate ships table with proper primary key"
            }
        
    except Exception as e:
        logger.error(f"Error fixing ships primary key: {str(e)}")
        if not dry_run:
            conn.rollback()
        raise
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Fix ships table primary key")
    parser.add_argument("--db-path", required=True, help="Path to SQLite database file")
    parser.add_argument("--dry-run", action="store_true", default=True, 
                       help="Only show what would be changed (default: True)")
    parser.add_argument("--execute", action="store_true", 
                       help="Actually perform the changes (overrides --dry-run)")
    
    args = parser.parse_args()
    
    # Determine if this is a dry run
    dry_run = not args.execute
    
    if not dry_run:
        response = input("Are you sure you want to recreate the ships table? This is a significant operation. (yes/no): ")
        if response.lower() != 'yes':
            print("Operation cancelled.")
            sys.exit(0)
    
    try:
        result = fix_ships_primary_key(
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

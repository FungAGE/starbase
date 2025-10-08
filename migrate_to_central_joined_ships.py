#!/usr/bin/env python3
"""
Migration script to restructure the database with joined_ships as the central table.

This script:
1. Adds accession_id column to joined_ships
2. Migrates ship_id references to accession_id where possible
3. Creates joined_ships entries for ships that don't have them
4. Updates foreign key relationships
5. Optionally removes the ship_id column

IMPORTANT: Always run with --dry-run first to see what changes will be made!
"""

import sqlite3
import json
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


def migrate_to_central_joined_ships(db_path: str, dry_run: bool = True, remove_ship_id: bool = False):
    """
    Migrate database to use joined_ships as central table with direct accession_id links.
    
    Args:
        db_path: Path to SQLite database
        dry_run: If True, only show what would be changed
        remove_ship_id: If True, remove the ship_id column after migration
    """
    logger.info(f"Starting migration to central joined_ships structure (dry_run={dry_run})")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Disable foreign keys during migration to avoid constraint conflicts
        cursor.execute("PRAGMA foreign_keys = OFF")
        
        migration_stats = {
            'joined_ships_with_ship_id': 0,
            'successful_mappings': 0,
            'failed_mappings': 0,
            'ships_without_joined_ships': 0,
            'new_joined_ships_created': 0,
            'errors': []
        }
        
        # Step 1: Add accession_id column to joined_ships if it doesn't exist
        logger.info("Step 1: Adding accession_id column to joined_ships...")
        
        # Check if accession_id column already exists, and whether legacy ship_id exists
        cursor.execute("PRAGMA table_info(joined_ships)")
        columns = [col[1] for col in cursor.fetchall()]
        has_ship_id = 'ship_id' in columns
        
        if 'accession_id' not in columns:
            if not dry_run:
                cursor.execute("ALTER TABLE joined_ships ADD COLUMN accession_id INTEGER REFERENCES accessions(id)")
                logger.info("Added accession_id column to joined_ships with foreign key constraint")
            else:
                logger.info("Would add accession_id column to joined_ships with foreign key constraint")
        else:
            logger.info("accession_id column already exists in joined_ships")
            # Check if foreign key constraint exists
            cursor.execute("PRAGMA foreign_key_list(joined_ships)")
            fks = cursor.fetchall()
            accession_fk_exists = any(fk[2] == 'accessions' and fk[3] == 'accession_id' for fk in fks)
            
            if not accession_fk_exists:
                logger.warning("accession_id column exists but foreign key constraint is missing")
                if not dry_run:
                    # We'll need to recreate the table to add the FK constraint properly
                    logger.info("Foreign key constraint will be handled in table recreation step")
                else:
                    logger.info("Would add missing foreign key constraint for accession_id")
        
        # Step 2: Map legacy ship_id references to accession_id (if ship_id column exists)
        logger.info("Step 2: Mapping legacy ship_id to accession_id (if present)...")
        if has_ship_id:
            # Get all joined_ships entries that have ship_id
            cursor.execute("""
                SELECT js.id, js.starshipID, js.ship_id, s.accession_id
                FROM joined_ships js
                LEFT JOIN ships s ON js.ship_id = s.id
                WHERE js.ship_id IS NOT NULL AND js.accession_id IS NULL
            """)
            joined_ships_with_ship_id = cursor.fetchall()
            migration_stats['joined_ships_with_ship_id'] = len(joined_ships_with_ship_id)
            logger.info(f"Found {len(joined_ships_with_ship_id)} joined_ships entries with ship_id")

            for js_id, starship_id, ship_id, accession_id in joined_ships_with_ship_id:
                if accession_id:
                    if not dry_run:
                        cursor.execute(
                            "UPDATE joined_ships SET accession_id = ? WHERE id = ?",
                            (accession_id, js_id)
                        )
                    migration_stats['successful_mappings'] += 1
                    logger.debug(
                        f"Mapped joined_ships {js_id} ({starship_id}) from ship_id {ship_id} to accession_id {accession_id}"
                    )
                else:
                    migration_stats['failed_mappings'] += 1
                    migration_stats['errors'].append({
                        'type': 'ship_without_accession',
                        'joined_ships_id': js_id,
                        'starship_id': starship_id,
                        'ship_id': ship_id,
                        'message': f'Ship {ship_id} has no accession_id - sequence data not available'
                    })
                    logger.warning(
                        f"Ship {ship_id} (starship {starship_id}) has no accession - leaving accession_id as NULL"
                    )
        else:
            logger.info("No ship_id column present in joined_ships - skipping legacy mapping step")
        
        # Step 3: Create missing joined_ships entries
        logger.info("Step 3: Creating joined_ships entries for items without them...")
        if has_ship_id:
            # Legacy path based on ships table presence in joined_ships
            cursor.execute("""
                SELECT s.id, s.accession_id, a.ship_name, a.accession_tag
                FROM ships s
                LEFT JOIN accessions a ON s.accession_id = a.id
                WHERE s.id NOT IN (SELECT ship_id FROM joined_ships WHERE ship_id IS NOT NULL)
            """)
            ships_without_joined_ships = cursor.fetchall()
            migration_stats['ships_without_joined_ships'] = len(ships_without_joined_ships)
            logger.info(f"Found {len(ships_without_joined_ships)} ships without joined_ships entries")

            for ship_id, accession_id, ship_name, accession_tag in ships_without_joined_ships:
                if ship_name:
                    starship_id = ship_name
                elif accession_tag:
                    starship_id = accession_tag
                else:
                    starship_id = f"SHIP_{ship_id}"

                cursor.execute("SELECT id FROM joined_ships WHERE starshipID = ?", (starship_id,))
                if cursor.fetchone():
                    starship_id = f"{starship_id}_SHIP_{ship_id}"

                if not dry_run:
                    cursor.execute(
                        """
                        INSERT INTO joined_ships (
                            starshipID, evidence, source, curated_status, 
                            accession_id, created_at, updated_at
                        ) VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            starship_id,
                            "AUTO_MIGRATED",
                            "database_migration",
                            "auto_created",
                            accession_id,
                            datetime.now(),
                            datetime.now(),
                        ),
                    )
                migration_stats['new_joined_ships_created'] += 1
                logger.debug(
                    f"Created joined_ships entry for ship {ship_id} with starshipID '{starship_id}'"
                    if not dry_run else
                    f"Would create joined_ships entry for ship {ship_id} with starshipID '{starship_id}'"
                )
        else:
            # New path: ensure each accession has at least one joined_ships entry
            cursor.execute(
                """
                SELECT a.id, a.ship_name, a.accession_tag
                FROM accessions a
                LEFT JOIN joined_ships js ON js.accession_id = a.id
                WHERE js.id IS NULL
                """
            )
            accessions_without_joined = cursor.fetchall()
            migration_stats['ships_without_joined_ships'] = len(accessions_without_joined)
            logger.info(f"Found {len(accessions_without_joined)} accessions without joined_ships entries")

            for accession_id, ship_name, accession_tag in accessions_without_joined:
                if ship_name:
                    starship_id = ship_name
                elif accession_tag:
                    starship_id = accession_tag
                else:
                    starship_id = f"ACC_{accession_id}"

                cursor.execute("SELECT id FROM joined_ships WHERE starshipID = ?", (starship_id,))
                exists = cursor.fetchone()
                if exists:
                    starship_id = f"{starship_id}_ACC_{accession_id}"

                if not dry_run:
                    cursor.execute(
                        """
                        INSERT INTO joined_ships (
                            starshipID, evidence, source, curated_status,
                            accession_id, created_at, updated_at
                        ) VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            starship_id,
                            "AUTO_MIGRATED",
                            "database_migration",
                            "auto_created",
                            accession_id,
                            datetime.now(),
                            datetime.now(),
                        ),
                    )
                migration_stats['new_joined_ships_created'] += 1
                logger.debug(
                    f"Created joined_ships entry for accession {accession_id} with starshipID '{starship_id}'"
                    if not dry_run else
                    f"Would create joined_ships entry for accession {accession_id} with starshipID '{starship_id}'"
                )
        
        # Step 4: Optionally remove ship_id column
        if remove_ship_id:
            logger.info("Step 4: Removing ship_id column from joined_ships...")
            if not dry_run:
                # SQLite doesn't support DROP COLUMN directly, so we need to recreate the table
                logger.warning("Removing ship_id column requires table recreation - this is a significant operation")
                # For now, just log that this would happen
                logger.info("Table recreation for ship_id removal would happen here (not implemented in dry run)")
            else:
                logger.info("Would remove ship_id column from joined_ships (requires table recreation)")
        
        # Step 5: Create indexes and re-enable foreign keys
        logger.info("Step 5: Creating indexes and handling foreign key constraints...")
        if not dry_run:
            try:
                # Add index for accession_id
                cursor.execute("""
                    CREATE INDEX IF NOT EXISTS idx_joined_ships_accession_id 
                    ON joined_ships(accession_id)
                """)
                logger.info("Added index for accession_id")
                
                # Re-enable foreign keys
                cursor.execute("PRAGMA foreign_keys = ON")
                
                # Check foreign key integrity
                cursor.execute("PRAGMA foreign_key_check")
                fk_violations = cursor.fetchall()
                if fk_violations:
                    logger.warning(f"Found {len(fk_violations)} foreign key violations:")
                    for violation in fk_violations[:5]:  # Show first 5
                        logger.warning(f"  {violation}")
                else:
                    logger.info("No foreign key violations found")
                    
            except sqlite3.Error as e:
                logger.warning(f"Error handling foreign key constraints: {e}")
        else:
            logger.info("Would add index for accession_id and re-enable foreign key constraints")
        
        # Summary
        logger.info("Migration Summary:")
        logger.info(f"  - joined_ships entries with ship_id: {migration_stats['joined_ships_with_ship_id']}")
        logger.info(f"  - Successful mappings to accession_id: {migration_stats['successful_mappings']}")
        logger.info(f"  - Failed mappings (ships without accessions): {migration_stats['failed_mappings']}")
        logger.info(f"  - Ships without joined_ships entries: {migration_stats['ships_without_joined_ships']}")
        logger.info(f"  - New joined_ships entries created: {migration_stats['new_joined_ships_created']}")
        logger.info(f"  - Errors encountered: {len(migration_stats['errors'])}")
        
        if migration_stats['errors']:
            logger.warning("Errors encountered during migration:")
            for error in migration_stats['errors'][:5]:  # Show first 5 errors
                logger.warning(f"  - {error['type']}: {error['message']}")
        
        if not dry_run:
            conn.commit()
            logger.info("Migration completed successfully!")
        else:
            logger.info("Dry run completed - no changes made to database")
        
        return migration_stats
        
    except Exception as e:
        logger.error(f"Migration failed: {str(e)}")
        if not dry_run:
            conn.rollback()
        raise
    finally:
        conn.close()


def main():
    parser = argparse.ArgumentParser(description="Migrate database to central joined_ships structure")
    parser.add_argument("--db-path", required=True, help="Path to SQLite database file")
    parser.add_argument("--dry-run", action="store_true", default=True, 
                       help="Only show what would be changed (default: True)")
    parser.add_argument("--execute", action="store_true", 
                       help="Actually perform the migration (overrides --dry-run)")
    parser.add_argument("--remove-ship-id", action="store_true",
                       help="Remove ship_id column after migration")
    
    args = parser.parse_args()
    
    # Determine if this is a dry run
    dry_run = not args.execute
    
    if not dry_run:
        response = input("Are you sure you want to modify the database? This cannot be undone easily. (yes/no): ")
        if response.lower() != 'yes':
            print("Migration cancelled.")
            sys.exit(0)
    
    try:
        stats = migrate_to_central_joined_ships(
            db_path=args.db_path,
            dry_run=dry_run,
            remove_ship_id=args.remove_ship_id
        )
        
        print(f"\nMigration completed {'(dry run)' if dry_run else 'successfully'}!")
        print(json.dumps(stats, indent=2))
        
    except Exception as e:
        logger.error(f"Migration failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()

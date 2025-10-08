#!/usr/bin/env python3
"""
Script to add foreign key constraints to SQLite database to match the schema.
SQLite requires recreating tables to add foreign key constraints.
"""

import sys
import os
import sqlite3
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from config.settings import DB_PATHS

def backup_database(db_path):
    """Create a backup of the database before migration"""
    backup_path = f"{db_path}.backup"
    print(f"Creating backup: {backup_path}")
    
    # Copy the database file
    import shutil
    shutil.copy2(db_path, backup_path)
    return backup_path

def migrate_foreign_keys(db_path, dry_run=False):
    """Add foreign key constraints to match the schema"""
    
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
    
    backup_path = backup_database(db_path)
    
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = OFF")  # Disable FK checks during migration
    
    try:
        cursor = conn.cursor()
        
        # Start transaction
        cursor.execute("BEGIN TRANSACTION")
        
        print("1. Updating joined_ships table to fix ship_id foreign key...")
        
        # Check current joined_ships structure
        cursor.execute("PRAGMA table_info(joined_ships)")
        current_columns = cursor.fetchall()
        print(f"Current joined_ships columns: {[col[1] for col in current_columns]}")
        
        # Create new joined_ships table with correct foreign keys
        create_new_joined_ships = """
        CREATE TABLE joined_ships_new (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            starshipID TEXT NOT NULL,
            evidence TEXT,
            source TEXT,
            curated_status TEXT,
            
            -- Foreign key columns
            ship_family_id INTEGER,
            tax_id INTEGER,
            ship_id INTEGER,  -- This should reference ships.id, not accessions.id
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
            CONSTRAINT fk_joined_ships_ship 
                FOREIGN KEY (ship_id) REFERENCES ships(id),
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
        
        if not dry_run:
            cursor.execute(create_new_joined_ships)
        print("Created new joined_ships table with foreign key constraints")
        
        # Copy data from old table to new table
        # First, we need to map accession IDs to ship IDs for existing records
        print("2. Mapping accession IDs to ship IDs...")
        
        # Get the mapping of accession.id -> ship.id
        cursor.execute("""
            SELECT a.id as accession_id, s.id as ship_id 
            FROM accessions a 
            LEFT JOIN ships s ON s.accession_id = a.id
        """)
        accession_to_ship = dict(cursor.fetchall())
        print(f"Found {len(accession_to_ship)} accession->ship mappings")
        
        # For accessions without ships, we need to create ship records
        cursor.execute("""
            SELECT a.id as accession_id 
            FROM accessions a 
            LEFT JOIN ships s ON s.accession_id = a.id 
            WHERE s.id IS NULL
        """)
        accessions_without_ships = [row[0] for row in cursor.fetchall()]
        
        if accessions_without_ships:
            print(f"3. Creating {len(accessions_without_ships)} missing ship records...")
            for acc_id in accessions_without_ships:
                if not dry_run:
                    cursor.execute("""
                        INSERT INTO ships (accession_id) VALUES (?)
                    """, (acc_id,))
                    ship_id = cursor.lastrowid
                    accession_to_ship[acc_id] = ship_id
        
        # Copy data to new joined_ships table, mapping accession_id to ship_id
        print("4. Copying data to new joined_ships table...")
        cursor.execute("SELECT * FROM joined_ships")
        old_records = cursor.fetchall()
        
        # Get column names to build insert statement
        cursor.execute("PRAGMA table_info(joined_ships)")
        old_columns = [col[1] for col in cursor.fetchall()]
        
        copied_count = 0
        for record in old_records:
            record_dict = dict(zip(old_columns, record))
            
            # Map the ship_id from accession to actual ship
            old_ship_id = record_dict.get('ship_id')  # This was accession_id
            new_ship_id = accession_to_ship.get(old_ship_id) if old_ship_id else None
            
            if new_ship_id:
                # Build insert statement for new table
                insert_data = {
                    'starshipID': record_dict.get('starshipID'),
                    'evidence': record_dict.get('evidence'),
                    'source': record_dict.get('source'),
                    'curated_status': record_dict.get('curated_status'),
                    'ship_family_id': record_dict.get('ship_family_id'),
                    'tax_id': record_dict.get('tax_id'),
                    'ship_id': new_ship_id,  # Use the actual ship ID
                    'genome_id': record_dict.get('genome_id'),
                    'captain_id': record_dict.get('captain_id'),
                    'ship_navis_id': record_dict.get('ship_navis_id'),
                    'ship_haplotype_id': record_dict.get('ship_haplotype_id'),
                    'created_at': record_dict.get('created_at'),
                    'updated_at': record_dict.get('updated_at'),
                }
                
                # Remove None values and build query
                clean_data = {k: v for k, v in insert_data.items() if v is not None}
                columns = list(clean_data.keys())
                values = list(clean_data.values())
                placeholders = ','.join(['?'] * len(values))
                
                if not dry_run:
                    cursor.execute(f"""
                        INSERT INTO joined_ships_new ({','.join(columns)}) 
                        VALUES ({placeholders})
                    """, values)
                copied_count += 1
                
        print(f"Copied {copied_count} records to new joined_ships table")
        
        # Add other foreign key constraints to existing tables
        print("5. Adding foreign key constraints to other tables...")
        
        # Update navis_names to add foreign key to family_names
        cursor.execute("""
            CREATE TABLE navis_names_new (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                navis_name TEXT,
                previous_navis_name TEXT,
                ship_family_id INTEGER,
                
                CONSTRAINT fk_navis_family 
                    FOREIGN KEY (ship_family_id) REFERENCES family_names(id)
            )
        """)
        
        if not dry_run:
            cursor.execute("INSERT INTO navis_names_new SELECT * FROM navis_names")
        
        # Update haplotype_names to add foreign keys
        cursor.execute("""
            CREATE TABLE haplotype_names_new (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                haplotype_name TEXT,
                previous_haplotype_name TEXT,
                navis_id INTEGER,
                ship_family_id INTEGER,
                
                CONSTRAINT fk_haplotype_navis 
                    FOREIGN KEY (navis_id) REFERENCES navis_names(id),
                CONSTRAINT fk_haplotype_family 
                    FOREIGN KEY (ship_family_id) REFERENCES family_names(id)
            )
        """)
        
        if not dry_run:
            cursor.execute("INSERT INTO haplotype_names_new SELECT * FROM haplotype_names")
        
        if not dry_run:
            # Drop old tables and rename new ones
            print("6. Replacing old tables with new ones...")
            cursor.execute("DROP TABLE joined_ships")
            cursor.execute("ALTER TABLE joined_ships_new RENAME TO joined_ships")
            
            cursor.execute("DROP TABLE navis_names")
            cursor.execute("ALTER TABLE navis_names_new RENAME TO navis_names")
            
            cursor.execute("DROP TABLE haplotype_names")
            cursor.execute("ALTER TABLE haplotype_names_new RENAME TO haplotype_names")
            
            # Recreate indexes
            print("7. Recreating indexes...")
            indexes = [
                "CREATE INDEX idx_joined_ships_final_starshipid ON joined_ships(starshipID)",
                "CREATE INDEX idx_joined_ships_final_family ON joined_ships(ship_family_id)",
                "CREATE INDEX idx_joined_ships_final_taxonomy ON joined_ships(tax_id)",
                "CREATE INDEX idx_joined_ships_final_ship ON joined_ships(ship_id)",
                "CREATE INDEX idx_joined_ships_final_genome ON joined_ships(genome_id)",
                "CREATE INDEX idx_joined_ships_final_captain ON joined_ships(captain_id)",
                "CREATE INDEX idx_joined_ships_final_navis ON joined_ships(ship_navis_id)",
                "CREATE INDEX idx_joined_ships_final_haplotype ON joined_ships(ship_haplotype_id)",
                "CREATE INDEX idx_joined_ships_curated ON joined_ships(curated_status)",
            ]
            
            for index_sql in indexes:
                try:
                    cursor.execute(index_sql)
                except sqlite3.Error as e:
                    print(f"Warning: Could not create index: {e}")
        
        # Commit transaction
        if not dry_run:
            cursor.execute("COMMIT")
            print("Migration completed successfully!")
        else:
            cursor.execute("ROLLBACK")
            print("Dry run completed - no changes made")
            
    except Exception as e:
        cursor.execute("ROLLBACK")
        print(f"Error during migration: {e}")
        raise
    finally:
        conn.execute("PRAGMA foreign_keys = ON")  # Re-enable FK checks
        conn.close()

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Migrate SQLite database to add foreign key constraints')
    parser.add_argument('--dry-run', action='store_true', 
                       help='Perform a dry run without making changes')
    parser.add_argument('--db-path', default=None,
                       help='Path to database file (default: use settings)')
    
    args = parser.parse_args()
    
    db_path = args.db_path or DB_PATHS["starbase"]
    
    if not os.path.exists(db_path):
        print(f"Error: Database file not found: {db_path}")
        sys.exit(1)
    
    print(f"Migrating database: {db_path}")
    migrate_foreign_keys(db_path, dry_run=args.dry_run)

if __name__ == "__main__":
    main()

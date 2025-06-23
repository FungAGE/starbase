#!/usr/bin/env python3
"""
Migration script to add and populate accession_id column in the gff table

Current gff table structure:
- accession (VARCHAR) - contains values like "SBS000528"
- ship_id (INTEGER) - foreign key to ships table
- Other GFF columns (source, type, start, end, etc.)

Goal:
- Add accession_id column to link to accessions.id
- Populate it by matching gff.accession with accessions.accession_tag
"""

import sqlite3
import pandas as pd
import os
from src.config.settings import DB_PATHS


def connect_to_database():
    """Connect to the starbase database"""
    db_path = DB_PATHS["starbase"]
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database not found at {db_path}")
    return sqlite3.connect(db_path)


def backup_database(conn):
    """Create a backup of the database before making changes"""
    print("Creating database backup...")
    backup_path = DB_PATHS["starbase"] + ".gff_migration_backup"

    if os.path.exists(backup_path):
        os.remove(backup_path)

    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def analyze_current_gff_state(conn):
    """Analyze the current state of the gff table"""
    print("\n=== GFF TABLE PRE-MIGRATION ANALYSIS ===")

    cursor = conn.cursor()

    # Get total gff records
    cursor.execute("SELECT COUNT(*) FROM gff")
    total_gff = cursor.fetchone()[0]
    print(f"Total gff records: {total_gff:,}")

    # Check if accession_id column already exists
    cursor.execute("PRAGMA table_info(gff)")
    columns = [col[1] for col in cursor.fetchall()]
    has_accession_id = "accession_id" in columns
    print(f"accession_id column exists: {has_accession_id}")

    # Check how many gff.accession values can match accessions.accession_tag
    cursor.execute("""
        SELECT COUNT(*) FROM gff g
        JOIN accessions a ON g.accession = a.accession_tag
    """)
    potential_matches = cursor.fetchone()[0]
    print(
        f"Potential gff.accession → accessions.accession_tag matches: {potential_matches:,}"
    )

    # Check unique accession values in gff
    cursor.execute("SELECT COUNT(DISTINCT accession) FROM gff")
    unique_accessions = cursor.fetchone()[0]
    print(f"Unique accession values in gff: {unique_accessions:,}")

    # Check for any gff records that can't be matched
    cursor.execute("""
        SELECT COUNT(*) FROM gff g
        LEFT JOIN accessions a ON g.accession = a.accession_tag
        WHERE a.accession_tag IS NULL
    """)
    unmatched_gff = cursor.fetchone()[0]
    print(f"GFF records that cannot be matched to accessions: {unmatched_gff:,}")

    if unmatched_gff > 0:
        print("\nSample unmatched gff.accession values:")
        cursor.execute("""
            SELECT DISTINCT g.accession
            FROM gff g
            LEFT JOIN accessions a ON g.accession = a.accession_tag
            WHERE a.accession_tag IS NULL
            LIMIT 10
        """)
        unmatched_samples = cursor.fetchall()
        for sample in unmatched_samples:
            print(f"  {sample[0]}")

    return total_gff, potential_matches, unmatched_gff


def step1_add_accession_id_column(conn):
    """Step 1: Add accession_id column to gff table"""
    print("\n=== STEP 1: ADD accession_id COLUMN TO gff ===")

    cursor = conn.cursor()

    # Check if column already exists
    cursor.execute("PRAGMA table_info(gff)")
    columns = [col[1] for col in cursor.fetchall()]

    if "accession_id" not in columns:
        cursor.execute("ALTER TABLE gff ADD COLUMN accession_id INTEGER")
        conn.commit()
        print("Added accession_id column to gff table")

        # Create index for better performance
        cursor.execute(
            "CREATE INDEX IF NOT EXISTS idx_gff_accession_id ON gff(accession_id)"
        )
        conn.commit()
        print("Created index on gff.accession_id")
    else:
        print("accession_id column already exists in gff table")


def step2_populate_accession_id(conn):
    """Step 2: Populate accession_id by matching gff.accession with accessions.accession_tag"""
    print("\n=== STEP 2: POPULATE gff.accession_id ===")

    # Find matches between gff.accession and accessions.accession_tag
    query = """
    SELECT 
        g.rowid,
        g.accession,
        a.id as accession_id
    FROM gff g
    JOIN accessions a ON g.accession = a.accession_tag
    WHERE g.accession_id IS NULL
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches):,} accession_id matches to populate")

    if len(matches) > 0:
        cursor = conn.cursor()
        batch_size = 5000  # Process in larger batches since this could be a lot of data
        updated_count = 0

        for i in range(0, len(matches), batch_size):
            batch = matches.iloc[i : i + batch_size]

            for _, row in batch.iterrows():
                cursor.execute(
                    "UPDATE gff SET accession_id = ? WHERE rowid = ?",
                    (row["accession_id"], row["rowid"]),
                )
                updated_count += 1

                if updated_count % 5000 == 0:
                    print(
                        f"Updated {updated_count:,} gff.accession_id relationships..."
                    )

        conn.commit()
        print(
            f"Successfully populated {updated_count:,} gff.accession_id relationships"
        )

        # Verify the update
        cursor.execute("SELECT COUNT(*) FROM gff WHERE accession_id IS NOT NULL")
        populated_count = cursor.fetchone()[0]
        print(f"Total gff records with accession_id: {populated_count:,}")

    else:
        print("No accession_id relationships found to populate")


def step3_validate_ship_id_consistency(conn):
    """Step 3: Validate that ship_id relationships are consistent with accession_id"""
    print("\n=== STEP 3: VALIDATE ship_id CONSISTENCY ===")

    cursor = conn.cursor()

    # Check if ship_id values are consistent with accession_id
    # (i.e., gff.ship_id should reference ships that belong to the same accession)
    cursor.execute("""
        SELECT COUNT(*) FROM gff g
        JOIN ships s ON g.ship_id = s.id
        WHERE g.accession_id IS NOT NULL 
        AND s.accession_id IS NOT NULL
        AND g.accession_id = s.accession_id
    """)
    consistent_relationships = cursor.fetchone()[0]
    print(
        f"GFF records with consistent ship_id ↔ accession_id relationships: {consistent_relationships:,}"
    )

    # Check for inconsistent relationships
    cursor.execute("""
        SELECT COUNT(*) FROM gff g
        JOIN ships s ON g.ship_id = s.id
        WHERE g.accession_id IS NOT NULL 
        AND s.accession_id IS NOT NULL
        AND g.accession_id != s.accession_id
    """)
    inconsistent_relationships = cursor.fetchone()[0]

    if inconsistent_relationships > 0:
        print(
            f"WARNING: {inconsistent_relationships:,} GFF records have inconsistent ship_id ↔ accession_id relationships"
        )

        # Show a sample of inconsistent records
        print("Sample of inconsistent records:")
        sample_query = """
        SELECT g.accession, g.ship_id, g.accession_id, s.accession_id as ship_accession_id
        FROM gff g
        JOIN ships s ON g.ship_id = s.id
        WHERE g.accession_id IS NOT NULL 
        AND s.accession_id IS NOT NULL
        AND g.accession_id != s.accession_id
        LIMIT 5
        """
        inconsistent_sample = pd.read_sql_query(sample_query, conn)
        print(inconsistent_sample.to_string(index=False))
    else:
        print("✅ All ship_id relationships are consistent with accession_id")


def validate_migration(conn):
    """Validate that the migration worked correctly"""
    print("\n=== MIGRATION VALIDATION ===")

    cursor = conn.cursor()

    # Check gff -> accessions relationship
    cursor.execute("""
        SELECT 
            COUNT(*) as total_gff,
            COUNT(g.accession_id) as with_accession_id,
            COUNT(*) - COUNT(g.accession_id) as missing_accession_id
        FROM gff g
    """)

    gff_result = cursor.fetchone()
    print(
        f"gff → accessions: {gff_result[1]:,}/{gff_result[0]:,} complete ({gff_result[2]:,} missing)"
    )

    # Check for valid relationships
    cursor.execute("""
        SELECT COUNT(*) FROM gff g
        JOIN accessions a ON g.accession_id = a.id
    """)

    valid_relationships = cursor.fetchone()[0]
    print(f"Valid gff.accession_id relationships: {valid_relationships:,}")

    # Summary statistics
    cursor.execute(
        "SELECT COUNT(DISTINCT accession_id) FROM gff WHERE accession_id IS NOT NULL"
    )
    unique_accession_ids = cursor.fetchone()[0]
    print(f"Unique accessions referenced by gff: {unique_accession_ids:,}")


def step4_add_primary_key_column(conn):
    """Step 4: Add id column as primary key"""
    print("\n=== STEP 4: ADD PRIMARY KEY COLUMN ===")

    cursor = conn.cursor()

    # Check if id column already exists
    cursor.execute("PRAGMA table_info(gff)")
    columns = [col[1] for col in cursor.fetchall()]

    if "id" not in columns:
        # Since SQLite doesn't allow adding a PRIMARY KEY to an existing table,
        # we need to recreate the table with the new structure
        print("Recreating table with id PRIMARY KEY column...")

        # Create new table with desired structure
        cursor.execute("""
            CREATE TABLE gff_new (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                accession_id INTEGER,
                source VARCHAR,
                type VARCHAR,
                start INTEGER,
                end INTEGER,
                phase VARCHAR,
                strand VARCHAR,
                score VARCHAR,
                attributes VARCHAR,
                ship_id INTEGER,
                FOREIGN KEY (accession_id) REFERENCES accessions(id),
                FOREIGN KEY (ship_id) REFERENCES ships(id)
            )
        """)

        # Copy data from old table to new table (excluding accession column)
        cursor.execute("""
            INSERT INTO gff_new (accession_id, source, type, start, end, phase, strand, score, attributes, ship_id)
            SELECT accession_id, source, type, start, end, phase, strand, score, attributes, ship_id
            FROM gff
        """)

        # Drop old table
        cursor.execute("DROP TABLE gff")

        # Rename new table
        cursor.execute("ALTER TABLE gff_new RENAME TO gff")

        # Recreate indexes
        cursor.execute(
            "CREATE INDEX IF NOT EXISTS idx_gff_accession_id ON gff(accession_id)"
        )
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_gff_ship_id ON gff(ship_id)")

        conn.commit()
        print(
            "✅ Successfully added id PRIMARY KEY column and dropped accession column"
        )

        # Verify the new structure
        cursor.execute("PRAGMA table_info(gff)")
        new_schema = cursor.fetchall()
        print("\nNew gff table structure:")
        for col in new_schema:
            pk_status = "PRIMARY KEY" if col[5] == 1 else ""
            print(f"  {col[1]} ({col[2]}) {pk_status}")

    else:
        print("id column already exists in gff table")


def step5_drop_accession_column(conn):
    """Step 5: Drop the accession column (this is handled in step4 during table recreation)"""
    print("\n=== STEP 5: ACCESSION COLUMN REMOVAL ===")

    cursor = conn.cursor()
    cursor.execute("PRAGMA table_info(gff)")
    columns = [col[1] for col in cursor.fetchall()]

    if "accession" in columns:
        print(
            "⚠️  WARNING: accession column still exists - this should have been removed in step 4"
        )
        print("Manual intervention may be required")
    else:
        print("✅ accession column successfully removed")


def main():
    """Main function to run the gff table migration"""

    try:
        conn = connect_to_database()

        # Create backup
        backup_database(conn)

        # Analyze current state
        analyze_current_gff_state(conn)

        # Run migration steps
        step1_add_accession_id_column(conn)
        step2_populate_accession_id(conn)
        step3_validate_ship_id_consistency(conn)
        step4_add_primary_key_column(conn)
        step5_drop_accession_column(conn)

        # Validate migration
        validate_migration(conn)

        conn.close()
        print("\n=== GFF TABLE MIGRATION COMPLETED ===")
        print("✅ Added id PRIMARY KEY column")
        print("✅ Added accession_id foreign key to accessions table")
        print("✅ Removed accession column")
        print("✅ The gff table now has proper structure and relationships!")

    except Exception as e:
        print(f"Error during migration: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

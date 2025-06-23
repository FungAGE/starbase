#!/usr/bin/env python3
"""
Schema Migration Plan: Fix relationships to make accessions the primary table

Current Issues:
1. joined_ships.ship_id links to ships.id, but relationship is broken
2. joined_ships.starshipID should link to accessions.accession_tag
3. ships.accession contains numeric IDs, not accession tags
4. No direct foreign key from joined_ships to accessions

Proposed Solution:
1. Add accession_id column to joined_ships table
2. Populate accession_id by matching starshipID to accessions.accession_tag
3. Fix ships.accession_id to properly reference accessions.id
4. Update joined_ships.ship_id to reference ships that are properly linked to accessions
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
    backup_path = DB_PATHS["starbase"] + ".schema_migration_backup"

    if os.path.exists(backup_path):
        os.remove(backup_path)

    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def analyze_current_state(conn):
    """Analyze the current state before migration"""
    print("\n=== PRE-MIGRATION ANALYSIS ===")

    cursor = conn.cursor()

    # Check joined_ships -> accessions potential matches
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js
        JOIN accessions a ON js.starshipID = a.ship_name
    """)
    potential_accession_matches = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM joined_ships")
    total_joined_ships = cursor.fetchone()[0]

    print(
        f"joined_ships records that can match accessions: {potential_accession_matches}/{total_joined_ships}"
    )

    # Check ships -> accessions potential matches
    cursor.execute("SELECT COUNT(*) FROM ships WHERE accession IS NOT NULL")
    ships_with_accession = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM ships")
    total_ships = cursor.fetchone()[0]

    print(f"ships records with accession values: {ships_with_accession}/{total_ships}")

    return potential_accession_matches, total_joined_ships


def step1_add_accession_id_to_joined_ships(conn):
    """Step 1: Add accession_id column to joined_ships table"""
    print("\n=== STEP 1: ADD accession_id TO joined_ships ===")

    cursor = conn.cursor()

    # Check if column already exists
    cursor.execute("PRAGMA table_info(joined_ships)")
    columns = [col[1] for col in cursor.fetchall()]

    if "accession_id" not in columns:
        cursor.execute("ALTER TABLE joined_ships ADD COLUMN accession_id INTEGER")
        conn.commit()
        print("Added accession_id column to joined_ships table")
    else:
        print("accession_id column already exists in joined_ships table")


def step2_populate_accession_id(conn):
    """Step 2: Populate accession_id by matching starshipID to accessions"""
    print("\n=== STEP 2: POPULATE accession_id IN joined_ships ===")

    # Find matches between joined_ships.starshipID and accessions.ship_name
    query = """
    SELECT 
        js.rowid,
        js.starshipID,
        a.id as accession_id
    FROM joined_ships js
    JOIN accessions a ON js.starshipID = a.ship_name
    WHERE js.accession_id IS NULL
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches)} accession_id matches to populate")

    if len(matches) > 0:
        cursor = conn.cursor()
        updated_count = 0

        for _, row in matches.iterrows():
            cursor.execute(
                "UPDATE joined_ships SET accession_id = ? WHERE rowid = ?",
                (row["accession_id"], row["rowid"]),
            )
            updated_count += 1

            if updated_count % 1000 == 0:
                print(f"Updated {updated_count} accession_id relationships...")

        conn.commit()
        print(f"Successfully populated {updated_count} accession_id relationships")

        # Verify
        cursor.execute(
            "SELECT COUNT(*) FROM joined_ships WHERE accession_id IS NOT NULL"
        )
        populated_count = cursor.fetchone()[0]
        print(f"Total joined_ships records with accession_id: {populated_count}")

    else:
        print("No accession_id relationships found to populate")


def step3_fix_ships_accession_relationship(conn):
    """Step 3: Fix ships table to properly reference accessions"""
    print("\n=== STEP 3: FIX ships.accession_id RELATIONSHIP ===")

    cursor = conn.cursor()

    # Check if ships table has accession_id column
    cursor.execute("PRAGMA table_info(ships)")
    columns = [col[1] for col in cursor.fetchall()]

    if "accession_id" not in columns:
        cursor.execute("ALTER TABLE ships ADD COLUMN accession_id INTEGER")
        conn.commit()
        print("Added accession_id column to ships table")

    # Check current state
    cursor.execute("SELECT COUNT(*) FROM ships WHERE header IS NULL")
    ships_with_null_header = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM ships WHERE header IS NOT NULL")
    ships_with_header = cursor.fetchone()[0]

    print(f"Ships with NULL header: {ships_with_null_header}")
    print(f"Ships with non-NULL header: {ships_with_header}")

    # Strategy 1: For ships with NULL header, match by ID (ships.id = accessions.id)
    print("\nStrategy 1: Matching ships with NULL header by ID correspondence...")

    null_header_query = """
    SELECT 
        s.id as ship_id,
        a.id as accession_id
    FROM ships s
    JOIN accessions a ON s.id = a.id
    WHERE s.header IS NULL 
    AND s.accession_id IS NULL
    """

    null_header_matches = pd.read_sql_query(null_header_query, conn)
    print(
        f"Found {len(null_header_matches)} ships with NULL header that can be matched by ID"
    )

    if len(null_header_matches) > 0:
        for _, row in null_header_matches.iterrows():
            cursor.execute(
                "UPDATE ships SET accession_id = ? WHERE id = ?",
                (row["accession_id"], row["ship_id"]),
            )
        conn.commit()
        print(
            f"Updated {len(null_header_matches)} ships.accession_id relationships (NULL header case)"
        )

    # Strategy 2: For ships with non-NULL header, match via joined_ships intermediary
    print("\nStrategy 2: Matching ships with non-NULL header via joined_ships...")

    header_match_query = """
    SELECT 
        s.id as ship_id,
        s.header,
        a.id as accession_id
    FROM ships s
    JOIN joined_ships js ON s.header = js.starshipID
    JOIN accessions a ON js.starshipID = a.ship_name
    WHERE s.header IS NOT NULL 
    AND s.accession_id IS NULL
    """

    header_matches = pd.read_sql_query(header_match_query, conn)
    print(
        f"Found {len(header_matches)} ships with non-NULL header that can be matched via joined_ships"
    )

    if len(header_matches) > 0:
        for _, row in header_matches.iterrows():
            cursor.execute(
                "UPDATE ships SET accession_id = ? WHERE id = ?",
                (row["accession_id"], row["ship_id"]),
            )
        conn.commit()
        print(
            f"Updated {len(header_matches)} ships.accession_id relationships (header via joined_ships case)"
        )

    # Check for orphaned ships (non-NULL header but no match)
    cursor.execute("""
        SELECT COUNT(*) FROM ships s
        LEFT JOIN joined_ships js ON s.header = js.starshipID
        WHERE s.header IS NOT NULL 
        AND s.accession_id IS NULL
        AND js.starshipID IS NULL
    """)
    orphaned_ships = cursor.fetchone()[0]
    if orphaned_ships > 0:
        print(
            f"Warning: {orphaned_ships} ships with non-NULL headers cannot be matched to any accession"
        )

    total_updated = len(null_header_matches) + len(header_matches)
    print(f"\nTotal ships.accession_id relationships updated: {total_updated}")


def step4_update_joined_ships_ship_id(conn):
    """Step 4: Update joined_ships.ship_id to reference properly linked ships"""
    print("\n=== STEP 4: UPDATE joined_ships.ship_id RELATIONSHIPS ===")

    # Now that we have accession_id in joined_ships and ships,
    # we can fix the ship_id relationships
    query = """
    SELECT 
        js.rowid,
        js.accession_id,
        s.id as correct_ship_id,
        js.ship_id as current_ship_id
    FROM joined_ships js
    JOIN ships s ON js.accession_id = s.accession_id
    WHERE js.ship_id != s.id OR js.ship_id IS NULL
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches)} ship_id relationships to fix")

    if len(matches) > 0:
        cursor = conn.cursor()
        updated_count = 0

        for _, row in matches.iterrows():
            cursor.execute(
                "UPDATE joined_ships SET ship_id = ? WHERE rowid = ?",
                (row["correct_ship_id"], row["rowid"]),
            )
            updated_count += 1

            if updated_count % 1000 == 0:
                print(f"Updated {updated_count} ship_id relationships...")

        conn.commit()
        print(f"Successfully updated {updated_count} ship_id relationships")


def validate_migration(conn):
    """Validate that the migration worked correctly"""
    print("\n=== MIGRATION VALIDATION ===")

    cursor = conn.cursor()

    # Check joined_ships -> accessions relationship
    cursor.execute("""
        SELECT 
            COUNT(*) as total_joined_ships,
            COUNT(js.accession_id) as with_accession_id,
            COUNT(*) - COUNT(js.accession_id) as missing_accession_id
        FROM joined_ships js
    """)

    accession_result = cursor.fetchone()
    print(
        f"joined_ships → accessions: {accession_result[1]}/{accession_result[0]} complete ({accession_result[2]} missing)"
    )

    # Check joined_ships -> ships relationship via accessions
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js
        JOIN accessions a ON js.accession_id = a.id
        JOIN ships s ON js.ship_id = s.id AND s.accession_id = a.id
    """)

    valid_relationships = cursor.fetchone()[0]
    print(
        f"Valid joined_ships → accessions → ships relationships: {valid_relationships}"
    )

    # Check ships -> accessions relationship
    cursor.execute("""
        SELECT 
            COUNT(*) as total_ships,
            COUNT(s.accession_id) as with_accession_id,
            COUNT(*) - COUNT(s.accession_id) as missing_accession_id
        FROM ships s
    """)

    ships_result = cursor.fetchone()
    print(
        f"ships → accessions: {ships_result[1]}/{ships_result[0]} complete ({ships_result[2]} missing)"
    )


def main():
    """Main function to run the schema migration"""

    try:
        conn = connect_to_database()

        # Create backup
        backup_database(conn)

        # Analyze current state
        analyze_current_state(conn)

        # Run migration steps
        step1_add_accession_id_to_joined_ships(conn)
        step2_populate_accession_id(conn)
        step3_fix_ships_accession_relationship(conn)
        step4_update_joined_ships_ship_id(conn)

        # Validate migration
        validate_migration(conn)

        conn.close()
        print("\n=== SCHEMA MIGRATION COMPLETED ===")
        print(
            "The accessions table is now the primary table with proper foreign key relationships!"
        )

    except Exception as e:
        print(f"Error during migration: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

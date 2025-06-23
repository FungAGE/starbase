#!/usr/bin/env python3
"""
Script to fix relationships between joined_ships, ships, and captains tables
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
    backup_path = DB_PATHS["starbase"] + ".backup"

    if os.path.exists(backup_path):
        os.remove(backup_path)

    # Create backup
    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def fix_ship_id_relationships(conn):
    """Fix missing ship_id relationships in joined_ships table"""

    print("\n=== FIXING SHIP_ID RELATIONSHIPS ===")

    # Find records that can be matched by header
    query = """
    SELECT 
        js.rowid,
        js.starshipID,
        s.id as correct_ship_id
    FROM joined_ships js
    JOIN ships s ON js.starshipID = s.header
    WHERE js.ship_id IS NULL
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches)} ship_id matches to fix")

    if len(matches) > 0:
        # Update in batches to avoid memory issues
        cursor = conn.cursor()
        batch_size = 1000
        updated_count = 0

        for i in range(0, len(matches), batch_size):
            batch = matches.iloc[i : i + batch_size]

            for _, row in batch.iterrows():
                cursor.execute(
                    "UPDATE joined_ships SET ship_id = ? WHERE rowid = ?",
                    (row["correct_ship_id"], row["rowid"]),
                )
                updated_count += 1

                if updated_count % 100 == 0:
                    print(f"Updated {updated_count} ship_id relationships...")

        conn.commit()
        print(f"Successfully updated {updated_count} ship_id relationships")

        # Verify the fix
        cursor.execute("SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NULL")
        remaining_null = cursor.fetchone()[0]
        print(f"Remaining NULL ship_id values: {remaining_null}")

    else:
        print("No ship_id relationships found to fix")


def fix_captain_id_relationships(conn):
    """Fix missing captain_id relationships in joined_ships table"""

    print("\n=== FIXING CAPTAIN_ID RELATIONSHIPS ===")

    # Find records that can be matched by captainID
    query = """
    SELECT 
        js.rowid,
        js.captainID,
        c.id as correct_captain_id
    FROM joined_ships js
    JOIN captains c ON js.captainID = c.captainID
    WHERE js.captain_id IS NULL 
    AND js.captainID IS NOT NULL 
    AND js.captainID != ''
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches)} captain_id matches to fix")

    if len(matches) > 0:
        cursor = conn.cursor()
        batch_size = 1000
        updated_count = 0

        for i in range(0, len(matches), batch_size):
            batch = matches.iloc[i : i + batch_size]

            for _, row in batch.iterrows():
                cursor.execute(
                    "UPDATE joined_ships SET captain_id = ? WHERE rowid = ?",
                    (row["correct_captain_id"], row["rowid"]),
                )
                updated_count += 1

                if updated_count % 100 == 0:
                    print(f"Updated {updated_count} captain_id relationships...")

        conn.commit()
        print(f"Successfully updated {updated_count} captain_id relationships")

        # Verify the fix
        cursor.execute("""
            SELECT COUNT(*) FROM joined_ships 
            WHERE captain_id IS NULL 
            AND captainID IS NOT NULL 
            AND captainID != ''
        """)
        remaining_null = cursor.fetchone()[0]
        print(
            f"Remaining NULL captain_id values (with valid captainID): {remaining_null}"
        )

    else:
        print("No captain_id relationships found to fix")


def validate_fixes(conn):
    """Validate that the fixes worked correctly"""

    print("\n=== VALIDATION RESULTS ===")

    # Check ship_id relationships
    cursor = conn.cursor()

    cursor.execute("""
        SELECT 
            COUNT(*) as total_joined_ships,
            COUNT(js.ship_id) as with_ship_id,
            COUNT(*) - COUNT(js.ship_id) as missing_ship_id
        FROM joined_ships js
    """)

    ship_result = cursor.fetchone()
    print(
        f"Ship relationships: {ship_result[1]}/{ship_result[0]} complete ({ship_result[2]} missing)"
    )

    # Check captain relationships
    cursor.execute("""
        SELECT 
            COUNT(*) as total_with_captainID,
            COUNT(js.captain_id) as with_captainID_new,
            COUNT(*) - COUNT(js.captain_id) as missing_captainID_new
        FROM joined_ships js
        WHERE js.captainID IS NOT NULL AND js.captainID != ''
    """)

    captain_result = cursor.fetchone()
    print(
        f"Captain relationships: {captain_result[1]}/{captain_result[0]} complete ({captain_result[2]} missing)"
    )

    # Check for successful joins
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js
        JOIN ships s ON js.ship_id = s.id
        JOIN captains c ON js.captain_id = c.id
        WHERE js.captainID IS NOT NULL AND js.captainID != ''
    """)

    successful_joins = cursor.fetchone()[0]
    print(
        f"Complete relationships (joined_ships → ships → captains): {successful_joins}"
    )


def analyze_remaining_issues(conn):
    """Analyze what issues remain after fixes"""

    print("\n=== REMAINING ISSUES ANALYSIS ===")

    # Find records that still can't be matched
    cursor = conn.cursor()

    # Ships that can't be matched
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js
        LEFT JOIN ships s ON js.starshipID = s.header
        WHERE js.ship_id IS NULL AND s.id IS NULL
    """)

    unmatched_ships = cursor.fetchone()[0]
    print(f"Ships that cannot be matched by starshipID: {unmatched_ships}")

    if unmatched_ships > 0:
        print("Sample of unmatched ships:")
        unmatched_sample = pd.read_sql_query(
            """
            SELECT DISTINCT js.starshipID
            FROM joined_ships js
            LEFT JOIN ships s ON js.starshipID = s.header
            WHERE js.ship_id IS NULL AND s.id IS NULL
            LIMIT 10
        """,
            conn,
        )
        print(unmatched_sample.to_string(index=False))

    # Captains that can't be matched
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js
        LEFT JOIN captains c ON js.captainID = c.captainID
        WHERE js.captain_id IS NULL 
        AND js.captainID IS NOT NULL 
        AND js.captainID != ''
        AND c.id IS NULL
    """)

    unmatched_captains = cursor.fetchone()[0]
    print(f"\nCaptains that cannot be matched by captainID: {unmatched_captains}")

    if unmatched_captains > 0:
        print("Sample of unmatched captains:")
        unmatched_captain_sample = pd.read_sql_query(
            """
            SELECT DISTINCT js.captainID
            FROM joined_ships js
            LEFT JOIN captains c ON js.captainID = c.captainID
            WHERE js.captain_id IS NULL 
            AND js.captainID IS NOT NULL 
            AND js.captainID != ''
            AND c.id IS NULL
            LIMIT 10
        """,
            conn,
        )
        print(unmatched_captain_sample.to_string(index=False))


def main():
    """Main function to run all fixes"""

    try:
        conn = connect_to_database()

        # Create backup before making changes
        backup_database(conn)

        # Fix relationships
        fix_ship_id_relationships(conn)
        fix_captain_id_relationships(conn)

        # Validate fixes
        validate_fixes(conn)

        # Analyze remaining issues
        analyze_remaining_issues(conn)

        conn.close()
        print("\n=== FIXES COMPLETED ===")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

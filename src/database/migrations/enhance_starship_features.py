#!/usr/bin/env python3
"""
Script to enhance starship_features table with additional columns from joined_ships

Changes to implement:
1. Add columns: dr, tir, target, spok, ars, other, hgt
2. If sf.upTIR and sf.downTIR are null, parse js.tir format "<upTIR>/<downTIR>"
3. If sf.upDR and sf.downDR are null, parse js.dr format "<upDR>/<downDR>"
4. Coalesce js.size and sf.elementLength (keep sf.elementLength if not null)

Format notes:
- TIR/DR format: "value1/value2" where "./." means null
- Size: prefer sf.elementLength over js.size
"""

import sqlite3
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
    backup_path = DB_PATHS["starbase"] + ".enhance_features_backup"

    if os.path.exists(backup_path):
        os.remove(backup_path)

    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def analyze_current_state(conn):
    """Analyze current state before making changes"""
    print("\n=== CURRENT STATE ANALYSIS ===")

    cursor = conn.cursor()

    # Check what columns will be added
    new_columns = ["dr", "tir", "target", "spok", "ars", "other", "hgt"]

    print("1. Data availability in joined_ships for new columns:")
    for col in new_columns:
        cursor.execute(f"""
            SELECT COUNT(*) as total, 
                   COUNT(CASE WHEN {col} IS NOT NULL AND {col} != '' THEN 1 END) as with_data
            FROM joined_ships
        """)
        total, with_data = cursor.fetchone()
        pct = (with_data / total * 100) if total > 0 else 0
        print(f"  {col:8}: {with_data:,}/{total:,} have data ({pct:.1f}%)")

    # Check TIR/DR parsing opportunities
    print("\n2. TIR/DR parsing opportunities:")
    cursor.execute("""
        SELECT COUNT(*) as total_tir_data,
               COUNT(CASE WHEN sf.upTIR IS NULL AND sf.downTIR IS NULL THEN 1 END) as can_parse_tir
        FROM joined_ships js
        JOIN starship_features sf ON js.starshipID = sf.starshipID
        WHERE js.tir IS NOT NULL AND js.tir != '' AND js.tir != './.'
    """)
    tir_total, tir_parseable = cursor.fetchone()
    print(f"  TIR: {tir_parseable:,}/{tir_total:,} records can be parsed")

    cursor.execute("""
        SELECT COUNT(*) as total_dr_data,
               COUNT(CASE WHEN sf.upDR IS NULL AND sf.downDR IS NULL THEN 1 END) as can_parse_dr
        FROM joined_ships js
        JOIN starship_features sf ON js.starshipID = sf.starshipID
        WHERE js.dr IS NOT NULL AND js.dr != '' AND js.dr != './.'
    """)
    dr_total, dr_parseable = cursor.fetchone()
    print(f"  DR:  {dr_parseable:,}/{dr_total:,} records can be parsed")

    # Check size coalescing
    print("\n3. Size coalescing opportunities:")
    cursor.execute("""
        SELECT COUNT(CASE WHEN sf.elementLength IS NULL AND js.size IS NOT NULL THEN 1 END) as can_use_size
        FROM joined_ships js
        JOIN starship_features sf ON js.starshipID = sf.starshipID
    """)
    size_updates = cursor.fetchone()[0]
    print(f"  Size: {size_updates:,} records can use js.size for NULL elementLength")


def add_new_columns(conn):
    """Add new columns to starship_features table"""
    print("\n=== ADDING NEW COLUMNS ===")

    cursor = conn.cursor()

    # Get current schema
    cursor.execute("PRAGMA table_info(starship_features)")
    current_columns = [col[1] for col in cursor.fetchall()]

    # Columns to add from joined_ships
    new_columns = {
        "dr": "TEXT",
        "tir": "TEXT",
        "target": "TEXT",
        "spok": "TEXT",
        "ars": "TEXT",
        "other": "TEXT",
        "hgt": "TEXT",
    }

    # Add each column if it doesn't exist
    added_columns = []
    for col_name, col_type in new_columns.items():
        if col_name not in current_columns:
            try:
                cursor.execute(
                    f"ALTER TABLE starship_features ADD COLUMN {col_name} {col_type}"
                )
                added_columns.append(col_name)
                print(f"✅ Added column: {col_name} ({col_type})")
            except Exception as e:
                print(f"❌ Error adding column {col_name}: {e}")
        else:
            print(f"⚠️  Column {col_name} already exists")

    conn.commit()
    print(f"Added {len(added_columns)} new columns")
    return added_columns


def populate_new_columns(conn):
    """Populate new columns with data from joined_ships"""
    print("\n=== POPULATING NEW COLUMNS ===")

    cursor = conn.cursor()

    # Update all new columns with data from joined_ships where available
    new_columns = ["dr", "tir", "target", "spok", "ars", "other", "hgt"]

    total_updates = 0
    for col in new_columns:
        print(f"Updating {col}...")

        update_query = f"""
            UPDATE starship_features 
            SET {col} = (
                SELECT js.{col}
                FROM joined_ships js
                WHERE js.starshipID = starship_features.starshipID
                AND js.{col} IS NOT NULL 
                AND js.{col} != ''
                LIMIT 1
            )
            WHERE starship_features.{col} IS NULL
        """

        cursor.execute(update_query)
        updated = cursor.rowcount
        total_updates += updated
        print(f"  ✅ Updated {updated:,} records for {col}")

    conn.commit()
    print(f"Total column updates: {total_updates:,}")


def parse_tir_dr_values(conn):
    """Parse TIR and DR values from js.tir and js.dr format"""
    print("\n=== PARSING TIR AND DR VALUES ===")

    cursor = conn.cursor()

    # Parse TIR values where sf.upTIR and sf.downTIR are null
    print("1. Parsing TIR values...")

    # Get records that need TIR parsing
    cursor.execute("""
        SELECT sf.id, js.tir
        FROM starship_features sf
        JOIN joined_ships js ON sf.starshipID = js.starshipID
        WHERE sf.upTIR IS NULL 
        AND sf.downTIR IS NULL
        AND js.tir IS NOT NULL 
        AND js.tir != '' 
        AND js.tir != './.'
    """)

    tir_records = cursor.fetchall()
    print(f"Found {len(tir_records):,} records to parse TIR values")

    tir_updates = 0
    for record_id, tir_value in tir_records:
        if "/" in tir_value:
            up_tir, down_tir = tir_value.split("/", 1)

            # Handle null values represented as "."
            up_tir = up_tir if up_tir != "." else None
            down_tir = down_tir if down_tir != "." else None

            cursor.execute(
                """
                UPDATE starship_features 
                SET upTIR = ?, downTIR = ? 
                WHERE id = ?
            """,
                (up_tir, down_tir, record_id),
            )
            tir_updates += 1

    print(f"  ✅ Updated {tir_updates:,} TIR records")

    # Parse DR values where sf.upDR and sf.downDR are null
    print("2. Parsing DR values...")

    cursor.execute("""
        SELECT sf.id, js.dr
        FROM starship_features sf
        JOIN joined_ships js ON sf.starshipID = js.starshipID
        WHERE sf.upDR IS NULL 
        AND sf.downDR IS NULL
        AND js.dr IS NOT NULL 
        AND js.dr != '' 
        AND js.dr != './.'
    """)

    dr_records = cursor.fetchall()
    print(f"Found {len(dr_records):,} records to parse DR values")

    dr_updates = 0
    for record_id, dr_value in dr_records:
        if "/" in dr_value:
            up_dr, down_dr = dr_value.split("/", 1)

            # Handle null values represented as "."
            up_dr = up_dr if up_dr != "." else None
            down_dr = down_dr if down_dr != "." else None

            cursor.execute(
                """
                UPDATE starship_features 
                SET upDR = ?, downDR = ? 
                WHERE id = ?
            """,
                (up_dr, down_dr, record_id),
            )
            dr_updates += 1

    print(f"  ✅ Updated {dr_updates:,} DR records")

    conn.commit()
    print(f"Total TIR/DR parsing updates: {tir_updates + dr_updates:,}")


def coalesce_size_data(conn):
    """Coalesce js.size with sf.elementLength (prefer sf.elementLength)"""
    print("\n=== COALESCING SIZE DATA ===")

    cursor = conn.cursor()

    # Update elementLength with js.size where elementLength is null
    print("Updating elementLength with js.size where elementLength is NULL...")

    update_query = """
        UPDATE starship_features 
        SET elementLength = (
            SELECT js.size
            FROM joined_ships js
            WHERE js.starshipID = starship_features.starshipID
            AND js.size IS NOT NULL
            LIMIT 1
        )
        WHERE starship_features.elementLength IS NULL
    """

    cursor.execute(update_query)
    updated = cursor.rowcount

    print(f"✅ Updated {updated:,} records with size from joined_ships")

    conn.commit()


def validate_changes(conn):
    """Validate the changes made to starship_features"""
    print("\n=== VALIDATION ===")

    cursor = conn.cursor()

    # Check new columns data
    new_columns = ["dr", "tir", "target", "spok", "ars", "other", "hgt"]

    print("1. New columns data availability:")
    for col in new_columns:
        cursor.execute(f"""
            SELECT COUNT(*) as total,
                   COUNT(CASE WHEN {col} IS NOT NULL AND {col} != '' THEN 1 END) as with_data
            FROM starship_features
        """)
        total, with_data = cursor.fetchone()
        pct = (with_data / total * 100) if total > 0 else 0
        print(f"  {col:8}: {with_data:,}/{total:,} have data ({pct:.1f}%)")

    # Check TIR/DR parsing results
    print("\n2. TIR/DR parsing results:")
    cursor.execute("""
        SELECT COUNT(CASE WHEN upTIR IS NOT NULL OR downTIR IS NOT NULL THEN 1 END) as with_tir_data,
               COUNT(CASE WHEN upDR IS NOT NULL OR downDR IS NOT NULL THEN 1 END) as with_dr_data
        FROM starship_features
    """)
    tir_data, dr_data = cursor.fetchone()
    print(f"  Records with TIR data: {tir_data:,}")
    print(f"  Records with DR data: {dr_data:,}")

    # Check elementLength coverage
    print("\n3. Element length coverage:")
    cursor.execute("""
        SELECT COUNT(*) as total,
               COUNT(CASE WHEN elementLength IS NOT NULL THEN 1 END) as with_length
        FROM starship_features
    """)
    total, with_length = cursor.fetchone()
    pct = (with_length / total * 100) if total > 0 else 0
    print(f"  Records with elementLength: {with_length:,}/{total:,} ({pct:.1f}%)")

    # Sample of enhanced data
    print("\n4. Sample of enhanced data:")
    cursor.execute("""
        SELECT starshipID, elementLength, upTIR, downTIR, upDR, downDR, target, spok
        FROM starship_features 
        WHERE (upTIR IS NOT NULL OR upDR IS NOT NULL OR target IS NOT NULL)
        LIMIT 5
    """)

    sample_data = cursor.fetchall()
    for row in sample_data:
        starship, length, up_tir, down_tir, up_dr, down_dr, target, spok = row
        print(
            f"  {starship[:25]:25}: len={length} upTIR={up_tir} downTIR={down_tir} upDR={up_dr} downDR={down_dr} target={target} spok={spok}"
        )


def main():
    """Main function to enhance starship_features table"""

    try:
        conn = connect_to_database()

        # Create backup
        backup_database(conn)

        # Analyze current state
        analyze_current_state(conn)

        # Add new columns
        added_columns = add_new_columns(conn)

        # Populate new columns with data from joined_ships
        populate_new_columns(conn)

        # Parse TIR and DR values
        parse_tir_dr_values(conn)

        # Coalesce size data
        coalesce_size_data(conn)

        # Validate changes
        validate_changes(conn)

        conn.close()

        print("\n=== COMPLETED ===")
        print("✅ Successfully enhanced starship_features table")
        print(f"✅ Added {len(added_columns)} new columns")
        print("✅ Parsed TIR/DR values from joined_ships format")
        print("✅ Coalesced size data with elementLength")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

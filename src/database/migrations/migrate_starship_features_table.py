#!/usr/bin/env python3
"""
Migration script to restructure the starship_features table

Current Issues:
- starship_features.ship_id currently points to what should be accession_id
- No proper link to ships table
- Each accession can have multiple features (duplicated accession_id is OK)

Proposed Changes:
1. Rename current ship_id → accession_id
2. Add new ship_id column that properly links to ships table
3. Populate new ship_id by finding the appropriate ship for each starshipID

Current schema pattern shows:
- starship_features.starshipID matches accessions.ship_name
- starship_features.ship_id matches accessions.id (should become accession_id)
- Need to find ships.id that corresponds to each starshipID
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
    backup_path = DB_PATHS["starbase"] + ".starship_features_migration_backup"

    if os.path.exists(backup_path):
        os.remove(backup_path)

    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def analyze_current_state(conn):
    """Analyze the current state of the starship_features table"""
    print("\n=== STARSHIP_FEATURES PRE-MIGRATION ANALYSIS ===")

    cursor = conn.cursor()

    # Get total records
    cursor.execute("SELECT COUNT(*) FROM starship_features")
    total_features = cursor.fetchone()[0]
    print(f"Total starship_features records: {total_features:,}")

    # Check relationships with joined_ships (our source of truth)
    cursor.execute("""
        SELECT COUNT(*) FROM starship_features sf
        JOIN joined_ships js ON sf.starshipID = js.starshipID
    """)
    starshipid_matches = cursor.fetchone()[0]
    print(f"starshipID → joined_ships.starshipID matches: {starshipid_matches:,}")

    # Check how many joined_ships records have the corrected relationships
    cursor.execute("""
        SELECT COUNT(*) FROM starship_features sf
        JOIN joined_ships js ON sf.starshipID = js.starshipID
        WHERE js.accession_id IS NOT NULL
    """)
    js_with_accession = cursor.fetchone()[0]
    print(f"Matches with valid joined_ships.accession_id: {js_with_accession:,}")

    cursor.execute("""
        SELECT COUNT(*) FROM starship_features sf
        JOIN joined_ships js ON sf.starshipID = js.starshipID
        WHERE js.ship_id IS NOT NULL
    """)
    js_with_ship = cursor.fetchone()[0]
    print(f"Matches with valid joined_ships.ship_id: {js_with_ship:,}")

    # Check how many unique accessions/ships are involved
    cursor.execute("SELECT COUNT(DISTINCT starshipID) FROM starship_features")
    unique_starships = cursor.fetchone()[0]
    print(f"Unique starshipID values: {unique_starships:,}")

    cursor.execute("SELECT COUNT(DISTINCT ship_id) FROM starship_features")
    unique_current_ship_ids = cursor.fetchone()[0]
    print(f"Unique current ship_id values: {unique_current_ship_ids:,}")

    # Check how many features per starship (to confirm duplicates are expected)
    cursor.execute("""
        SELECT AVG(feature_count) as avg_features_per_starship
        FROM (
            SELECT starshipID, COUNT(*) as feature_count
            FROM starship_features 
            GROUP BY starshipID
        )
    """)
    avg_features = cursor.fetchone()[0]
    print(f"Average features per starshipID: {avg_features:.2f}")

    return total_features, starshipid_matches


def step1_add_accession_id_column(conn):
    """Step 1: Add accession_id column to starship_features table"""
    print("\n=== STEP 1: ADD accession_id COLUMN ===")

    cursor = conn.cursor()

    # Check if column already exists
    cursor.execute("PRAGMA table_info(starship_features)")
    columns = [col[1] for col in cursor.fetchall()]

    if "accession_id" not in columns:
        cursor.execute("ALTER TABLE starship_features ADD COLUMN accession_id INTEGER")
        conn.commit()
        print("Added accession_id column to starship_features table")
    else:
        print("accession_id column already exists in starship_features table")


def step2_populate_accession_id(conn):
    """Step 2: Populate accession_id by inheriting from joined_ships"""
    print("\n=== STEP 2: POPULATE accession_id FROM joined_ships ===")

    # Use joined_ships as source of truth for accession_id
    # Since starship_features.ship_id currently contains what should be accession_id,
    # and we've confirmed this matches joined_ships.accession_id
    query = """
    SELECT 
        sf.rowid,
        sf.starshipID,
        js.accession_id
    FROM starship_features sf
    JOIN joined_ships js ON sf.starshipID = js.starshipID
    WHERE sf.accession_id IS NULL
    AND js.accession_id IS NOT NULL
    """

    matches = pd.read_sql_query(query, conn)
    print(f"Found {len(matches):,} accession_id matches to populate from joined_ships")

    if len(matches) > 0:
        cursor = conn.cursor()
        batch_size = 1000
        updated_count = 0

        # Group by starshipID and accession_id to handle duplicates efficiently
        unique_matches = matches.drop_duplicates(subset=["rowid"])

        for i in range(0, len(unique_matches), batch_size):
            batch = unique_matches.iloc[i : i + batch_size]

            for _, row in batch.iterrows():
                cursor.execute(
                    "UPDATE starship_features SET accession_id = ? WHERE rowid = ?",
                    (row["accession_id"], row["rowid"]),
                )
                updated_count += 1

                if updated_count % 1000 == 0:
                    print(f"Updated {updated_count:,} accession_id relationships...")

        conn.commit()
        print(f"Successfully populated {updated_count:,} accession_id relationships")

        # Verify the update
        cursor.execute(
            "SELECT COUNT(*) FROM starship_features WHERE accession_id IS NOT NULL"
        )
        populated_count = cursor.fetchone()[0]
        print(f"Total starship_features records with accession_id: {populated_count:,}")

    else:
        print("No accession_id relationships found to populate")


def step3_add_new_ship_id_column(conn):
    """Step 3: Add new ship_id column that will properly link to ships table"""
    print("\n=== STEP 3: ADD NEW ship_id COLUMN ===")

    cursor = conn.cursor()

    # Check current table structure
    cursor.execute("PRAGMA table_info(starship_features)")
    columns = [col[1] for col in cursor.fetchall()]

    if "new_ship_id" not in columns:
        cursor.execute("ALTER TABLE starship_features ADD COLUMN new_ship_id INTEGER")
        conn.commit()
        print("Added new_ship_id column to starship_features table")
    else:
        print("new_ship_id column already exists in starship_features table")


def step4_populate_new_ship_id(conn):
    """Step 4: Populate new ship_id by inheriting from joined_ships and finding unique ships"""
    print("\n=== STEP 4: POPULATE NEW ship_id FROM joined_ships ===")

    cursor = conn.cursor()

    # Strategy 1: Use joined_ships.ship_id where available
    print("Strategy 1: Inheriting ship_id from joined_ships...")
    query1 = """
    SELECT 
        sf.rowid,
        sf.starshipID,
        js.ship_id
    FROM starship_features sf
    JOIN joined_ships js ON sf.starshipID = js.starshipID
    WHERE sf.new_ship_id IS NULL
    AND js.ship_id IS NOT NULL
    """

    matches1 = pd.read_sql_query(query1, conn)
    print(f"Found {len(matches1):,} ship_id matches from joined_ships")

    # Strategy 2: For records without ship_id in joined_ships, find ships via accession
    print("Strategy 2: Finding ships via accession relationship...")
    query2 = """
    SELECT 
        sf.rowid,
        sf.starshipID,
        s.id as ship_id
    FROM starship_features sf
    JOIN joined_ships js ON sf.starshipID = js.starshipID
    JOIN ships s ON js.accession_id = s.accession_id
    WHERE sf.new_ship_id IS NULL
    AND js.ship_id IS NULL
    AND js.accession_id IS NOT NULL
    AND s.accession_id IS NOT NULL
    """

    matches2 = pd.read_sql_query(query2, conn)
    print(f"Found {len(matches2):,} ship_id matches via accession")

    # Strategy 3: Direct header matching as fallback
    print("Strategy 3: Direct header matching as fallback...")
    query3 = """
    SELECT 
        sf.rowid,
        sf.starshipID,
        s.id as ship_id
    FROM starship_features sf
    JOIN ships s ON sf.starshipID = s.header
    WHERE sf.new_ship_id IS NULL
    AND sf.rowid NOT IN (
        SELECT sf2.rowid FROM starship_features sf2
        JOIN joined_ships js2 ON sf2.starshipID = js2.starshipID
        WHERE js2.ship_id IS NOT NULL OR js2.accession_id IS NOT NULL
    )
    """

    matches3 = pd.read_sql_query(query3, conn)
    print(f"Found {len(matches3):,} ship_id matches via direct header")

    # Combine all strategies, prioritizing joined_ships data
    all_matches = pd.concat([matches1, matches2, matches3], ignore_index=True)

    # Remove duplicates, keeping first occurrence (prioritizes Strategy 1)
    all_matches = all_matches.drop_duplicates(subset=["rowid"], keep="first")
    print(f"Total unique new ship_id matches to populate: {len(all_matches):,}")

    if len(all_matches) > 0:
        updated_count = 0

        for _, row in all_matches.iterrows():
            cursor.execute(
                "UPDATE starship_features SET new_ship_id = ? WHERE rowid = ?",
                (row["ship_id"], row["rowid"]),
            )
            updated_count += 1

            if updated_count % 1000 == 0:
                print(f"Updated {updated_count:,} new ship_id relationships...")

        conn.commit()
        print(f"Successfully populated {updated_count:,} new ship_id relationships")

        # Verify the update
        cursor.execute(
            "SELECT COUNT(*) FROM starship_features WHERE new_ship_id IS NOT NULL"
        )
        populated_count = cursor.fetchone()[0]
        print(f"Total starship_features records with new_ship_id: {populated_count:,}")

        # Check for unique ships
        cursor.execute(
            "SELECT COUNT(DISTINCT new_ship_id) FROM starship_features WHERE new_ship_id IS NOT NULL"
        )
        unique_ships = cursor.fetchone()[0]
        print(f"Unique ships referenced in starship_features: {unique_ships:,}")

    else:
        print("No new ship_id relationships found to populate")


def step5_restructure_table(conn):
    """Step 5: Restructure table to rename columns properly"""
    print("\n=== STEP 5: RESTRUCTURE TABLE ===")

    cursor = conn.cursor()

    # Create new table with the correct structure
    print("Creating new table with proper column names...")

    cursor.execute("""
        CREATE TABLE starship_features_new (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession_id INTEGER,
            ship_id INTEGER,
            contigID VARCHAR,
            starshipID VARCHAR,
            captainID VARCHAR,
            elementBegin VARCHAR,
            elementEnd VARCHAR,
            elementLength VARCHAR,
            strand VARCHAR,
            boundaryType VARCHAR,
            emptySiteID VARCHAR,
            emptyContig VARCHAR,
            emptyBegin VARCHAR,
            emptyEnd VARCHAR,
            emptySeq VARCHAR,
            upDR VARCHAR,
            downDR VARCHAR,
            DRedit VARCHAR,
            upTIR VARCHAR,
            downTIR VARCHAR,
            TIRedit VARCHAR,
            nestedInside VARCHAR,
            containNested VARCHAR,
            FOREIGN KEY (accession_id) REFERENCES accessions(id),
            FOREIGN KEY (ship_id) REFERENCES ships(id)
        )
    """)

    # Copy data from old table to new table
    print("Copying data to new table structure...")
    cursor.execute("""
        INSERT INTO starship_features_new (
            accession_id, ship_id, contigID, starshipID, captainID, 
            elementBegin, elementEnd, elementLength, strand, boundaryType,
            emptySiteID, emptyContig, emptyBegin, emptyEnd, emptySeq,
            upDR, downDR, DRedit, upTIR, downTIR, TIRedit, 
            nestedInside, containNested
        )
        SELECT 
            accession_id, new_ship_id, contigID, starshipID, captainID,
            elementBegin, elementEnd, elementLength, strand, boundaryType,
            emptySiteID, emptyContig, emptyBegin, emptyEnd, emptySeq,
            upDR, downDR, DRedit, upTIR, downTIR, TIRedit,
            nestedInside, containNested
        FROM starship_features
    """)

    # Drop old table and rename new one
    cursor.execute("DROP TABLE starship_features")
    cursor.execute("ALTER TABLE starship_features_new RENAME TO starship_features")

    # Create indexes for performance
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_starship_features_accession_id ON starship_features(accession_id)"
    )
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_starship_features_ship_id ON starship_features(ship_id)"
    )
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_starship_features_starshipID ON starship_features(starshipID)"
    )

    conn.commit()
    print("✅ Successfully restructured starship_features table")

    # Verify the new structure
    cursor.execute("PRAGMA table_info(starship_features)")
    new_schema = cursor.fetchall()
    print("\nNew starship_features table structure:")
    for col in new_schema:
        pk_status = "PRIMARY KEY" if col[5] == 1 else ""
        print(f"  {col[1]} ({col[2]}) {pk_status}")


def validate_migration(conn):
    """Validate that the migration worked correctly"""
    print("\n=== MIGRATION VALIDATION ===")

    cursor = conn.cursor()

    # Check starship_features -> accessions relationship
    cursor.execute("""
        SELECT 
            COUNT(*) as total_features,
            COUNT(sf.accession_id) as with_accession_id,
            COUNT(*) - COUNT(sf.accession_id) as missing_accession_id
        FROM starship_features sf
    """)

    accession_result = cursor.fetchone()
    print(
        f"starship_features → accessions: {accession_result[1]:,}/{accession_result[0]:,} complete ({accession_result[2]:,} missing)"
    )

    # Check starship_features -> ships relationship
    cursor.execute("""
        SELECT 
            COUNT(*) as total_features,
            COUNT(sf.ship_id) as with_ship_id,
            COUNT(*) - COUNT(sf.ship_id) as missing_ship_id
        FROM starship_features sf
    """)

    ship_result = cursor.fetchone()
    print(
        f"starship_features → ships: {ship_result[1]:,}/{ship_result[0]:,} complete ({ship_result[2]:,} missing)"
    )

    # Check for valid relationships
    cursor.execute("""
        SELECT COUNT(*) FROM starship_features sf
        JOIN accessions a ON sf.accession_id = a.id
    """)
    valid_accession_relationships = cursor.fetchone()[0]
    print(f"Valid accession_id relationships: {valid_accession_relationships:,}")

    cursor.execute("""
        SELECT COUNT(*) FROM starship_features sf
        JOIN ships s ON sf.ship_id = s.id
    """)
    valid_ship_relationships = cursor.fetchone()[0]
    print(f"Valid ship_id relationships: {valid_ship_relationships:,}")

    # Summary statistics
    cursor.execute(
        "SELECT COUNT(DISTINCT accession_id) FROM starship_features WHERE accession_id IS NOT NULL"
    )
    unique_accessions = cursor.fetchone()[0]
    print(f"Unique accessions with features: {unique_accessions:,}")

    cursor.execute(
        "SELECT COUNT(DISTINCT ship_id) FROM starship_features WHERE ship_id IS NOT NULL"
    )
    unique_ships = cursor.fetchone()[0]
    print(f"Unique ships with features: {unique_ships:,}")


def main():
    """Main function to run the starship_features table migration"""

    try:
        conn = connect_to_database()

        # Create backup
        backup_database(conn)

        # Analyze current state
        analyze_current_state(conn)

        # Run migration steps
        step1_add_accession_id_column(conn)
        step2_populate_accession_id(conn)
        step3_add_new_ship_id_column(conn)
        step4_populate_new_ship_id(conn)
        step5_restructure_table(conn)

        # Validate migration
        validate_migration(conn)

        conn.close()
        print("\n=== STARSHIP_FEATURES MIGRATION COMPLETED ===")
        print("✅ Added accession_id column linking to accessions table")
        print("✅ Updated ship_id column to properly link to ships table")
        print("✅ Maintained support for multiple features per accession")
        print("✅ Added proper primary key and foreign key constraints")

    except Exception as e:
        print(f"Error during migration: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

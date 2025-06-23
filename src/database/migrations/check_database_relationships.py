#!/usr/bin/env python3
"""
Script to check relationships and completeness between joined_ships, ships, and captains tables
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


def analyze_table_completeness(conn):
    """Analyze the completeness of joined_ships, ships, and captains tables"""

    print("=" * 80)
    print("DATABASE RELATIONSHIP ANALYSIS")
    print("=" * 80)

    # Get basic counts
    cursor = conn.cursor()

    cursor.execute("SELECT COUNT(*) FROM joined_ships")
    joined_ships_count = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM ships")
    ships_count = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM captains")
    captains_count = cursor.fetchone()[0]

    print("Table Counts:")
    print(f"  joined_ships: {joined_ships_count:,}")
    print(f"  ships: {ships_count:,}")
    print(f"  captains: {captains_count:,}")
    print()

    return joined_ships_count, ships_count, captains_count


def check_joined_ships_relationships(conn):
    """Check relationships from joined_ships to other tables"""

    print("JOINED_SHIPS → SHIPS RELATIONSHIP ANALYSIS")
    print("-" * 50)

    # Check ship_id relationships
    query = """
    SELECT 
        COUNT(*) as total_joined_ships,
        COUNT(DISTINCT js.ship_id) as unique_ship_ids_in_joined,
        COUNT(DISTINCT s.id) as matching_ships,
        COUNT(*) - COUNT(s.id) as missing_ships
    FROM joined_ships js
    LEFT JOIN ships s ON js.ship_id = s.id
    """

    result = pd.read_sql_query(query, conn)
    print(result.to_string(index=False))
    print()

    # Find missing ship_ids
    missing_ships_query = """
    SELECT DISTINCT js.ship_id, js.starshipID
    FROM joined_ships js
    LEFT JOIN ships s ON js.ship_id = s.id
    WHERE s.id IS NULL
    LIMIT 10
    """

    missing_ships = pd.read_sql_query(missing_ships_query, conn)
    if not missing_ships.empty:
        print("Sample of missing ship_ids in ships table:")
        print(missing_ships.to_string(index=False))
        print()

    print("JOINED_SHIPS → CAPTAINS RELATIONSHIP ANALYSIS")
    print("-" * 50)

    # Check captain_id relationships
    query = """
    SELECT 
        COUNT(*) as total_joined_ships,
        COUNT(DISTINCT js.captain_id) as unique_captain_ids_in_joined,
        COUNT(DISTINCT c.id) as matching_captains,
        COUNT(*) - COUNT(c.id) as missing_captains
    FROM joined_ships js
    LEFT JOIN captains c ON js.captain_id = c.id
    WHERE js.captain_id IS NOT NULL AND js.captain_id != ''
    """

    result = pd.read_sql_query(query, conn)
    print(result.to_string(index=False))
    print()

    # Find missing captain_ids
    missing_captains_query = """
    SELECT DISTINCT js.captain_id, js.captainID, js.starshipID
    FROM joined_ships js
    LEFT JOIN captains c ON js.captain_id = c.id
    WHERE js.captain_id IS NOT NULL 
    AND js.captain_id != ''
    AND c.id IS NULL
    LIMIT 10
    """

    missing_captains = pd.read_sql_query(missing_captains_query, conn)
    if not missing_captains.empty:
        print("Sample of missing captain_id in captains table:")
        print(missing_captains.to_string(index=False))
        print()


def check_ships_captains_relationship(conn):
    """Check relationship between ships and captains tables"""

    print("SHIPS ↔ CAPTAINS RELATIONSHIP ANALYSIS")
    print("-" * 50)

    # Check ship_id in captains table
    query = """
    SELECT 
        COUNT(*) as total_captains,
        COUNT(DISTINCT c.ship_id) as unique_ship_ids_in_captains,
        COUNT(DISTINCT s.id) as matching_ships,
        COUNT(*) - COUNT(s.id) as captains_missing_ships
    FROM captains c
    LEFT JOIN ships s ON c.ship_id = s.id
    WHERE c.ship_id IS NOT NULL
    """

    result = pd.read_sql_query(query, conn)
    print(result.to_string(index=False))
    print()

    # Reverse check - ships without captains
    ships_without_captains_query = """
    SELECT 
        COUNT(*) as total_ships,
        COUNT(DISTINCT c.ship_id) as ships_with_captains,
        COUNT(*) - COUNT(DISTINCT c.ship_id) as ships_without_captains
    FROM ships s
    LEFT JOIN captains c ON s.id = c.ship_id
    """

    result = pd.read_sql_query(ships_without_captains_query, conn)
    print("Ships without captains:")
    print(result.to_string(index=False))
    print()


def analyze_data_quality(conn):
    """Analyze data quality issues"""

    print("DATA QUALITY ANALYSIS")
    print("-" * 50)

    # Check for NULL/empty values in key fields
    quality_checks = [
        ("joined_ships", "starshipID", "NULL or empty starshipID"),
        ("joined_ships", "captainID", "NULL or empty captainID"),
        ("joined_ships", "ship_id", "NULL ship_id"),
        ("ships", "sequence", "NULL sequences"),
        ("captains", "sequence", "NULL captain sequences"),
        ("captains", "captainID", "NULL or empty captainID"),
    ]

    for table, column, description in quality_checks:
        if column == "ship_id":
            query = f"SELECT COUNT(*) FROM {table} WHERE {column} IS NULL"
        else:
            query = (
                f"SELECT COUNT(*) FROM {table} WHERE {column} IS NULL OR {column} = ''"
            )

        cursor = conn.cursor()
        cursor.execute(query)
        count = cursor.fetchone()[0]
        print(f"{description}: {count}")

    print()


def suggest_fixes(conn):
    """Suggest potential fixes for missing relationships"""

    print("SUGGESTED FIXES")
    print("-" * 50)

    # 1. Check if we can match ships by accession/starshipID
    print("1. Potential ship_id matches by starshipID:")

    match_query = """
    SELECT 
        js.starshipID,
        js.ship_id as current_ship_id,
        s.id as potential_ship_id,
        s.accession
    FROM joined_ships js
    LEFT JOIN ships s_current ON js.ship_id = s_current.id
    LEFT JOIN ships s ON js.starshipID = s.accession OR js.starshipID = s.header
    WHERE s_current.id IS NULL 
    AND s.id IS NOT NULL
    LIMIT 10
    """

    matches = pd.read_sql_query(match_query, conn)
    if not matches.empty:
        print(matches.to_string(index=False))
        print(f"Found {len(matches)} potential matches (showing first 10)")
    else:
        print("No obvious matches found by starshipID")
    print()

    # 2. Check if we can match captains by captainID
    print("2. Potential captain matches by captainID:")

    captain_match_query = """
    SELECT 
        js.captainID,
        js.captain_id as current_captain_id,
        c.id as potential_captain_id,
        c.captainID
    FROM joined_ships js
    LEFT JOIN captains c_current ON js.captain_id = c_current.id
    LEFT JOIN captains c ON js.captainID = c.captainID
    WHERE js.captain_id IS NOT NULL 
    AND c_current.id IS NULL 
    AND c.id IS NOT NULL
    LIMIT 10
    """

    captain_matches = pd.read_sql_query(captain_match_query, conn)
    if not captain_matches.empty:
        print(captain_matches.to_string(index=False))
        print(
            f"Found {len(captain_matches)} potential captain matches (showing first 10)"
        )
    else:
        print("No obvious captain matches found by captainID")
    print()


def generate_update_queries(conn):
    """Generate SQL update queries to fix relationships"""

    print("GENERATED UPDATE QUERIES")
    print("-" * 50)

    # Generate queries to fix ship_id relationships
    print("-- Fix ship_id relationships based on starshipID matches:")

    update_ship_query = """
    SELECT 
        'UPDATE joined_ships SET ship_id = ' || s.id || 
        ' WHERE starshipID = ''' || js.starshipID || ''' AND ship_id = ' || js.ship_id || ';' as update_query
    FROM joined_ships js
    LEFT JOIN ships s_current ON js.ship_id = s_current.id
    LEFT JOIN ships s ON js.starshipID = s.accession OR js.starshipID = s.header
    WHERE s_current.id IS NULL 
    AND s.id IS NOT NULL
    LIMIT 5
    """

    update_queries = pd.read_sql_query(update_ship_query, conn)
    if not update_queries.empty:
        for query in update_queries["update_query"]:
            print(query)
    print()

    # Generate queries to fix captain relationships
    print("-- Fix captain_id relationships based on captainID matches:")

    update_captain_query = """
    SELECT 
        'UPDATE joined_ships SET captain_id = ' || c.id || 
        ' WHERE captainID = ''' || js.captainID || ''' AND captain_id = ' || COALESCE(js.captain_id, 'NULL') || ';' as update_query
    FROM joined_ships js
    LEFT JOIN captains c_current ON js.captain_id = c_current.id
    LEFT JOIN captains c ON js.captainID = c.captainID
    WHERE (js.captain_id IS NULL OR c_current.id IS NULL)
    AND c.id IS NOT NULL
    LIMIT 5
    """

    captain_update_queries = pd.read_sql_query(update_captain_query, conn)
    if not captain_update_queries.empty:
        for query in captain_update_queries["update_query"]:
            print(query)
    print()


def main():
    """Main function to run all analyses"""

    try:
        conn = connect_to_database()

        # Run all analyses
        analyze_table_completeness(conn)
        check_joined_ships_relationships(conn)
        check_ships_captains_relationship(conn)
        analyze_data_quality(conn)
        suggest_fixes(conn)
        generate_update_queries(conn)

        conn.close()

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

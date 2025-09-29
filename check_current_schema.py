#!/usr/bin/env python3
"""
Check the current database schema to understand the existing structure
"""

import sqlite3
import sys

def check_schema(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Get joined_ships table structure
        print("=== JOINED_SHIPS TABLE STRUCTURE ===")
        cursor.execute("PRAGMA table_info(joined_ships)")
        columns = cursor.fetchall()
        for col in columns:
            print(f"  {col[1]} {col[2]} {'NOT NULL' if col[3] else 'NULL'} {'PK' if col[5] else ''}")
        
        print("\n=== JOINED_SHIPS FOREIGN KEYS ===")
        cursor.execute("PRAGMA foreign_key_list(joined_ships)")
        fks = cursor.fetchall()
        for fk in fks:
            print(f"  {fk[3]} -> {fk[2]}.{fk[4]}")
        
        print("\n=== JOINED_SHIPS INDEXES ===")
        cursor.execute("SELECT name, sql FROM sqlite_master WHERE type='index' AND tbl_name='joined_ships'")
        indexes = cursor.fetchall()
        for idx in indexes:
            print(f"  {idx[0]}: {idx[1]}")
        
        print("\n=== SAMPLE DATA ===")
        cursor.execute("SELECT id, starshipID, accession_id FROM joined_ships LIMIT 5")
        rows = cursor.fetchall()
        print("id | starshipID | accession_id")
        print("-" * 30)
        for row in rows:
            print(f"{row[0]} | {row[1]} | {row[2] if row[2] is not None else 'NULL'}")
        
        print(f"\n=== COUNTS ===")
        cursor.execute("SELECT COUNT(*) FROM joined_ships")
        total = cursor.fetchone()[0]
        print(f"Total joined_ships: {total}")
        
        cursor.execute("SELECT COUNT(*) FROM joined_ships WHERE accession_id IS NOT NULL")
        with_accession_id = cursor.fetchone()[0]
        print(f"With accession_id: {with_accession_id}")
        
        cursor.execute("SELECT COUNT(*) FROM joined_ships WHERE accession_id IS NULL")
        without_accession_id = cursor.fetchone()[0]
        print(f"Without accession_id: {without_accession_id}")
        
    finally:
        conn.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_current_schema.py <database_path>")
        sys.exit(1)
    
    check_schema(sys.argv[1])

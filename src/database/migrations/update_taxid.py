#!/usr/bin/env python3
"""
TaxID Update Script
==================

This script updates the taxID field in the taxonomy table by querying NCBI's taxonomy database
using esearch for records where taxID is NULL but name is not NULL.

Requirements:
- NCBI E-utilities (esearch, esummary, xtract) must be installed
- Internet connection for NCBI queries
- Optional: NCBI API key for better rate limits

Usage:
- Set NCBI_API_KEY environment variable or modify the API_KEY variable below
- python3 update_taxid.py
"""

import sqlite3
import subprocess
import sys
import time
import os
from pathlib import Path

# NCBI API Key - Set your API key here or use environment variable
API_KEY = os.getenv("NCBI_API_KEY", "")


def connect_to_database(db_path):
    """Connect to SQLite database and return connection."""
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error as e:
        print(f"Error connecting to database: {e}")
        sys.exit(1)


def get_records_needing_taxid(conn):
    """Get all records where taxID is NULL but name is not NULL."""
    cursor = conn.cursor()
    cursor.execute("""
        SELECT id, name 
        FROM taxonomy 
        WHERE taxID IS NULL 
        AND name IS NOT NULL 
        AND TRIM(name) != ''
        ORDER BY id
    """)
    return cursor.fetchall()


def parse_name_components(name):
    """Parse a taxonomic name into genus, species components."""
    parts = name.strip().split()

    if len(parts) >= 2:
        genus = parts[0]
        species = parts[1]
        genus_species = f"{genus} {species}"
        return genus, species, genus_species
    elif len(parts) == 1:
        genus = parts[0]
        return genus, None, genus
    else:
        return None, None, None


def search_ncbi_taxid_single(query, api_key=None):
    """Search NCBI taxonomy database for a single query using shell command."""
    try:
        # Don't pass API key as parameter - it should be set as environment variable
        # The NCBI tools automatically use NCBI_API_KEY when it's set in environment
        cmd = f'esearch -query "{query}" -db taxonomy | esummary | xtract -pattern DocumentSummary -element TaxId'

        # Execute the command
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, timeout=30
        )

        if result.returncode == 0:
            output = result.stdout.strip()
            if output:
                # Handle multiple results - take the first valid TaxID
                lines = [line.strip() for line in output.split("\n") if line.strip()]
                for line in lines:
                    if line.isdigit():
                        return line
        else:
            # Print error for debugging
            if result.stderr:
                print(f"    Command error: {result.stderr.strip()}")

        return None

    except subprocess.TimeoutExpired:
        print(f"    Timeout while searching for '{query}'")
        return None
    except Exception as e:
        print(f"    Exception while searching for '{query}': {e}")
        return None


def search_ncbi_taxid(name, api=None):
    """Search NCBI taxonomy database with fallback strategy."""
    # First try the full name
    print(f"    Trying full name: '{name}'")
    taxid = search_ncbi_taxid_single(name, API_KEY)
    if taxid:
        print(f"    ✓ Found with full name: {taxid}")
        return taxid

    # Parse name components
    genus, species, genus_species = parse_name_components(name)

    if not genus:
        print("    ❌ Could not parse genus from name")
        return None

    # If full name didn't work and we have both genus and species, try genus + species
    if species and genus_species != name:
        print(f"    Trying genus + species: '{genus_species}'")
        taxid = search_ncbi_taxid_single(genus_species, API_KEY)
        if taxid:
            print(f"    ✓ Found with genus + species: {taxid}")
            return taxid

    # Finally, try just the genus
    if genus != name and genus != genus_species:
        print(f"    Trying genus only: '{genus}'")
        taxid = search_ncbi_taxid_single(genus, API_KEY)
        if taxid:
            print(f"    ✓ Found with genus only: {taxid}")
            return taxid

    print("    ❌ No taxID found with any search strategy")
    return None


def update_taxid_in_database(conn, record_id, taxid):
    """Update the taxID for a specific record."""
    cursor = conn.cursor()
    cursor.execute(
        """
        UPDATE taxonomy 
        SET taxID = ? 
        WHERE id = ?
    """,
        (taxid, record_id),
    )
    conn.commit()


def main():
    """Main function to update taxID fields."""
    db_path = "src/database/db/starbase.sqlite"

    if not Path(db_path).exists():
        print(f"Database file not found: {db_path}")
        sys.exit(1)

    # Check if NCBI tools are available
    try:
        subprocess.run(["esearch", "-help"], capture_output=True, check=True)
        print("✓ NCBI esearch tool found")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ NCBI esearch tool not found. Please install NCBI E-utilities first.")
        print("Visit: https://www.ncbi.nlm.nih.gov/books/NBK179288/")
        sys.exit(1)

    # Check API key
    if API_KEY:
        print(f"✓ Using NCBI API key: {API_KEY[:8]}...")
    else:
        print(
            "⚠️  No NCBI API key found. Consider setting NCBI_API_KEY environment variable"
        )
        print("   for better rate limits (10 requests/second vs 3 requests/second)")

    print("\nStarting taxID update process...")
    print("=" * 50)

    # Connect to database
    conn = connect_to_database(db_path)

    try:
        # Get records needing updates
        records = get_records_needing_taxid(conn)
        total_records = len(records)

        print(f"Found {total_records} records needing taxID updates")

        if total_records == 0:
            print("No updates needed!")
            return

        print("\nStarting NCBI searches with fallback strategy...")
        print("Strategy: Full name → Genus + Species → Genus only")
        print("This may take a while due to NCBI rate limiting...")

        successful_updates = 0
        failed_searches = 0

        # Adjust delay based on API key availability
        delay = 0.1 if API_KEY else 0.5  # 10 req/sec with API key, 2 req/sec without

        for i, record in enumerate(records, 1):
            record_id = record["id"]
            name = record["name"]

            print(f"\n[{i}/{total_records}] Searching for: {name}")

            # Search NCBI with fallback strategy
            taxid = search_ncbi_taxid(name, API_KEY)

            if taxid:
                # Update database
                update_taxid_in_database(conn, record_id, taxid)
                print(f"  ✅ Updated record ID {record_id} with taxID: {taxid}")
                successful_updates += 1
            else:
                print(f"  ❌ No taxID found for: {name}")
                failed_searches += 1

            # Rate limiting: wait between requests
            if i < total_records:  # Don't wait after the last request
                time.sleep(delay)

        print("\n" + "=" * 60)
        print("Update Summary:")
        print(f"  Total records processed: {total_records}")
        print(f"  Successful updates: {successful_updates}")
        print(f"  Failed searches: {failed_searches}")
        print(f"  Success rate: {(successful_updates / total_records) * 100:.1f}%")

        if successful_updates > 0:
            print(f"\n✅ Successfully updated {successful_updates} taxID fields!")

        if failed_searches > 0:
            print(f"\n⚠️  {failed_searches} records could not be found in NCBI taxonomy")
            print("This may be due to:")
            print("  - Strain-specific names not in NCBI")
            print("  - Variations in naming conventions")
            print("  - Outdated or deprecated names")
            print("  - Names that don't exist in NCBI taxonomy database")

    except KeyboardInterrupt:
        print("\n\nProcess interrupted by user")
    except Exception as e:
        print(f"\nError during update process: {e}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()

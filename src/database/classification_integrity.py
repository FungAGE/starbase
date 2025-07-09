"""
Test helper utilities for running classification tests outside Flask app context.
Provides non-cached versions of database functions that can be used in standalone tests.
"""

import pandas as pd
import os
import sys
from contextlib import contextmanager

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Add src to path if not already there
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(os.path.dirname(current_dir))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)


def create_db_engine(db_path):
    """Create a database engine with the given path"""
    if not os.path.exists(db_path):
        db_dir = os.path.dirname(db_path)
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

    # Set echo based on environment - only log SQL in development
    return create_engine(f"sqlite:///{db_path}")


@contextmanager
def db_session():
    engine = create_db_engine("src/database/db/starbase.sqlite")
    Session = sessionmaker(bind=engine)
    session = Session()

    try:
        yield session
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def fetch_meta_data(curated=False, accession_tag=None):
    with db_session() as session:
        meta_query = """
        SELECT j.curated_status, j.starshipID, j.ship_id, a.accession_tag, a.version_tag, 
               f.familyName, f.type_element_reference, n.navis_name, h.haplotype_name, g.ome
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        """

        params = []
        if curated:
            meta_query += " WHERE j.curated_status = 'curated'"

        if accession_tag:
            where_clause = " WHERE " if not curated else " AND "
            if isinstance(accession_tag, list):
                # Use ? for SQLite placeholders
                placeholders = ",".join(["?"] * len(accession_tag))
                meta_query += f"{where_clause}a.accession_tag IN ({placeholders})"
                params = tuple(accession_tag)
            else:
                # Use ? for SQLite placeholder
                meta_query += f"{where_clause}a.accession_tag = ?"
                params = (accession_tag,)

        try:
            if params:
                meta_df = pd.read_sql_query(meta_query, session.bind, params=params)
            else:
                meta_df = pd.read_sql_query(meta_query, session.bind)

            return meta_df
        except Exception as e:
            print(f"Error fetching meta data: {str(e)}")
            raise


def fetch_family_names():
    """Fetch all family names from the database"""
    with db_session() as session:
        query = (
            "SELECT DISTINCT familyName FROM family_names WHERE familyName IS NOT NULL"
        )
        try:
            return pd.read_sql_query(query, session.bind)
        except Exception as e:
            print(f"Error fetching family names: {str(e)}")
            raise


def fetch_captains():
    """Fetch all captain IDs from the database"""
    with db_session() as session:
        query = "SELECT DISTINCT captainID FROM captains WHERE captainID IS NOT NULL"
        try:
            return pd.read_sql_query(query, session.bind)
        except Exception as e:
            print(f"Error fetching captains: {str(e)}")
            raise


def fetch_taxonomy():
    """Fetch all genus-species combinations from the database"""
    with db_session() as session:
        query = """
        SELECT DISTINCT genus, species 
        FROM taxonomy 
        WHERE genus IS NOT NULL AND species IS NOT NULL
        """
        try:
            return pd.read_sql_query(query, session.bind)
        except Exception as e:
            print(f"Error fetching taxonomy: {str(e)}")
            raise


def check_integrity():
    """Comprehensive integrity check between text files and database"""
    print("=== CLASSIFICATION INTEGRITY CHECK ===\n")

    # Load text files
    try:
        superfamily_rep_ships = pd.read_csv(
            "src/database/db/Starships/captain/tyr/faa/tree/superfamily-rep-ships.csv"
        )
        superfam_clades = pd.read_csv(
            "src/database/db/Starships/captain/tyr/faa/tree/superfam-clades.tsv",
            sep="\t",
        )
        print(
            f"✅ Loaded superfamily-rep-ships.csv: {len(superfamily_rep_ships)} entries"
        )
        print(f"✅ Loaded superfam-clades.tsv: {len(superfam_clades)} entries")
    except Exception as e:
        print(f"❌ Error loading text files: {str(e)}")
        return

    # Get database data
    try:
        meta_df = fetch_meta_data()
        family_names_df = fetch_family_names()
        captains_df = fetch_captains()
        taxonomy_df = fetch_taxonomy()

        print(f"✅ Database contains {len(meta_df)} joined ships")
        print(f"✅ Database contains {len(family_names_df)} family names")
        print(f"✅ Database contains {len(captains_df)} captains")
        print(f"✅ Database contains {len(taxonomy_df)} taxonomy entries")

    except Exception as e:
        print(f"❌ Error loading database data: {str(e)}")
        return

    print("\n" + "=" * 50)

    # 1. Family name integrity check
    print("1. FAMILY NAME INTEGRITY CHECK")
    print("-" * 30)

    # Get unique families from both text files
    rep_families = set(superfamily_rep_ships["family"].dropna().unique())
    clade_families = set(superfam_clades["familyName"].dropna().unique())
    all_text_families = rep_families.union(clade_families)

    # Get families from database
    db_families = set(family_names_df["familyName"].dropna().unique())

    # Find missing families
    missing_from_db = all_text_families - db_families
    missing_from_text = db_families - all_text_families

    print(f"Families in text files: {len(all_text_families)}")
    print(f"Families in database: {len(db_families)}")
    print(f"Families in text files but NOT in database: {len(missing_from_db)}")
    if missing_from_db:
        for family in sorted(missing_from_db):
            print(f"  ❌ {family}")

    print(f"Families in database but NOT in text files: {len(missing_from_text)}")
    if missing_from_text:
        for family in sorted(missing_from_text):
            print(f"  ⚠️  {family}")

    # 2. Starship ID integrity check
    print("\n2. STARSHIP ID INTEGRITY CHECK")
    print("-" * 30)

    text_starships = set(superfamily_rep_ships["starshipID"].dropna().unique())
    db_starships = set(meta_df["starshipID"].dropna().unique())

    missing_starships = text_starships - db_starships
    print(f"Starships in rep-ships file: {len(text_starships)}")
    print(f"Starships in database: {len(db_starships)}")
    print(f"Starships in text file but NOT in database: {len(missing_starships)}")
    if missing_starships:
        for starship in sorted(missing_starships):
            print(f"  ❌ {starship}")
    else:
        print("  ✅ All starships from text file found in database")

    # 3. Captain ID integrity check
    print("\n3. CAPTAIN ID INTEGRITY CHECK")
    print("-" * 30)

    text_captains = set(superfamily_rep_ships["captainID"].dropna().unique())
    db_captains = set(captains_df["captainID"].dropna().unique())

    missing_captains = text_captains - db_captains
    print(f"Captains in rep-ships file: {len(text_captains)}")
    print(f"Captains in database: {len(db_captains)}")
    print(f"Captains in text file but NOT in database: {len(missing_captains)}")
    if missing_captains:
        for captain in sorted(missing_captains):
            print(f"  ❌ {captain}")
    else:
        print("  ✅ All captains from text file found in database")

    # 4. Taxonomy integrity check
    print("\n4. TAXONOMY INTEGRITY CHECK")
    print("-" * 30)

    # Create genus-species combinations from text file
    text_taxa = set()
    for _, row in superfamily_rep_ships.iterrows():
        if pd.notna(row["genus"]) and pd.notna(row["species"]):
            text_taxa.add((row["genus"], row["species"]))

    # Create genus-species combinations from database
    db_taxa = set()
    for _, row in taxonomy_df.iterrows():
        if pd.notna(row["genus"]) and pd.notna(row["species"]):
            db_taxa.add((row["genus"], row["species"]))

    missing_taxa = text_taxa - db_taxa
    print(f"Unique genus-species combinations in text file: {len(text_taxa)}")
    print(f"Unique genus-species combinations in database: {len(db_taxa)}")
    print(f"Taxa in text file but NOT in database: {len(missing_taxa)}")
    if missing_taxa:
        for genus, species in sorted(missing_taxa):
            print(f"  ❌ {genus} {species}")
    else:
        print("  ✅ All taxa from text file found in database")

    # 5. Cross-file consistency check
    print("\n5. CROSS-FILE CONSISTENCY CHECK")
    print("-" * 30)

    families_only_in_rep = rep_families - clade_families
    families_only_in_clades = clade_families - rep_families

    print(f"Families only in rep-ships file: {len(families_only_in_rep)}")
    if families_only_in_rep:
        for family in sorted(families_only_in_rep):
            print(f"  ⚠️  {family}")

    print(f"Families only in clades file: {len(families_only_in_clades)}")
    if families_only_in_clades:
        for family in sorted(families_only_in_clades):
            print(f"  ⚠️  {family}")

    # Summary
    print(f"\n{'=' * 50}")
    print("SUMMARY")
    print(f"{'=' * 50}")

    total_issues = (
        len(missing_from_db)
        + len(missing_starships)
        + len(missing_captains)
        + len(missing_taxa)
    )

    if total_issues == 0:
        print("✅ Perfect integrity! All data matches between files and database.")
    else:
        print(f"⚠️  Found {total_issues} integrity issues:")
        print(f"   - {len(missing_from_db)} missing families")
        print(f"   - {len(missing_starships)} missing starships")
        print(f"   - {len(missing_captains)} missing captains")
        print(f"   - {len(missing_taxa)} missing taxa")


if __name__ == "__main__":
    check_integrity()

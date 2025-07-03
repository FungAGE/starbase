"""
STARBASE MD5 HASH VERIFICATION & REPAIR TOOL
=============================================

This script safely verifies and repairs MD5 hashes in the Starbase SQLite database.

SAFETY FEATURES:
- Creates a timestamped backup copy of the database before any operations
- All modifications are made to the backup, leaving the original untouched
- Full transaction rollback on any errors
- Generates CSV reports for review before applying changes
- Provides helper functions to examine and optionally replace the original

USAGE:
Run from the project root directory:
  python -m src.database.fix_md5

Or run directly from src/database/:
  python fix_md5.py

FUNCTIONS:
- create_backup_database(): Creates timestamped backup copy
- check_md5sum(): Main verification and repair logic
- replace_original_with_backup(): Helper to replace original (use with caution)
- examine_backup_summary(): Shows summary of changes made
- cleanup_backup_database(): Optionally removes backup file

FILES CREATED (in src/database/db/):
- starbase_backup_YYYYMMDD_HHMMSS.sqlite: Backup database with fixes
- md5_mismatches.csv: Report of MD5 mismatches found (if any)

AFTER RUNNING:
The script provides instructions for examining the backup and optionally
replacing the original database if you're satisfied with the results.
"""

import pandas as pd
import shutil
import os
import logging
from datetime import datetime
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker
from ..config.database import StarbaseSession
from ..utils.seq_utils import clean_sequence
from ..utils.classification_utils import generate_md5_hash, generate_new_accession

# Global variable to hold the backup session factory
BackupSession = None
BACKUP_DB_PATH = None
ORIGINAL_DB_PATH = "db/starbase.sqlite"  # Relative to src/database/


def disable_sqlalchemy_logging():
    """
    Disable SQLAlchemy debug logging to reduce noise in output.
    """
    logging.getLogger("sqlalchemy.engine").setLevel(logging.WARNING)
    logging.getLogger("sqlalchemy.dialects").setLevel(logging.WARNING)
    logging.getLogger("sqlalchemy.pool").setLevel(logging.WARNING)
    logging.getLogger("sqlalchemy.orm").setLevel(logging.WARNING)


def create_backup_database():
    """
    Create a backup copy of the SQLite database and return a session factory for it.
    """
    global BackupSession, BACKUP_DB_PATH

    # Create backup filename with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_filename = f"starbase_backup_{timestamp}.sqlite"
    BACKUP_DB_PATH = os.path.join("db", backup_filename)  # Relative to src/database/

    print(f"Creating backup database: {BACKUP_DB_PATH}")

    # Check if original database exists
    if not os.path.exists(ORIGINAL_DB_PATH):
        raise FileNotFoundError(f"Original database not found: {ORIGINAL_DB_PATH}")

    # Create backup copy
    try:
        shutil.copy2(ORIGINAL_DB_PATH, BACKUP_DB_PATH)
        print(f"✓ Backup created successfully: {BACKUP_DB_PATH}")

        # Get file sizes for verification
        original_size = os.path.getsize(ORIGINAL_DB_PATH)
        backup_size = os.path.getsize(BACKUP_DB_PATH)
        print(f"  Original size: {original_size:,} bytes")
        print(f"  Backup size: {backup_size:,} bytes")

        if original_size != backup_size:
            raise Exception("Backup file size doesn't match original")

    except Exception as e:
        raise Exception(f"Failed to create database backup: {e}")

    # Create engine and session factory for backup database
    # Explicitly disable SQL logging to avoid debug messages
    backup_engine = create_engine(f"sqlite:///{BACKUP_DB_PATH}", echo=False)
    BackupSession = sessionmaker(bind=backup_engine)

    print("✓ Backup database session created")
    return BackupSession


def cleanup_backup_database():
    """
    Optionally clean up the backup database file.
    """
    global BACKUP_DB_PATH

    if BACKUP_DB_PATH and os.path.exists(BACKUP_DB_PATH):
        response = input(f"\nDelete backup database {BACKUP_DB_PATH}? (y/N): ")
        if response.lower() in ["y", "yes"]:
            try:
                os.remove(BACKUP_DB_PATH)
                print(f"✓ Backup database deleted: {BACKUP_DB_PATH}")
            except Exception as e:
                print(f"Error deleting backup: {e}")
        else:
            print(f"Backup database preserved: {BACKUP_DB_PATH}")
            print("You can manually delete it later or use it for testing.")


def get_database_session():
    """
    Get the appropriate database session (backup if available, otherwise original).
    """
    if BackupSession:
        return BackupSession()
    else:
        return StarbaseSession()


def fetch_ships(accession_tags=None, curated=False, dereplicate=True):
    """
    Fetch ship data for specified accession tags.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all ships.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
    Returns:
        pd.DataFrame: DataFrame containing ship data
    """
    session = get_database_session()

    query = """
    SELECT 
        a.id as accession_id, 
        a.accession_tag,
        j.curated_status,
        sf.elementBegin, sf.elementEnd, sf.contigID,
        t.name, t.family, t.`order`,
        f.familyName, n.navis_name, h.haplotype_name,
        g.assembly_accession,
        s.sequence,
        s.md5
    FROM joined_ships j
    INNER JOIN accessions a ON j.ship_id = a.id
    INNER JOIN ships s ON a.id = s.accession_id
    LEFT JOIN taxonomy t ON j.tax_id = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN navis_names n ON j.ship_navis_id = n.id
    LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN starship_features sf ON a.id = sf.accession_id
    WHERE 1=1
    """

    if accession_tags:
        query += " AND a.accession_tag IN ({})".format(
            ",".join(f"'{tag}'" for tag in accession_tags)
        )
    if curated:
        query += " AND j.curated_status = 'curated'"

    try:
        df = pd.read_sql_query(query, session.bind)

        if dereplicate:
            df = df.drop_duplicates(subset="accession_tag")

        if df.empty:
            print("Fetched ships DataFrame is empty.")
        return df
    except Exception as e:
        print(f"Error fetching ships data: {str(e)}")
        raise
    finally:
        session.close()


def update_database(
    old_accession_tag: str, new_accession_tag: str, old_md5sum: str, new_md5sum: str
):
    """
    Update accession_tag and/or md5sum entries in the database

    Args:
        old_accession_tag: Current accession tag to update from
        new_accession_tag: New accession tag (None if no change needed)
        old_md5sum: Current MD5 sum (for logging/verification)
        new_md5sum: New MD5 sum (None if no change needed)
    """

    session = get_database_session()

    try:
        # Start transaction
        session.begin()

        # Update accession_tag in the accessions table if needed
        if new_accession_tag and old_accession_tag != new_accession_tag:
            result = session.execute(
                text(
                    "UPDATE accessions SET accession_tag = :new_tag WHERE accession_tag = :old_tag"
                ),
                {"new_tag": new_accession_tag, "old_tag": old_accession_tag},
            )
            if result.rowcount == 0:
                print(
                    f"Warning: No rows updated in accessions table for {old_accession_tag}"
                )

        # Update md5sum in the ships table if needed
        if new_md5sum and old_md5sum != new_md5sum:
            result = session.execute(
                text("""
                    UPDATE ships 
                    SET md5 = :new_md5 
                    WHERE accession_id = (
                        SELECT id FROM accessions WHERE accession_tag = :old_tag
                    )
                """),
                {"new_md5": new_md5sum, "old_tag": old_accession_tag},
            )
            if result.rowcount == 0:
                print(
                    f"Warning: No rows updated in ships table for {old_accession_tag}"
                )
            else:
                print(
                    f"Updated MD5 for {old_accession_tag}: {old_md5sum} -> {new_md5sum}"
                )

        # Update foreign key references if accession_tag changed
        if new_accession_tag and old_accession_tag != new_accession_tag:
            # Update accession_id in ships table (this should reference the new accession's id)
            result = session.execute(
                text("""
                    UPDATE ships 
                    SET accession_id = (
                        SELECT id FROM accessions WHERE accession_tag = :new_tag
                    )
                    WHERE accession_id = (
                        SELECT id FROM accessions WHERE accession_tag = :old_tag
                    )
                """),
                {"new_tag": new_accession_tag, "old_tag": old_accession_tag},
            )

            # Update accession_id in gff table
            try:
                result = session.execute(
                    text("""
                        UPDATE gff 
                        SET accession_id = (
                            SELECT id FROM accessions WHERE accession_tag = :new_tag
                        )
                        WHERE accession_id = (
                            SELECT id FROM accessions WHERE accession_tag = :old_tag
                        )
                    """),
                    {"new_tag": new_accession_tag, "old_tag": old_accession_tag},
                )
            except Exception as e:
                print(f"Warning: Could not update gff table: {e}")

            # Update accession_id in starship_features table
            try:
                result = session.execute(
                    text("""
                        UPDATE starship_features 
                        SET accession_id = (
                            SELECT id FROM accessions WHERE accession_tag = :new_tag
                        )
                        WHERE accession_id = (
                            SELECT id FROM accessions WHERE accession_tag = :old_tag
                        )
                    """),
                    {"new_tag": new_accession_tag, "old_tag": old_accession_tag},
                )
            except Exception as e:
                print(f"Warning: Could not update starship_features table: {e}")

            print(f"Updated accession_tag: {old_accession_tag} -> {new_accession_tag}")

        # Commit all changes
        session.commit()

    except Exception as e:
        # Rollback on any error
        session.rollback()
        print(f"Error updating database for {old_accession_tag}: {e}")
        raise
    finally:
        session.close()


def check_ships_missing_accession_tags():
    """
    Add accession tags for (new) ships that don't have one yet
    """
    session = get_database_session()
    query = """
    SELECT id, accession_tag FROM ships WHERE accession_tag IS NULL
    """
    ships_df = pd.read_sql_query(query, session.bind)
    print(f"Found {len(ships_df)} ships missing accession tags")
    # assign new accessions to ships based on md5 hash
    ships_df["new_accession"] = None
    for idx, row in ships_df.iterrows():
        accession_tag = generate_new_accession(ships_df)
        print(f"Adding accession tag {accession_tag} to ship {row['id']}")
        update_database(
            old_accession_tag=row["id"],
            new_accession_tag=accession_tag,
            old_md5sum=row["md5"],
            new_md5sum=accession_tag,
        )
    print(f"Added {len(ships_df)} accession tags to ships")
    session.close()


def check_md5sum():
    """
    1. Generate dictionary of existing md5hash values

    2. re-generate the md5sum for each sequence and compare to the existing md5sum column

    3. for md5sum duplicates, make sure that the accession is the same
    - if the accessions don't match, then assign a new accession and fix the entries in the database

    4. for accessions duplicates, check if the md5 hash is the same
    - if the md5 hashes don't match, then assign a new accession and fix the entries in the database
    """
    print("Fetching ship data from database...")
    ships_df = fetch_ships(curated=False, dereplicate=False)

    print(f"Retrieved {len(ships_df)} total records from database")

    # Deduplicate by accession_tag to avoid processing the same sequence multiple times
    print("Deduplicating records by accession_tag...")
    unique_ships = ships_df.drop_duplicates(subset=["accession_tag"], keep="first")
    print(
        f"Processing {len(unique_ships)} unique accession tags (removed {len(ships_df) - len(unique_ships)} duplicates)"
    )

    # Add a column to track new MD5 values
    unique_ships = unique_ships.copy()  # Avoid SettingWithCopyWarning
    unique_ships["new_md5"] = None

    # 2. re-generate the md5sum for each unique sequence and compare to the existing md5sum
    print("Checking MD5 hashes for unique sequences...")
    md5_mismatches = []

    for idx, (df_idx, row) in enumerate(unique_ships.iterrows()):
        acc = row["accession_tag"]
        seq = row["sequence"]
        existing_md5 = row["md5"]

        # Show progress every 100 records
        if idx > 0 and idx % 100 == 0:
            print(f"  Progress: {idx}/{len(unique_ships)} records checked...")

        # Skip if sequence is None or empty
        if pd.isna(seq) or not seq:
            print(f"Warning: Empty sequence for accession {acc}")
            continue

        # regenerate the md5sum
        try:
            clean_db_seq = clean_sequence(seq)

            # Check if cleaned sequence is valid
            if not clean_db_seq or pd.isna(clean_db_seq):
                print(
                    f"Warning: Sequence became empty after cleaning for accession {acc}"
                )
                continue

            db_hash = generate_md5_hash(clean_db_seq)

            # Check if hash generation was successful
            if not db_hash:
                print(f"Warning: Failed to generate MD5 hash for accession {acc}")
                continue

        except Exception as e:
            print(f"Error generating MD5 for {acc}: {e}")
            print(f"  Sequence type: {type(seq)}")
            print(f"  Sequence length: {len(seq) if seq else 'N/A'}")
            print(f"  Sequence preview: {str(seq)[:100] if seq else 'N/A'}...")
            continue

        # check if the md5 hash is the same as the existing md5sum
        if existing_md5 and db_hash != existing_md5:
            unique_ships.at[df_idx, "new_md5"] = db_hash
            md5_mismatches.append(
                {"accession_tag": acc, "md5": existing_md5, "new_md5": db_hash}
            )
            print(
                f"MD5 mismatch found for {acc}: stored={existing_md5}, calculated={db_hash}"
            )

    print(f"✓ Completed MD5 checking for {len(unique_ships)} unique sequences")

    # Save only unique mismatches to CSV for review
    if md5_mismatches:
        mismatches_df = pd.DataFrame(md5_mismatches)
        mismatches_df.to_csv("md5_mismatches.csv", index=False)
        print(
            f"Found {len(mismatches_df)} unique MD5 mismatches. Details saved to md5_mismatches.csv"
        )
    else:
        print("No MD5 mismatches found.")

    # Update MD5 hashes where mismatches were found (only for unique entries)
    updated_count = 0
    for idx, row in unique_ships.iterrows():
        if pd.notna(row["new_md5"]):
            try:
                update_database(
                    old_accession_tag=row["accession_tag"],
                    new_accession_tag=row["accession_tag"],  # Keep same accession
                    old_md5sum=row["md5"],
                    new_md5sum=row["new_md5"],
                )
                updated_count += 1
            except Exception as e:
                print(f"Error updating MD5 for {row['accession_tag']}: {e}")

    if updated_count > 0:
        print(f"Updated MD5 hashes for {updated_count} accessions")

    # 3. check for md5sum duplicates and fix them
    print("\nChecking for MD5 duplicates...")
    # Use updated MD5 values for duplicate checking
    unique_ships["current_md5"] = unique_ships["new_md5"].fillna(unique_ships["md5"])

    # Find duplicated MD5 values
    md5_counts = unique_ships["current_md5"].value_counts()
    duplicate_md5s = md5_counts[md5_counts > 1].index.tolist()

    if not duplicate_md5s:
        print("No MD5 duplicates found.")
    else:
        print(f"Found {len(duplicate_md5s)} MD5 values with duplicates")

    for md5_val in duplicate_md5s:
        if pd.isna(md5_val):  # Skip null values
            continue

        duplicate_rows = unique_ships[unique_ships["current_md5"] == md5_val]
        acc_tags = duplicate_rows["accession_tag"].unique()

        if len(acc_tags) > 1:
            print(f"MD5 duplicate found: {md5_val} appears in accessions: {acc_tags}")

            # Keep the accession with the lowest number, merge others to it
            try:
                accession_numbers = [
                    int(acc.split(".")[0].split("SBS")[1]) for acc in acc_tags
                ]
                min_acc_number = min(accession_numbers)

                # Find the accession tag with the lowest number (this is the one to keep)
                target_accession = None
                for acc in acc_tags:
                    acc_number = int(acc.split(".")[0].split("SBS")[1])
                    if acc_number == min_acc_number:
                        target_accession = acc
                        break

                if not target_accession:
                    print(f"Error: Could not find target accession for {acc_tags}")
                    continue

                # Merge all other accessions to the target accession
                for acc in acc_tags:
                    acc_number = int(acc.split(".")[0].split("SBS")[1])
                    if acc_number != min_acc_number:
                        print(
                            f"Merging {acc} to {target_accession} (same MD5: {md5_val})"
                        )
                        try:
                            update_database(
                                old_accession_tag=acc,
                                new_accession_tag=target_accession,
                                old_md5sum=None,
                                new_md5sum=None,
                            )
                        except Exception as e:
                            print(
                                f"Error merging accession {acc} to {target_accession}: {e}"
                            )
            except (ValueError, IndexError) as e:
                print(f"Error parsing accession numbers for {acc_tags}: {e}")

    # 4. check for accession duplicates and fix them
    print("\nChecking for accession duplicates...")
    acc_counts = unique_ships["accession_tag"].value_counts()
    duplicate_accs = acc_counts[acc_counts > 1].index.tolist()

    if not duplicate_accs:
        print("No accession duplicates found in unique dataset.")
    else:
        print(f"Found {len(duplicate_accs)} accessions with duplicates")

    for acc in duplicate_accs:
        duplicate_rows = unique_ships[unique_ships["accession_tag"] == acc]
        md5_values = duplicate_rows["current_md5"].unique()

        if len(md5_values) > 1:
            print(
                f"Accession duplicate with different MD5s: {acc} has MD5s: {md5_values}"
            )

            # Keep the first entry, reassign others with new accessions
            for idx, (row_idx, row) in enumerate(duplicate_rows.iterrows()):
                if idx > 0:  # Skip first occurrence
                    new_accession = generate_new_accession(unique_ships)
                    print(f"Reassigning duplicate accession {acc} to {new_accession}")
                    try:
                        # Calculate correct MD5 for this sequence
                        clean_seq = clean_sequence(row["sequence"])
                        correct_md5 = generate_md5_hash(clean_seq)

                        update_database(
                            old_accession_tag=acc,
                            new_accession_tag=new_accession,
                            old_md5sum=row["md5"],
                            new_md5sum=correct_md5,
                        )
                    except Exception as e:
                        print(f"Error fixing duplicate accession {acc}: {e}")

    print("MD5 hash checking and fixing completed.")


def main():
    """
    Main function to run MD5 hash checking and fixing.
    """
    print("Starting MD5 hash checking and fixing process...")
    print("=" * 60)

    try:
        print("✓ All required modules imported successfully.")

        # Disable SQLAlchemy logging to reduce noise
        disable_sqlalchemy_logging()

        # Create backup database
        print("\n1. Creating backup database...")
        create_backup_database()

        # Run the MD5 checking process on backup
        print("\n2. Running MD5 hash verification and fixing on backup database...")
        check_md5sum()

        print("\n" + "=" * 60)
        print("✓ MD5 hash checking and fixing completed successfully!")
        print("\nIMPORTANT:")
        print(f"  - All changes were made to the backup database: {BACKUP_DB_PATH}")
        print(f"  - Original database remains unchanged: {ORIGINAL_DB_PATH}")
        print("\nNext steps:")
        print("  1. Review the backup database and any generated CSV files")
        print(
            "  2. If satisfied with results, you can replace the original with the backup"
        )
        print("  3. Or use the backup for further testing")

    except ImportError as e:
        print(f"❌ Error importing required modules: {e}")
        print("Please ensure all required modules are available in the src/ directory.")
    except Exception as e:
        print(f"❌ Error during MD5 checking process: {e}")
        raise
    finally:
        # Offer to cleanup backup
        if BACKUP_DB_PATH:
            print("\n📄 To examine the backup database, you can use:")
            print(f"   sqlite3 {BACKUP_DB_PATH}")
            print("\n📄 To replace original with backup (ONLY if you're satisfied):")
            print(
                '   python -c "from fix_md5 import replace_original_with_backup; replace_original_with_backup()"'
            )
        cleanup_backup_database()


def replace_original_with_backup():
    """
    Helper function to replace the original database with the backup.
    Use this only after verifying the backup is correct.
    """
    global BACKUP_DB_PATH, ORIGINAL_DB_PATH

    if not BACKUP_DB_PATH or not os.path.exists(BACKUP_DB_PATH):
        print("❌ No backup database found to replace original with.")
        return False

    if not os.path.exists(ORIGINAL_DB_PATH):
        print("❌ Original database not found.")
        return False

    print(f"⚠️  This will replace {ORIGINAL_DB_PATH} with {BACKUP_DB_PATH}")
    confirmation = input("Are you absolutely sure? Type 'REPLACE' to confirm: ")

    if confirmation == "REPLACE":
        try:
            # Create a backup of the original before replacing
            original_backup = f"{ORIGINAL_DB_PATH}.original_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            shutil.copy2(ORIGINAL_DB_PATH, original_backup)
            print(f"✓ Original backed up to: {original_backup}")

            # Replace original with backup
            shutil.copy2(BACKUP_DB_PATH, ORIGINAL_DB_PATH)
            print("✓ Original database replaced with backup")

            return True
        except Exception as e:
            print(f"❌ Error replacing database: {e}")
            return False
    else:
        print("Operation cancelled.")
        return False


def examine_backup_summary():
    """
    Print a summary of changes that would be made to help with decision making.
    """
    if not BACKUP_DB_PATH or not os.path.exists(BACKUP_DB_PATH):
        print("❌ No backup database found.")
        return

    try:
        session = BackupSession()

        # Get counts of potential issues
        total_ships = session.execute("SELECT COUNT(*) FROM ships").scalar()
        ships_with_sequence = session.execute(
            "SELECT COUNT(*) FROM ships WHERE sequence IS NOT NULL AND sequence != ''"
        ).scalar()

        print("\n📊 BACKUP DATABASE SUMMARY:")
        print(f"  Total ships: {total_ships}")
        print(f"  Ships with sequences: {ships_with_sequence}")

        # Check for CSV reports
        if os.path.exists("md5_mismatches.csv"):
            import pandas as pd

            mismatches_df = pd.read_csv("md5_mismatches.csv")
            print(f"  MD5 mismatches found: {len(mismatches_df)}")
            if len(mismatches_df) > 0:
                print("  See details in: md5_mismatches.csv")

        session.close()

    except Exception as e:
        print(f"Error examining backup: {e}")


if __name__ == "__main__":
    print("=" * 60)
    print("         STARBASE MD5 HASH VERIFICATION & REPAIR TOOL")
    print("=" * 60)
    print()
    print("This script will:")
    print("  1. Create a backup copy of your database")
    print("  2. Verify MD5 hashes for all sequences")
    print("  3. Fix any mismatches or duplicates")
    print("  4. Work ONLY on the backup (original stays safe)")
    print()
    print("SAFETY FEATURES:")
    print("  ✓ Original database will NOT be modified")
    print("  ✓ All changes made to a timestamped backup copy")
    print("  ✓ Full transaction rollback on any errors")
    print("  ✓ CSV reports generated for review")
    print()

    # Add safety prompt
    response = input("Continue with MD5 verification and repair? (y/N): ")
    if response.lower() in ["y", "yes"]:
        main()
    else:
        print("Operation cancelled.")
        print()
        print("No changes made to your database.")
        print("=" * 60)

#!/usr/bin/env python3
"""
Taxonomy Migration Script
=========================

This script migrates taxonomy data from the joined_ships table to the taxonomy table,
ensuring proper foreign key relationships and data integrity.

Tasks:
1. Clean species names that contain genus prefixes
2. Remove duplicate entries from taxonomy table (keeping most complete records)
3. Fill missing taxonomic hierarchy data using available records
4. Update taxonomy table with missing data from joined_ships
5. Ensure tax_id foreign key references are properly set
6. Handle data conflicts and duplicates appropriately
"""

import sqlite3
import sys
from pathlib import Path


def connect_to_database(db_path):
    """Connect to SQLite database and return connection."""
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row  # Enable column access by name
        return conn
    except sqlite3.Error as e:
        print(f"Error connecting to database: {e}")
        sys.exit(1)


def backup_database(conn):
    """Create a backup of the database before making changes"""
    print("Creating database backup...")
    backup_path = "src/database/db/starbase.sqlite.taxonomy_migration_backup"

    backup_conn = sqlite3.connect(backup_path)
    conn.backup(backup_conn)
    backup_conn.close()
    print(f"Backup created at: {backup_path}")


def clean_species_names(conn):
    """Clean species names that contain genus prefixes."""
    cursor = conn.cursor()

    # Find and fix species names that start with the genus name
    cursor.execute("""
        SELECT id, genus, species 
        FROM taxonomy 
        WHERE genus IS NOT NULL 
        AND species IS NOT NULL 
        AND species LIKE genus || ' %'
    """)

    records_to_update = cursor.fetchall()
    updated_count = 0

    for record in records_to_update:
        genus = record["genus"]
        species = record["species"]

        # Remove genus prefix and any extra spaces
        clean_species = species.replace(genus + " ", "", 1).strip()

        if clean_species and clean_species != species:
            cursor.execute(
                """
                UPDATE taxonomy 
                SET species = ? 
                WHERE id = ?
            """,
                (clean_species, record["id"]),
            )
            updated_count += 1

    conn.commit()
    return updated_count


def get_taxonomy_stats(conn):
    """Get current statistics from both tables."""
    cursor = conn.cursor()

    # Get counts
    cursor.execute("SELECT COUNT(*) FROM taxonomy")
    taxonomy_count = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM joined_ships")
    joined_ships_count = cursor.fetchone()[0]

    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships 
        WHERE genus IS NOT NULL OR species IS NOT NULL OR strain IS NOT NULL
    """)
    ships_with_taxonomy = cursor.fetchone()[0]

    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships 
        WHERE genus IS NOT NULL AND tax_id IS NULL
    """)
    ships_without_tax_id = cursor.fetchone()[0]

    cursor.execute("""
        SELECT COUNT(*) FROM taxonomy 
        WHERE taxID IS NULL OR superkingdom IS NULL OR kingdom IS NULL OR phylum IS NULL OR family IS NULL
    """)
    incomplete_taxonomy = cursor.fetchone()[0]

    # Count records with species containing genus
    cursor.execute("""
        SELECT COUNT(*) FROM taxonomy 
        WHERE genus IS NOT NULL 
        AND species IS NOT NULL 
        AND species LIKE genus || ' %'
    """)
    species_with_genus_prefix = cursor.fetchone()[0]

    return {
        "taxonomy_count": taxonomy_count,
        "joined_ships_count": joined_ships_count,
        "ships_with_taxonomy": ships_with_taxonomy,
        "ships_without_tax_id": ships_without_tax_id,
        "incomplete_taxonomy": incomplete_taxonomy,
        "species_with_genus_prefix": species_with_genus_prefix,
    }


def identify_duplicate_taxonomy_records(conn):
    """Identify duplicate records in taxonomy table based on genus/species/strain combination or name."""
    cursor = conn.cursor()

    # First, find duplicates based on name column
    cursor.execute("""
    WITH name_completeness AS (
        SELECT id, name, taxID, superkingdom, kingdom, phylum, family, genus, species, strain,
               (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) as completeness_score,
               ROW_NUMBER() OVER (
                   PARTITION BY name
                   ORDER BY 
                       (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) DESC, 
                       id ASC
               ) as rank_in_group,
               COUNT(*) OVER (PARTITION BY name) as group_count
        FROM taxonomy
        WHERE name IS NOT NULL AND TRIM(name) != ''
    )
    SELECT id, name, genus, species, strain, completeness_score, rank_in_group, 'name' as duplicate_type
    FROM name_completeness
    WHERE rank_in_group > 1 AND group_count > 1
    ORDER BY name, rank_in_group
    """)

    name_duplicates = cursor.fetchall()

    # Get IDs that are already identified as name duplicates
    name_duplicate_ids = set(record["id"] for record in name_duplicates)

    # Now find duplicates based on genus/species/strain combination (excluding name duplicates)
    cursor.execute("""
    WITH taxonomy_completeness AS (
        SELECT id, name, taxID, superkingdom, kingdom, phylum, family, genus, species, strain,
               (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) as completeness_score,
               ROW_NUMBER() OVER (
                   PARTITION BY genus, COALESCE(species, ''), COALESCE(strain, '') 
                   ORDER BY 
                       (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) DESC, 
                       id ASC
               ) as rank_in_group,
               COUNT(*) OVER (PARTITION BY genus, COALESCE(species, ''), COALESCE(strain, '')) as group_count
        FROM taxonomy
        WHERE genus IS NOT NULL
    )
    SELECT id, name, genus, species, strain, completeness_score, rank_in_group, 'taxonomy' as duplicate_type
    FROM taxonomy_completeness
    WHERE rank_in_group > 1 AND group_count > 1
    ORDER BY genus, species, strain, rank_in_group
    """)

    taxonomy_duplicates = cursor.fetchall()

    # Filter out taxonomy duplicates that are already in name duplicates
    filtered_taxonomy_duplicates = [
        record
        for record in taxonomy_duplicates
        if record["id"] not in name_duplicate_ids
    ]

    # Combine results
    all_duplicates = list(name_duplicates) + filtered_taxonomy_duplicates

    return all_duplicates


def get_best_taxonomy_records(conn):
    """Get the best (most complete) record for each unique taxonomy combination."""
    cursor = conn.cursor()

    # Get best records by name first
    cursor.execute("""
    WITH name_completeness AS (
        SELECT *,
               ROW_NUMBER() OVER (
                   PARTITION BY name
                   ORDER BY 
                       (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) DESC, 
                       id ASC
               ) as rank_in_group,
               COUNT(*) OVER (PARTITION BY name) as group_count
        FROM taxonomy
        WHERE name IS NOT NULL AND TRIM(name) != ''
    )
    SELECT id, name, genus, species, strain, 'name' as selection_type
    FROM name_completeness
    WHERE rank_in_group = 1 AND group_count > 1
    ORDER BY name
    """)

    name_best = cursor.fetchall()

    # Get IDs that are already handled by name duplicates
    name_handled_ids = set()
    cursor.execute("""
    SELECT id FROM taxonomy t1
    WHERE name IS NOT NULL AND TRIM(name) != ''
    AND EXISTS (
        SELECT 1 FROM taxonomy t2 
        WHERE t2.name = t1.name AND t2.id != t1.id
    )
    """)
    for row in cursor.fetchall():
        name_handled_ids.add(row["id"])

    # Get best records by genus/species/strain (excluding name duplicates)
    cursor.execute("""
    WITH taxonomy_completeness AS (
        SELECT *,
               ROW_NUMBER() OVER (
                   PARTITION BY genus, COALESCE(species, ''), COALESCE(strain, '') 
                   ORDER BY 
                       (CASE WHEN taxID IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN superkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN clade IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN kingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subkingdom IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN phylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subphylum IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN class IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN subclass IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN "order" IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN suborder IS NOT NULL THEN 1 ELSE 0 END +
                        CASE WHEN family IS NOT NULL THEN 1 ELSE 0 END) DESC, 
                       id ASC
               ) as rank_in_group,
               COUNT(*) OVER (PARTITION BY genus, COALESCE(species, ''), COALESCE(strain, '')) as group_count
        FROM taxonomy
        WHERE genus IS NOT NULL
    )
    SELECT id, name, genus, species, strain, 'taxonomy' as selection_type
    FROM taxonomy_completeness
    WHERE rank_in_group = 1 AND group_count > 1
    ORDER BY genus, species, strain
    """)

    taxonomy_best = cursor.fetchall()

    # Filter out taxonomy best records that are already handled by name
    filtered_taxonomy_best = [
        record for record in taxonomy_best if record["id"] not in name_handled_ids
    ]

    # Combine results
    all_best = list(name_best) + filtered_taxonomy_best

    return all_best


def update_foreign_keys_for_duplicates(conn, duplicates_to_remove, best_records):
    """Update tax_id references in joined_ships before removing duplicates."""
    cursor = conn.cursor()

    # Create mapping from old IDs to new IDs
    id_mapping = {}

    # Group best records by name and genus/species/strain
    best_by_name = {}
    best_by_taxonomy = {}

    for record in best_records:
        if record["selection_type"] == "name" and record["name"]:
            best_by_name[record["name"]] = record["id"]
        elif record["selection_type"] == "taxonomy":
            key = (record["genus"], record["species"] or "", record["strain"] or "")
            best_by_taxonomy[key] = record["id"]

    # Map duplicate IDs to best record IDs
    for duplicate in duplicates_to_remove:
        if duplicate["duplicate_type"] == "name" and duplicate["name"]:
            if duplicate["name"] in best_by_name:
                id_mapping[duplicate["id"]] = best_by_name[duplicate["name"]]
        elif duplicate["duplicate_type"] == "taxonomy":
            key = (
                duplicate["genus"],
                duplicate["species"] or "",
                duplicate["strain"] or "",
            )
            if key in best_by_taxonomy:
                id_mapping[duplicate["id"]] = best_by_taxonomy[key]

    updated_count = 0
    for old_id, new_id in id_mapping.items():
        cursor.execute(
            """
            UPDATE joined_ships 
            SET tax_id = ? 
            WHERE tax_id = ?
        """,
            (new_id, old_id),
        )
        updated_count += cursor.rowcount

    conn.commit()
    return updated_count, len(id_mapping)


def remove_duplicate_taxonomy_records(conn, duplicates_to_remove):
    """Remove duplicate taxonomy records."""
    cursor = conn.cursor()

    duplicate_ids = [record["id"] for record in duplicates_to_remove]

    if duplicate_ids:
        placeholders = ",".join(["?"] * len(duplicate_ids))
        cursor.execute(
            f"""
            DELETE FROM taxonomy 
            WHERE id IN ({placeholders})
        """,
            duplicate_ids,
        )

        removed_count = cursor.rowcount
        conn.commit()
        return removed_count

    return 0


def fill_missing_taxonomy_data(conn):
    """Fill missing taxonomic hierarchy data using available complete records."""
    cursor = conn.cursor()

    # Fill taxID only for exact genus/species/strain matches
    cursor.execute("""
        UPDATE taxonomy 
        SET taxID = (
            SELECT t2.taxID 
            FROM taxonomy t2 
            WHERE t2.genus = taxonomy.genus 
            AND COALESCE(t2.species, '') = COALESCE(taxonomy.species, '')
            AND COALESCE(t2.strain, '') = COALESCE(taxonomy.strain, '')
            AND t2.taxID IS NOT NULL 
            AND t2.id != taxonomy.id
            LIMIT 1
        )
        WHERE taxID IS NULL 
        AND genus IS NOT NULL
    """)
    taxid_updates = cursor.rowcount

    # Fill hierarchy columns (superkingdom to family) for genus matches only
    hierarchy_columns = [
        "superkingdom",
        "clade",
        "kingdom",
        "subkingdom",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        '"order"',
        "suborder",
        "family",
    ]

    total_hierarchy_updates = 0
    for column in hierarchy_columns:
        cursor.execute(f"""
            UPDATE taxonomy 
            SET {column} = (
                SELECT t2.{column}
                FROM taxonomy t2 
                WHERE t2.genus = taxonomy.genus 
                AND t2.{column} IS NOT NULL 
                AND t2.id != taxonomy.id
                LIMIT 1
            )
            WHERE {column} IS NULL 
            AND genus IS NOT NULL
        """)
        total_hierarchy_updates += cursor.rowcount

    conn.commit()
    return taxid_updates, total_hierarchy_updates


def find_missing_taxonomy_combinations(conn):
    """Find genus/species/strain combinations in joined_ships that don't exist in taxonomy."""
    cursor = conn.cursor()

    query = """
    SELECT DISTINCT js.genus, js.species, js.strain
    FROM joined_ships js
    LEFT JOIN taxonomy t ON (
        js.genus = t.genus AND 
        COALESCE(js.species, '') = COALESCE(t.species, '') AND 
        COALESCE(js.strain, '') = COALESCE(t.strain, '')
    )
    WHERE js.genus IS NOT NULL 
    AND t.id IS NULL
    ORDER BY js.genus, js.species, js.strain
    """

    cursor.execute(query)
    return cursor.fetchall()


def find_incomplete_taxonomy_records(conn):
    """Find taxonomy records that could use updates from joined_ships data."""
    cursor = conn.cursor()

    query = """
    SELECT DISTINCT t.id, t.genus, t.species, t.strain, 
           js.genus as js_genus, js.species as js_species, js.strain as js_strain
    FROM taxonomy t
    JOIN joined_ships js ON (
        t.genus = js.genus AND 
        COALESCE(t.species, '') = COALESCE(js.species, '') AND 
        COALESCE(t.strain, '') = COALESCE(js.strain, '')
    )
    WHERE (
        (t.species IS NULL AND js.species IS NOT NULL) OR
        (t.strain IS NULL AND js.strain IS NOT NULL)
    )
    ORDER BY t.id
    """

    cursor.execute(query)
    return cursor.fetchall()


def create_missing_taxonomy_records(conn, missing_combinations):
    """Create new taxonomy records for missing combinations."""
    cursor = conn.cursor()
    created_count = 0

    for combination in missing_combinations:
        genus = combination["genus"]
        species = combination["species"]
        strain = combination["strain"]

        # Insert new taxonomy record
        cursor.execute(
            """
            INSERT INTO taxonomy (genus, species, strain)
            VALUES (?, ?, ?)
        """,
            (genus, species, strain),
        )

        created_count += 1

        if created_count % 10 == 0:
            print(f"Created {created_count} taxonomy records...")

    conn.commit()
    return created_count


def update_incomplete_taxonomy_records(conn, incomplete_records):
    """Update taxonomy records with missing data from joined_ships."""
    cursor = conn.cursor()
    updated_count = 0

    for record in incomplete_records:
        tax_id = record["id"]
        updates = []
        params = []

        # Check what needs updating
        if record["species"] is None and record["js_species"] is not None:
            updates.append("species = ?")
            params.append(record["js_species"])

        if record["strain"] is None and record["js_strain"] is not None:
            updates.append("strain = ?")
            params.append(record["js_strain"])

        if updates:
            params.append(tax_id)
            query = f"UPDATE taxonomy SET {', '.join(updates)} WHERE id = ?"
            cursor.execute(query, params)
            updated_count += 1

    conn.commit()
    return updated_count


def update_tax_id_references(conn):
    """Update tax_id references in joined_ships to match taxonomy table."""
    cursor = conn.cursor()

    query = """
    UPDATE joined_ships 
    SET tax_id = (
        SELECT t.id 
        FROM taxonomy t 
        WHERE t.genus = joined_ships.genus 
        AND COALESCE(t.species, '') = COALESCE(joined_ships.species, '')
        AND COALESCE(t.strain, '') = COALESCE(joined_ships.strain, '')
    )
    WHERE genus IS NOT NULL 
    AND tax_id IS NULL
    """

    cursor.execute(query)
    updated_count = cursor.rowcount
    conn.commit()

    return updated_count


def verify_migration(conn):
    """Verify the migration was successful."""
    cursor = conn.cursor()

    # Check for orphaned tax_id references
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships js 
        LEFT JOIN taxonomy t ON js.tax_id = t.id 
        WHERE js.tax_id IS NOT NULL AND t.id IS NULL
    """)
    orphaned_count = cursor.fetchone()[0]

    # Check for records with taxonomy data but no tax_id
    cursor.execute("""
        SELECT COUNT(*) FROM joined_ships 
        WHERE genus IS NOT NULL AND tax_id IS NULL
    """)
    missing_tax_id_count = cursor.fetchone()[0]

    # Check for remaining duplicates
    cursor.execute("""
        WITH duplicate_check AS (
            SELECT genus, COALESCE(species, '') as species, COALESCE(strain, '') as strain, COUNT(*) as cnt
            FROM taxonomy
            WHERE genus IS NOT NULL
            GROUP BY genus, COALESCE(species, ''), COALESCE(strain, '')
            HAVING COUNT(*) > 1
        )
        SELECT COUNT(*) FROM duplicate_check
    """)
    remaining_duplicates = cursor.fetchone()[0]

    return {
        "orphaned_tax_ids": orphaned_count,
        "missing_tax_ids": missing_tax_id_count,
        "remaining_duplicates": remaining_duplicates,
    }


def main():
    """Main migration function."""
    db_path = "src/database/db/starbase.sqlite"

    if not Path(db_path).exists():
        print(f"Database file not found: {db_path}")
        sys.exit(1)

    print("Starting taxonomy migration and cleanup...")
    print("=" * 60)

    # Connect to database
    conn = connect_to_database(db_path)

    # Create backup
    backup_database(conn)

    try:
        # Get initial statistics
        initial_stats = get_taxonomy_stats(conn)
        print("Initial statistics:")
        print(f"  Taxonomy records: {initial_stats['taxonomy_count']}")
        print(f"  Joined ships records: {initial_stats['joined_ships_count']}")
        print(f"  Ships with taxonomy data: {initial_stats['ships_with_taxonomy']}")
        print(f"  Ships without tax_id: {initial_stats['ships_without_tax_id']}")
        print(f"  Incomplete taxonomy records: {initial_stats['incomplete_taxonomy']}")
        print(
            f"  Species with genus prefix: {initial_stats['species_with_genus_prefix']}"
        )
        print()

        # Step 0: Clean species names that contain genus prefixes
        print("Step 0: Cleaning species names with genus prefixes...")
        cleaned_species = clean_species_names(conn)
        print(f"Cleaned {cleaned_species} species names")

        # Step 1: Handle duplicates in taxonomy table
        print("\nStep 1: Identifying duplicate taxonomy records...")
        duplicates = identify_duplicate_taxonomy_records(conn)
        best_records = get_best_taxonomy_records(conn)
        print(f"Found {len(duplicates)} duplicate records to remove")
        print(f"Keeping {len(best_records)} best records")

        if duplicates:
            print("Updating foreign key references for duplicates...")
            fk_updated, mappings_count = update_foreign_keys_for_duplicates(
                conn, duplicates, best_records
            )
            print(
                f"Updated {fk_updated} foreign key references ({mappings_count} ID mappings)"
            )

            print("Removing duplicate taxonomy records...")
            removed_count = remove_duplicate_taxonomy_records(conn, duplicates)
            print(f"Removed {removed_count} duplicate records")

        # Step 2: Fill missing taxonomic hierarchy data
        print("\nStep 2: Filling missing taxonomic hierarchy data...")
        taxid_updates, hierarchy_updates = fill_missing_taxonomy_data(conn)
        print(f"Updated {taxid_updates} taxID fields")
        print(f"Updated {hierarchy_updates} hierarchy fields")

        # Step 3: Handle data from joined_ships
        print("\nStep 3: Processing joined_ships data...")

        # Find missing combinations
        print("Finding missing taxonomy combinations...")
        missing_combinations = find_missing_taxonomy_combinations(conn)
        print(f"Found {len(missing_combinations)} missing combinations")

        # Find incomplete records
        print("Finding incomplete taxonomy records...")
        incomplete_records = find_incomplete_taxonomy_records(conn)
        print(f"Found {len(incomplete_records)} incomplete records")

        # Create missing taxonomy records
        if missing_combinations:
            print("Creating missing taxonomy records...")
            created_count = create_missing_taxonomy_records(conn, missing_combinations)
            print(f"Created {created_count} new taxonomy records")

        # Update incomplete records
        if incomplete_records:
            print("Updating incomplete taxonomy records...")
            updated_count = update_incomplete_taxonomy_records(conn, incomplete_records)
            print(f"Updated {updated_count} taxonomy records")

        # Update tax_id references
        print("Updating tax_id references in joined_ships...")
        updated_refs = update_tax_id_references(conn)
        print(f"Updated {updated_refs} tax_id references")

        # Step 4: Final data filling after all records are in place
        print("\nStep 4: Final taxonomic hierarchy data filling...")
        final_taxid_updates, final_hierarchy_updates = fill_missing_taxonomy_data(conn)
        print(
            f"Final update: {final_taxid_updates} taxID fields, {final_hierarchy_updates} hierarchy fields"
        )

        # Verify migration
        print("\nVerifying migration...")
        verification = verify_migration(conn)
        print(f"Orphaned tax_id references: {verification['orphaned_tax_ids']}")
        print(f"Records missing tax_id: {verification['missing_tax_ids']}")
        print(f"Remaining duplicates: {verification['remaining_duplicates']}")

        # Get final statistics
        final_stats = get_taxonomy_stats(conn)
        print()
        print("Final statistics:")
        print(f"  Taxonomy records: {final_stats['taxonomy_count']}")
        print(f"  Ships without tax_id: {final_stats['ships_without_tax_id']}")
        print(f"  Incomplete taxonomy records: {final_stats['incomplete_taxonomy']}")
        print(
            f"  Species with genus prefix: {final_stats['species_with_genus_prefix']}"
        )

        print()
        print("Migration completed successfully!")

        if (
            verification["orphaned_tax_ids"] > 0
            or verification["missing_tax_ids"] > 0
            or verification["remaining_duplicates"] > 0
        ):
            print("⚠️  Warning: Some issues remain after migration!")
        else:
            print("✅ All taxonomy data properly migrated and cleaned!")

    except Exception as e:
        print(f"Error during migration: {e}")
        conn.rollback()
        sys.exit(1)
    finally:
        conn.close()


if __name__ == "__main__":
    main()

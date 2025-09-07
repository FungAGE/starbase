#!/usr/bin/env python3
"""
Taxonomy SQLite Updater

This utility updates an SQLite taxonomy table with taxonomic information
from a TSV file containing JSON taxonomy data.
"""

import sqlite3
import csv
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


class TaxonomyUpdater:
    def __init__(self, db_path: str):
        """
        Initialize the TaxonomyUpdater with database connection.
        
        Args:
            db_path: Path to the SQLite database file
        """
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Enable column access by name
        
        # Taxonomy hierarchy order (from highest to lowest)
        self.hierarchy_order = [
            'superkingdom', 'clade', 'kingdom', 'subkingdom', 'phylum', 
            'subphylum', 'class', 'subclass', '`order`', 'suborder', 
            'family', 'genus', 'species', 'section', 'species_group', 
            'subgenus', 'strain'
        ]
        
        # Mapping from JSON keys to database columns
        self.json_to_db_mapping = {
            'kingdom': 'kingdom',
            'phylum': 'phylum',
            'subphylum': 'subphylum',
            'class': 'class',
            '`order`': '`order`',
            'family': 'family',
            'subfamily': 'family',  # Note: subfamily maps to family for simplicity
            'genus': 'genus',
            'species': 'species',
            'strain': 'strain'
        }

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.conn.close()

    def get_existing_taxonomies(self) -> Set[Tuple]:
        """
        Get all existing taxonomy entries as tuples for duplicate detection.
        
        Returns:
            Set of tuples representing existing taxonomy combinations
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT kingdom, phylum, subphylum, class, `order`, family, 
                   genus, species, strain
            FROM taxonomy
        """)
        
        existing = set()
        for row in cursor.fetchall():
            # Convert None values to empty strings for consistent comparison
            taxonomy_tuple = tuple(val if val is not None else '' for val in row)
            existing.add(taxonomy_tuple)
        
        return existing

    def parse_taxonomy_json(self, json_str: str) -> Dict[str, str]:
        """
        Parse taxonomy JSON string and map to database columns.
        
        Args:
            json_str: JSON string containing taxonomy information
            
        Returns:
            Dictionary mapping database columns to values
        """
        try:
            taxonomy_data = json.loads(json_str)
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON: {e}")
            return {}
        
        # Map JSON fields to database columns
        mapped_taxonomy = {}
        for json_key, db_column in self.json_to_db_mapping.items():
            value = taxonomy_data.get(json_key, '')
            # Clean up the value
            if value and value.strip():
                mapped_taxonomy[db_column] = value.strip()
            else:
                mapped_taxonomy[db_column] = ''
        
        return mapped_taxonomy

    def extract_species_name(self, species_full: str) -> str:
        """
        Extract species epithet from full species name.
        
        Args:
            species_full: Full species name (e.g., "Alternaria alternata")
            
        Returns:
            Species epithet (e.g., "alternata")
        """
        if not species_full:
            return ''
        
        parts = species_full.strip().split()
        if len(parts) >= 2:
            return parts[1]  # Second part is the species epithet
        return species_full

    def create_taxonomy_entry(self, taxonomy: Dict[str, str]) -> Dict[str, str]:
        """
        Create a complete taxonomy entry with proper formatting.
        
        Args:
            taxonomy: Dictionary with taxonomy information
            
        Returns:
            Complete taxonomy entry ready for database insertion
        """
        # Extract species epithet from full species name
        if 'species' in taxonomy and taxonomy['species']:
            species_full = taxonomy['species']
            taxonomy['species'] = self.extract_species_name(species_full)
        
        # Create full taxonomy entry with all required fields
        full_entry = {
            'name': '',  # Will be constructed from genus + species + strain
            'taxID': None,  # Not provided in source data
            'superkingdom': 'Eukaryota',  # Assuming all fungi are eukaryotic
            'clade': 'Opisthokonta',  # Assuming all fungi belong to this clade
            'kingdom': taxonomy.get('kingdom', ''),
            'subkingdom': 'Dikarya',  # Common for most fungi, can be adjusted
            'phylum': taxonomy.get('phylum', ''),
            'subphylum': taxonomy.get('subphylum', ''),
            'class': taxonomy.get('class', ''),
            'subclass': '',  # Not provided in source
            '`order`': taxonomy.get('`order`', ''),
            'suborder': '',  # Not provided in source
            'family': taxonomy.get('family', ''),
            'genus': taxonomy.get('genus', ''),
            'species': taxonomy.get('species', ''),
            'section': '',  # Not provided in source
            'species_group': '',  # Not provided in source
            'subgenus': '',  # Not provided in source
            'strain': taxonomy.get('strain', '')
        }
        
        # Construct the name field
        name_parts = []
        if full_entry['genus']:
            name_parts.append(full_entry['genus'])
        if full_entry['species']:
            name_parts.append(full_entry['species'])
        if full_entry['strain']:
            name_parts.append(full_entry['strain'])
        
        full_entry['name'] = ' '.join(name_parts)
        
        return full_entry

    def is_duplicate(self, taxonomy: Dict[str, str], existing: Set[Tuple]) -> bool:
        """
        Check if taxonomy entry already exists in the database.
        
        Args:
            taxonomy: Taxonomy entry to check
            existing: Set of existing taxonomy tuples
            
        Returns:
            True if duplicate exists, False otherwise
        """
        # Create tuple in same order as existing taxonomies query
        taxonomy_tuple = (
            taxonomy.get('kingdom', '') or '',
            taxonomy.get('phylum', '') or '',
            taxonomy.get('subphylum', '') or '',
            taxonomy.get('class', '') or '',
            taxonomy.get('`order`', '') or '',
            taxonomy.get('family', '') or '',
            taxonomy.get('genus', '') or '',
            taxonomy.get('species', '') or '',
            taxonomy.get('strain', '') or ''
        )
        
        return taxonomy_tuple in existing

    def insert_taxonomy(self, taxonomy: Dict[str, str]) -> int:
        """
        Insert a new taxonomy entry into the database.
        
        Args:
            taxonomy: Complete taxonomy entry
            
        Returns:
            ID of the inserted record
        """
        cursor = self.conn.cursor()
        
        # Prepare the insert statement
        columns = list(taxonomy.keys())
        placeholders = ', '.join(['?' for _ in columns])
        column_names = ', '.join(columns)
        
        query = f"""
            INSERT INTO taxonomy ({column_names})
            VALUES ({placeholders})
        """
        
        values = [taxonomy[col] if taxonomy[col] else None for col in columns]
        
        cursor.execute(query, values)
        self.conn.commit()
        
        return cursor.lastrowid

    def process_tsv_file(self, tsv_path: str, dry_run: bool = False) -> Tuple[int, int]:
        """
        Process the TSV file and update the taxonomy table.
        
        Args:
            tsv_path: Path to the TSV file
            dry_run: If True, don't actually insert records
            
        Returns:
            Tuple of (processed_count, inserted_count)
        """
        if not Path(tsv_path).exists():
            raise FileNotFoundError(f"TSV file not found: {tsv_path}")
        
        existing_taxonomies = self.get_existing_taxonomies()
        processed_count = 0
        inserted_count = 0
        
        with open(tsv_path, 'r', encoding='utf-8') as file:
            # Use csv.reader with tab delimiter
            reader = csv.reader(file, delimiter='\t')
            
            for row_num, row in enumerate(reader, 1):
                try:
                    if len(row) < 5:  # Ensure we have enough columns
                        print(f"Row {row_num}: Insufficient columns, skipping")
                        continue
                    
                    # Extract taxonomy JSON (assumed to be in column 4, 0-indexed)
                    taxonomy_json = row[4]
                    
                    # Parse the taxonomy data
                    taxonomy_data = self.parse_taxonomy_json(taxonomy_json)
                    if not taxonomy_data:
                        print(f"Row {row_num}: Failed to parse taxonomy JSON, skipping")
                        continue
                    
                    # Create complete taxonomy entry
                    full_taxonomy = self.create_taxonomy_entry(taxonomy_data)
                    
                    processed_count += 1
                    
                    # Check for duplicates
                    if self.is_duplicate(full_taxonomy, existing_taxonomies):
                        print(f"Row {row_num}: Duplicate entry found, skipping: {full_taxonomy['name']}")
                        continue
                    
                    # Create tuple for tracking (same format as is_duplicate uses)
                    taxonomy_tuple = (
                        full_taxonomy.get('kingdom', '') or '',
                        full_taxonomy.get('phylum', '') or '',
                        full_taxonomy.get('subphylum', '') or '',
                        full_taxonomy.get('class', '') or '',
                        full_taxonomy.get('`order`', '') or '',
                        full_taxonomy.get('family', '') or '',
                        full_taxonomy.get('genus', '') or '',
                        full_taxonomy.get('species', '') or '',
                        full_taxonomy.get('strain', '') or ''
                    )
                    
                    if not dry_run:
                        # Insert the new taxonomy entry
                        new_id = self.insert_taxonomy(full_taxonomy)
                        print(f"Row {row_num}: Inserted new taxonomy entry (ID: {new_id}): {full_taxonomy['name']}")
                    else:
                        print(f"Row {row_num}: Would insert: {full_taxonomy['name']}")
                    
                    # Add to existing set to prevent duplicates within the same run
                    existing_taxonomies.add(taxonomy_tuple)
                    inserted_count += 1
                    
                except Exception as e:
                    print(f"Row {row_num}: Error processing row: {e}")
                    continue
        
        return processed_count, inserted_count


def main():
    parser = argparse.ArgumentParser(description='Update SQLite taxonomy table from TSV file')
    parser.add_argument('database', help='Path to SQLite database file')
    parser.add_argument('tsv_file', help='Path to TSV file with taxonomy data')
    parser.add_argument('--dry-run', action='store_true', 
                        help='Show what would be inserted without actually doing it')
    
    args = parser.parse_args()
    
    try:
        with TaxonomyUpdater(args.database) as updater:
            processed, inserted = updater.process_tsv_file(args.tsv_file, args.dry_run)
            
            if args.dry_run:
                print(f"\nDry run completed:")
                print(f"Processed: {processed} entries")
                print(f"Would insert: {inserted} new entries")
            else:
                print(f"\nUpdate completed:")
                print(f"Processed: {processed} entries")
                print(f"Inserted: {inserted} new entries")
                
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
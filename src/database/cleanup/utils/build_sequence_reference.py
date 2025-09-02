#!/usr/bin/env python3
"""
Standalone tool to build a sequence reference table from FASTA files.

This tool will:
1. Create a new table in the SQLite database to store sequence reference data
2. Parse FASTA files to populate this table with:
   - starshipID (from FASTA header)
   - sequence
   - MD5 hash of original sequence
   - MD5 hash of reverse complement sequence

This reference table can then be used to validate ship_id relationships in existing tables.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from difflib import SequenceMatcher

# Add the project root to the path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))

from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.utils.seq_utils import clean_sequence, revcomp
from src.utils.classification_utils import generate_md5_hash

logger = get_logger(__name__)


def calculate_similarity(str1: str, str2: str) -> float:
    """
    Calculate similarity between two strings using SequenceMatcher.
    
    Args:
        str1: First string
        str2: Second string
        
    Returns:
        Similarity ratio (0.0 to 1.0, where 1.0 is identical)
    """
    return SequenceMatcher(None, str1.lower(), str2.lower()).ratio()


def find_fuzzy_match(header: str, existing_starshipids: List[str], threshold: float = 0.85) -> Optional[str]:
    """
    Find a fuzzy match for the header among existing starshipIDs.
    
    Args:
        header: The FASTA header to match
        existing_starshipids: List of existing starshipIDs to search
        threshold: Minimum similarity ratio (default: 0.85)
        
    Returns:
        The best matching starshipID or None if no good match found
    """
    if not existing_starshipids:
        return None
    
    best_match = None
    best_similarity = 0.0
    
    for starshipid in existing_starshipids:
        similarity = calculate_similarity(header, starshipid)
        if similarity > best_similarity and similarity >= threshold:
            best_similarity = similarity
            best_match = starshipid
    
    return best_match


def get_existing_starshipids() -> List[str]:
    """
    Get all existing starshipIDs from the database for fuzzy matching.
    
    Returns:
        List of existing starshipIDs
    """
    session = StarbaseSession()
    
    try:
        # Get starshipIDs from joined_ships table
        try:
            joined_starshipids = session.execute(
                "SELECT DISTINCT starshipID FROM joined_ships WHERE starshipID IS NOT NULL"
            ).fetchall()
        except Exception as e:
            logger.debug(f"Could not query joined_ships table: {str(e)}")
            joined_starshipids = []
        
        # Get starshipIDs from sequence_reference table (if it exists)
        try:
            ref_starshipids = session.execute(
                "SELECT DISTINCT starshipID FROM sequence_reference WHERE starshipID IS NOT NULL"
            ).fetchall()
        except Exception as e:
            logger.debug(f"Could not query sequence_reference table (may not exist yet): {str(e)}")
            ref_starshipids = []
        
        # Combine and deduplicate
        all_starshipids = set()
        for row in joined_starshipids + ref_starshipids:
            if row[0]:
                all_starshipids.add(row[0])
        
        logger.debug(f"Found {len(all_starshipids)} existing starshipIDs for fuzzy matching")
        return list(all_starshipids)
        
    except Exception as e:
        logger.warning(f"Could not fetch existing starshipIDs: {str(e)}")
        return []
    finally:
        session.close()


def transform_fasta_header_to_starshipid(header: str, use_fuzzy_match: bool = True, fuzzy_threshold: float = 0.85) -> str:
    """
    Transform FASTA header to match the expected starshipID format.
    
    This function applies known transformations and optionally uses fuzzy matching
    to convert FASTA headers to the correct starshipID format used in the database.
    
    Args:
        header: The raw FASTA header
        use_fuzzy_match: Whether to use fuzzy matching for unknown headers
        fuzzy_threshold: Minimum similarity ratio for fuzzy matching
        
    Returns:
        The transformed starshipID
    """
    # Start with the header
    starshipID = header.strip()
    
    # Debug: Log if header becomes empty after stripping
    if not starshipID:
        logger.warning(f"Header became empty after stripping: original='{header}' (length: {len(header)})")
        return ""
    
    # Known transformations
    transformations = [
        # Replace "fam_" with "_" (e.g., "Prometheusfam_Mrub" -> "Prometheus_Mrub")
        ("fam_", "_"),
        
        # Add more transformations here as needed
        # ("old_pattern", "new_pattern"),
    ]
    
    # Apply all transformations
    for old_pattern, new_pattern in transformations:
        starshipID = starshipID.replace(old_pattern, new_pattern)
    
    # If fuzzy matching is enabled, try to find a better match
    if use_fuzzy_match:
        existing_starshipids = get_existing_starshipids()
        if existing_starshipids:
            fuzzy_match = find_fuzzy_match(starshipID, existing_starshipids, fuzzy_threshold)
            if fuzzy_match:
                logger.debug(f"Fuzzy match found: '{starshipID}' -> '{fuzzy_match}'")
                return fuzzy_match
    
    # Debug: Log the final result
    logger.debug(f"Final starshipID: '{starshipID}' (type: {type(starshipID)})")
    
    # Ensure we never return None
    if starshipID is None:
        logger.error(f"starshipID became None! Original header: '{header}'")
        return ""
    
    return starshipID


def create_sequence_reference_table():
    """
    Create the sequence_reference table in the database.
    
    This table will store:
    - id: Primary key
    - starshipID: From FASTA header
    - sequence: Cleaned sequence
    - md5_original: MD5 hash of original sequence
    - md5_revcomp: MD5 hash of reverse complement sequence
    - canonical_md5: Lexicographically smaller of the two MD5s
    - created_at: Timestamp
    """
    logger.info("Creating sequence_reference table...")
    session = StarbaseSession()
    
    try:
        # Create the table using raw SQL
        create_table_sql = """
        CREATE TABLE IF NOT EXISTS sequence_reference (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            starshipID TEXT NOT NULL UNIQUE,
            sequence TEXT NOT NULL,
            md5_original TEXT NOT NULL,
            md5_revcomp TEXT NOT NULL,
            canonical_md5 TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        """
        
        session.execute(create_table_sql)
        session.commit()
        
        # Create indexes for efficient lookups
        indexes = [
            "CREATE INDEX IF NOT EXISTS idx_sequence_ref_starshipID ON sequence_reference(starshipID)",
            "CREATE INDEX IF NOT EXISTS idx_sequence_ref_canonical_md5 ON sequence_reference(canonical_md5)",
            "CREATE INDEX IF NOT EXISTS idx_sequence_ref_md5_original ON sequence_reference(md5_original)",
            "CREATE INDEX IF NOT EXISTS idx_sequence_ref_md5_revcomp ON sequence_reference(md5_revcomp)"
        ]
        
        for index_sql in indexes:
            session.execute(index_sql)
        
        session.commit()
        logger.info("Sequence reference table created successfully with indexes")
        
    except Exception as e:
        logger.error(f"Error creating sequence_reference table: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def parse_fasta_file(fasta_path: str, use_fuzzy_match: bool = True, fuzzy_threshold: float = 0.85, debug_headers: bool = False) -> List[Tuple[str, str]]:
    """
    Parse a FASTA file and extract (starshipID, sequence) pairs.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        List of (starshipID, sequence) tuples
    """
    logger.info(f"Parsing FASTA file: {fasta_path}")
    
    sequences = []
    current_header = None
    current_sequence = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_header and current_sequence:
                        sequence = ''.join(current_sequence)
                        # Transform the header to starshipID format
                        starshipID = transform_fasta_header_to_starshipid(current_header, use_fuzzy_match, fuzzy_threshold)
                        # Only add if starshipID is not empty
                        if starshipID and starshipID.strip():
                            sequences.append((starshipID, sequence))
                        else:
                            logger.warning(f"Skipping sequence with empty header: '{current_header}'")
                    
                    # Start new sequence
                    current_header = line[1:]  # Remove '>' character
                    current_sequence = []
                    
                    # Debug: Log problematic headers
                    if debug_headers:
                        if not current_header or current_header.strip() == '':
                            logger.warning(f"Found empty header at line {line_num}: '{line}' (length: {len(line)})")
                        elif len(current_header.strip()) < 3:  # Very short headers
                            logger.info(f"Found short header: '{current_header}' (length: {len(current_header)})")
                    
                elif line and current_header:
                    # Add sequence line
                    current_sequence.append(line)
                elif line and not current_header:
                    logger.warning(f"Line {line_num}: Found sequence data without header, skipping")
        
        # Don't forget the last sequence
        if current_header and current_sequence:
            sequence = ''.join(current_sequence)
            # Transform the header to starshipID format
            starshipID = transform_fasta_header_to_starshipid(current_header, use_fuzzy_match, fuzzy_threshold)
            # Only add if starshipID is not empty
            if starshipID and starshipID.strip():
                sequences.append((starshipID, sequence))
            else:
                logger.warning(f"Skipping final sequence with empty header: '{current_header}'")
        
        logger.info(f"Parsed {len(sequences)} sequences from {fasta_path}")
        return sequences
        
    except Exception as e:
        logger.error(f"Error parsing FASTA file {fasta_path}: {str(e)}")
        raise


def process_sequence(starshipID: str, sequence: str) -> Optional[Dict]:
    """
    Process a single sequence: clean it and generate MD5 hashes.
    
    Args:
        starshipID: The starshipID from FASTA header
        sequence: The raw sequence
        
    Returns:
        Dictionary with processed sequence data or None if processing failed
    """
    try:
        # Validate starshipID is not empty
        if not starshipID or starshipID.strip() == '':
            logger.warning(f"Skipping sequence with empty starshipID")
            return None
        
        # Clean the sequence
        clean_seq = clean_sequence(sequence)
        if clean_seq is None:
            logger.warning(f"Could not clean sequence for {starshipID}")
            return None
        
        # Generate MD5 hashes
        md5_original = generate_md5_hash(clean_seq)
        md5_revcomp = generate_md5_hash(revcomp(clean_seq))
        
        if md5_original is None or md5_revcomp is None:
            logger.warning(f"Could not generate MD5 for {starshipID}")
            return None
        
        # Determine canonical MD5 (lexicographically smaller)
        canonical_md5 = min(md5_original, md5_revcomp)
        
        result = {
            'starshipID': starshipID,
            'sequence': clean_seq,
            'md5_original': md5_original,
            'md5_revcomp': md5_revcomp,
            'canonical_md5': canonical_md5
        }
        
        # Debug: Log the processed data (truncate sequence)
        debug_result = result.copy()
        if len(debug_result['sequence']) > 50:
            debug_result['sequence'] = debug_result['sequence'][:50] + "..."
        logger.debug(f"Processed sequence for {starshipID}: {debug_result}")
        
        return result
        
    except Exception as e:
        logger.error(f"Error processing sequence for {starshipID}: {str(e)}")
        return None


def insert_sequence_data(sequence_data: List[Dict], batch_size: int = 1000):
    """
    Insert sequence data into the sequence_reference table.
    
    Args:
        sequence_data: List of dictionaries with sequence data
        batch_size: Number of records to insert in each batch
    """
    logger.info(f"Inserting {len(sequence_data)} sequences into sequence_reference table...")
    
    total_inserted = 0
    total_skipped = 0
    
    # Process in batches
    for i in range(0, len(sequence_data), batch_size):
        batch = sequence_data[i:i + batch_size]
        
        # Create a new session for each batch
        session = StarbaseSession()
        logger.debug(f"Created new session for batch {i//batch_size + 1}")
        
        try:
            for seq_data in batch:
                try:
                    # Debug: Log the actual data we're trying to insert (truncate sequence)
                    debug_seq_data = seq_data.copy()
                    if len(debug_seq_data['sequence']) > 50:
                        debug_seq_data['sequence'] = debug_seq_data['sequence'][:50] + "..."
                    logger.debug(f"Processing sequence data: {debug_seq_data}")
                    
                    # Validate starshipID is not None or empty
                    if seq_data['starshipID'] is None:
                        logger.warning(f"Skipping sequence with None starshipID: {seq_data}")
                        total_skipped += 1
                        continue
                    
                    if not seq_data['starshipID'] or seq_data['starshipID'].strip() == '':
                        logger.warning(f"Skipping sequence with empty starshipID: {seq_data}")
                        total_skipped += 1
                        continue
                    
                    # Check if starshipID already exists
                    starship_id_value = seq_data['starshipID']
                    logger.debug(f"About to execute SQL with starshipID: '{starship_id_value}' (type: {type(starship_id_value)})")
                    
                    # Create the parameter tuple
                    param_tuple = (starship_id_value,)
                    logger.debug(f"SQL parameter tuple: {param_tuple} (type: {type(param_tuple)}, length: {len(param_tuple)})")
                    
                    # Try using raw SQL instead of SQLAlchemy session
                    try:
                        existing = session.execute(
                            "SELECT id FROM sequence_reference WHERE starshipID = ?",
                            param_tuple
                        ).fetchone()
                    except Exception as e:
                        logger.error(f"SQLAlchemy session.execute failed: {str(e)}")
                        # Try alternative approach using raw sqlite3
                        logger.debug("Trying raw sqlite3 execution...")
                        try:
                            # Get the raw connection from SQLAlchemy session
                            raw_conn = session.connection().connection
                            cursor = raw_conn.cursor()
                            cursor.execute("SELECT id FROM sequence_reference WHERE starshipID = ?", param_tuple)
                            existing = cursor.fetchone()
                            cursor.close()
                        except Exception as raw_e:
                            logger.error(f"Raw sqlite3 also failed: {str(raw_e)}")
                            existing = None
                    
                    if existing:
                        logger.debug(f"StarshipID {seq_data['starshipID']} already exists, skipping")
                        total_skipped += 1
                        continue
                    
                    # Insert new record using raw sqlite3
                    try:
                        raw_conn = session.connection().connection
                        cursor = raw_conn.cursor()
                        
                        insert_sql = """
                        INSERT INTO sequence_reference 
                        (starshipID, sequence, md5_original, md5_revcomp, canonical_md5)
                        VALUES (?, ?, ?, ?, ?)
                        """
                        
                        insert_params = (
                            seq_data['starshipID'],
                            seq_data['sequence'],
                            seq_data['md5_original'],
                            seq_data['md5_revcomp'],
                            seq_data['canonical_md5']
                        )
                        
                        logger.debug(f"Inserting with raw sqlite3: {insert_params}")
                        cursor.execute(insert_sql, insert_params)
                        cursor.close()
                        
                    except Exception as insert_e:
                        logger.error(f"Raw sqlite3 insert failed: {str(insert_e)}")
                        total_skipped += 1
                        continue
                    
                    total_inserted += 1
                    
                except Exception as e:
                    logger.error(f"Error inserting {seq_data.get('starshipID', 'UNKNOWN')}: {str(e)}")
                    total_skipped += 1
            
            # Commit batch
            session.commit()
            logger.info(f"Processed batch {i//batch_size + 1}, inserted {total_inserted} so far")
            
        except Exception as e:
            logger.error(f"Error in batch {i//batch_size + 1}: {str(e)}")
            session.rollback()
        finally:
            session.close()
    
    logger.info(f"Insertion complete: {total_inserted} inserted, {total_skipped} skipped")


def get_table_stats():
    """
    Get statistics about the sequence_reference table.
    
    Returns:
        Dictionary with table statistics
    """
    session = StarbaseSession()
    
    try:
        # Get total count
        total_count = session.execute("SELECT COUNT(*) FROM sequence_reference").fetchone()[0]
        
        # Get unique canonical MD5 count
        unique_canonical = session.execute("SELECT COUNT(DISTINCT canonical_md5) FROM sequence_reference").fetchone()[0]
        
        # Get duplicate canonical MD5 groups
        duplicate_groups = session.execute("""
            SELECT canonical_md5, COUNT(*) as count 
            FROM sequence_reference 
            GROUP BY canonical_md5 
            HAVING COUNT(*) > 1
        """).fetchall()
        
        stats = {
            'total_sequences': total_count,
            'unique_canonical_md5': unique_canonical,
            'duplicate_groups': len(duplicate_groups),
            'duplicate_details': [
                {'canonical_md5': md5, 'count': count} 
                for md5, count in duplicate_groups
            ]
        }
        
        return stats
        
    except Exception as e:
        logger.error(f"Error getting table stats: {str(e)}")
        raise
    finally:
        session.close()


def test_single_sequence(fasta_files=None):
    """
    Test processing of a single sequence to debug the issue.
    """
    logger.info("Testing single sequence processing...")
    
    # Test with a simple sequence first
    test_header = "Test_Sequence_123"
    test_sequence = "ATCGATCGATCG"
    
    logger.info(f"Testing with header: '{test_header}'")
    logger.info(f"Testing with sequence: '{test_sequence}'")
    
    # Test transformation
    transformed = transform_fasta_header_to_starshipid(test_header, use_fuzzy_match=False)
    logger.info(f"Transformed header: '{transformed}'")
    
    # Test processing
    processed = process_sequence(transformed, test_sequence)
    logger.info(f"Processed result: {processed}")
    
    if processed:
        logger.info("Single sequence test PASSED")
    else:
        logger.error("Single sequence test FAILED")
    
    # Now test with a real sequence from the FASTA files if provided
    if fasta_files:
        logger.info("Testing with real FASTA file...")
        for fasta_file in fasta_files:
            if os.path.exists(fasta_file):
                logger.info(f"Testing with file: {fasta_file}")
                
                # Parse just the first sequence
                sequences = parse_fasta_file(fasta_file, use_fuzzy_match=False, fuzzy_threshold=0.85, debug_headers=True)
                
                if sequences:
                    first_starshipID, first_sequence = sequences[0]
                    logger.info(f"First sequence from file:")
                    logger.info(f"  starshipID: '{first_starshipID}' (type: {type(first_starshipID)})")
                    logger.info(f"  sequence length: {len(first_sequence)}")
                    
                    # Process it
                    processed = process_sequence(first_starshipID, first_sequence)
                    logger.info(f"Processed result: {processed}")
                    
                    if processed and processed.get('starshipID'):
                        logger.info("Real sequence test PASSED")
                    else:
                        logger.error("Real sequence test FAILED")
                        logger.error(f"Processed result: {processed}")
                else:
                    logger.warning(f"No sequences found in {fasta_file}")
                break  # Only test the first file


def main():
    """
    Main function to build the sequence reference table.
    """
    parser = argparse.ArgumentParser(description="Build sequence reference table from FASTA files")
    parser.add_argument("fasta_files", nargs="+", help="FASTA files to process")
    parser.add_argument("--create-table", action="store_true", 
                       help="Create the sequence_reference table (if not exists)")
    parser.add_argument("--batch-size", type=int, default=1000,
                       help="Batch size for database inserts (default: 1000)")
    parser.add_argument("--stats", action="store_true",
                       help="Show table statistics after processing")
    parser.add_argument("--dry-run", action="store_true",
                       help="Parse files but don't insert into database")
    parser.add_argument("--show-transformations", action="store_true",
                       help="Show header transformations that would be applied")
    parser.add_argument("--no-fuzzy-match", action="store_true",
                       help="Disable fuzzy matching for header transformations")
    parser.add_argument("--fuzzy-threshold", type=float, default=0.85,
                       help="Minimum similarity ratio for fuzzy matching (default: 0.85)")
    parser.add_argument("--debug-headers", action="store_true",
                       help="Show detailed debugging information about headers")
    parser.add_argument("--test-single", action="store_true",
                       help="Test processing of a single sequence for debugging")
    
    args = parser.parse_args()
    
    try:
        # Run single sequence test if requested
        if args.test_single:
            test_single_sequence(args.fasta_files)
            return
        
        # Create table if requested or if we're not in dry-run mode
        if args.create_table or (not args.dry_run):
            create_sequence_reference_table()
        
        # Process all FASTA files
        all_sequences = []
        header_transformations = {}  # Track transformations for display
        
        for fasta_file in args.fasta_files:
            if not os.path.exists(fasta_file):
                logger.error(f"FASTA file not found: {fasta_file}")
                continue
            
            # Parse FASTA file
            use_fuzzy = not args.no_fuzzy_match
            fasta_sequences = parse_fasta_file(fasta_file, use_fuzzy, args.fuzzy_threshold, args.debug_headers)
            
            # Process sequences
            for starshipID, sequence in fasta_sequences:
                # Debug: Log what we're processing
                logger.debug(f"Processing starshipID: '{starshipID}' (type: {type(starshipID)})")
                
                # Track transformations if requested
                if args.show_transformations:
                    # Re-parse to get original header for comparison
                    with open(fasta_file, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                original_header = line[1:].strip()
                                transformed = transform_fasta_header_to_starshipid(original_header, use_fuzzy, args.fuzzy_threshold)
                                if original_header != transformed:
                                    header_transformations[original_header] = transformed
                                break
                
                processed = process_sequence(starshipID, sequence)
                if processed:
                    # Debug: Log the processed data before adding to list (truncate sequence)
                    debug_processed = processed.copy()
                    if len(debug_processed['sequence']) > 50:
                        debug_processed['sequence'] = debug_processed['sequence'][:50] + "..."
                    logger.debug(f"Adding processed data: {debug_processed}")
                    all_sequences.append(processed)
                else:
                    logger.warning(f"Failed to process sequence for {starshipID}")
        
        logger.info(f"Successfully processed {len(all_sequences)} sequences from {len(args.fasta_files)} files")
        
        # Debug: Check the first few sequences (truncate sequence)
        if all_sequences:
            debug_first = all_sequences[0].copy()
            if len(debug_first['sequence']) > 50:
                debug_first['sequence'] = debug_first['sequence'][:50] + "..."
            logger.debug(f"First sequence in list: {debug_first}")
            logger.debug(f"First sequence starshipID: '{all_sequences[0].get('starshipID')}' (type: {type(all_sequences[0].get('starshipID'))})")
        
        # Show header transformations if requested
        if args.show_transformations and header_transformations:
            print("\n" + "=" * 60)
            print("FASTA HEADER TRANSFORMATIONS")
            print("=" * 60)
            for original, transformed in header_transformations.items():
                print(f"  '{original}' -> '{transformed}'")
            print(f"\nTotal transformations: {len(header_transformations)}")
        
        # Insert into database (unless dry run)
        if not args.dry_run and all_sequences:
            # Ensure table exists before inserting
            create_sequence_reference_table()
            insert_sequence_data(all_sequences, args.batch_size)
        
        # Show statistics
        if args.stats and not args.dry_run:
            stats = get_table_stats()
            print("\n" + "=" * 60)
            print("SEQUENCE REFERENCE TABLE STATISTICS")
            print("=" * 60)
            print(f"Total sequences: {stats['total_sequences']}")
            print(f"Unique canonical MD5: {stats['unique_canonical_md5']}")
            print(f"Duplicate groups: {stats['duplicate_groups']}")
            
            if stats['duplicate_details']:
                print("\nDuplicate groups (first 10):")
                for dup in stats['duplicate_details'][:10]:
                    print(f"  Canonical MD5: {dup['canonical_md5']} - {dup['count']} sequences")
                if len(stats['duplicate_details']) > 10:
                    print(f"  ... and {len(stats['duplicate_details']) - 10} more")
        
        logger.info("Sequence reference table build completed successfully")
        
    except Exception as e:
        logger.error(f"Error building sequence reference table: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()

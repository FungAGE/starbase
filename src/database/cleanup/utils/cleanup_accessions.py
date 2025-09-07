#!/usr/bin/env python3
"""
Optimized database cleaning script for accession tag correction.

This script performs the following tasks more efficiently:
1. Generate MD5 hashes for both normal and reverse complement sequences
2. Identify sequences with identical MD5 hashes for normal and reverse complement
3. Check for nested sequences
4. Update accession versions based on sequence relationships
"""

import hashlib
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
from sqlalchemy import text
from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.database.models.schema import Accessions, Ships, JoinedShips, Captains, StarshipFeatures, Gff
from src.utils.seq_utils import revcomp
from src.database.sql_manager import fetch_ships
from src.utils.classification_utils import generate_new_accession
import time
import sys

logger = get_logger(__name__)


def generate_sequence_hashes(sequence: str) -> Tuple[str, str]:
    """
    Generate MD5 hashes for both normal and reverse complement sequences.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        Tuple[str, str]: (normal_hash, reverse_complement_hash)
    """
    # Normalize sequence (uppercase, remove whitespace)
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    
    # Generate hash for normal sequence
    normal_hash = hashlib.md5(seq.encode()).hexdigest()
    
    # Generate hash for reverse complement
    rev_comp_seq = revcomp(seq)
    # Convert BioPython Seq object to string
    rev_comp_str = str(rev_comp_seq)
    rev_comp_hash = hashlib.md5(rev_comp_str.encode()).hexdigest()
    
    return normal_hash, rev_comp_hash


def fetch_all_sequences() -> pd.DataFrame:
    """
    Fetch all sequences from the database with their accession information.
    
    Returns:
        pd.DataFrame: DataFrame with accession_id, accession_tag, sequence, md5, and rev_comp_md5
    """
    session = StarbaseSession()
    
    query = """
    SELECT 
        a.id as accession_id,
        a.accession_tag,
        a.version_tag,
        s.sequence,
        s.md5,
        s.rev_comp_md5
    FROM accessions a
    INNER JOIN ships s ON a.id = s.accession_id
    WHERE s.sequence IS NOT NULL AND s.sequence != ''
    """
    
    try:
        df = pd.read_sql_query(query, session.bind)
        logger.info(f"Fetched {len(df)} sequences from database")
        return df
    except Exception as e:
        logger.error(f"Error fetching sequences: {str(e)}")
        raise
    finally:
        session.close()


def analyze_sequence_hashes(sequences_df: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Analyze sequences to find those with identical MD5 hashes for normal and reverse complement.
    Uses database-stored MD5 values for efficiency.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame with sequences, md5, and rev_comp_md5
        
    Returns:
        Dict[str, List[str]]: Dictionary mapping hash to list of accession tags
    """
    hash_groups = defaultdict(list)
    
    logger.info("Analyzing sequence hashes using database-stored MD5 values...")
    for idx, row in sequences_df.iterrows():
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{len(sequences_df)} sequences")
            
        accession_tag = row['accession_tag']
        normal_hash = row['md5']
        rev_comp_hash = row['rev_comp_md5']
        
        # Skip if either hash is missing
        if pd.isna(normal_hash) or pd.isna(rev_comp_hash):
            logger.warning(f"Skipping {accession_tag}: missing MD5 values")
            continue
        
        # Check if normal and reverse complement hashes are identical
        if normal_hash == rev_comp_hash:
            hash_groups[normal_hash].append(accession_tag)
            logger.debug(f"Found self-complementary sequence: {accession_tag}")
    
    # Filter to only groups with more than one accession
    duplicate_groups = {hash_val: accessions for hash_val, accessions in hash_groups.items() 
                       if len(accessions) > 1}
    
    logger.info(f"Found {len(duplicate_groups)} groups of sequences with identical normal/reverse complement hashes")
    
    return duplicate_groups


def find_reverse_complement_pairs(sequences_df: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Find sequences that are reverse complements of each other.
    Uses database-stored MD5 values for efficiency.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame with sequences, md5, and rev_comp_md5
        
    Returns:
        List[Tuple[str, str]]: List of (accession1, accession2) pairs that are reverse complements
    """
    rev_comp_pairs = []
    
    logger.info("Finding reverse complement pairs using database-stored MD5 values...")
    
    # Create dictionaries mapping hashes to lists of accessions (to handle duplicates)
    normal_hash_to_accessions = defaultdict(list)
    rev_comp_hash_to_accessions = defaultdict(list)
    
    for idx, row in sequences_df.iterrows():
        if idx % 100 == 0:
            logger.info(f"Processed {idx}/{len(sequences_df)} sequences for reverse complement analysis")
            
        accession_tag = row['accession_tag']
        normal_hash = row['md5']
        rev_comp_hash = row['rev_comp_md5']
        
        # Skip if either hash is missing
        if pd.isna(normal_hash) or pd.isna(rev_comp_hash):
            logger.warning(f"Skipping {accession_tag}: missing MD5 values")
            continue
        
        # Store mappings (append to lists to handle duplicates)
        normal_hash_to_accessions[normal_hash].append(accession_tag)
        rev_comp_hash_to_accessions[rev_comp_hash].append(accession_tag)
    
    # Find pairs where one sequence's normal hash matches another's reverse complement hash
    processed_pairs = set()
    
    for normal_hash, accs1 in normal_hash_to_accessions.items():
        if normal_hash in rev_comp_hash_to_accessions:
            accs2 = rev_comp_hash_to_accessions[normal_hash]
            
            # Create all possible pairs between the two groups
            for acc1 in accs1:
                for acc2 in accs2:
                    # Avoid self-pairs and duplicate pairs
                    if acc1 != acc2:
                        pair = tuple(sorted([acc1, acc2]))  # Sort to ensure consistent ordering
                        if pair not in processed_pairs:
                            rev_comp_pairs.append((acc1, acc2))
                            processed_pairs.add(pair)
                            logger.debug(f"Found reverse complement pair: {acc1} <-> {acc2}")
    
    logger.info(f"Found {len(rev_comp_pairs)} reverse complement pairs")
    return rev_comp_pairs


def find_nested_sequences(sequences_df: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Find sequences that are nested within other sequences using efficient batch alignment.
    
    Args:
        sequences_df (pd.DataFrame): DataFrame with sequences
        
    Returns:
        List[Tuple[str, str]]: List of (nested_accession, containing_accession) tuples
    """
    nested_pairs = []
    
    # Sort sequences by length (longest first) to optimize nested detection
    sequences_df_sorted = sequences_df.copy()
    sequences_df_sorted['seq_length'] = sequences_df_sorted['sequence'].str.len()
    sequences_df_sorted = sequences_df_sorted.sort_values('seq_length', ascending=False)
    
    # Create a list of (accession_tag, sequence) tuples, removing exact duplicates
    sequences = []
    seen_sequences = set()
    for _, row in sequences_df_sorted.iterrows():
        seq_upper = row['sequence'].upper()
        if seq_upper not in seen_sequences:
            sequences.append((row['accession_tag'], seq_upper))
            seen_sequences.add(seq_upper)
        else:
            logger.debug(f"Skipping duplicate sequence for accession: {row['accession_tag']}")
    
    logger.info(f"After removing exact duplicates: {len(sequences)} unique sequences")
    
    # Use efficient batch alignment-based detection
    logger.info("Using efficient batch alignment-based nested sequence detection")
    
    # Parameters for nested sequence detection
    min_coverage = 0.95  # 95% coverage threshold
    min_identity = 0.95  # 95% identity threshold
    
    logger.info(f"Parameters: min_coverage={min_coverage}, min_identity={min_identity}")
    
    start_time = time.time()
    
    # Process sequences in batches for efficiency
    batch_size = 50  # Process 50 sequences at a time
    total_batches = (len(sequences) + batch_size - 1) // batch_size
    
    for batch_idx in range(0, len(sequences), batch_size):
        batch_end = min(batch_idx + batch_size, len(sequences))
        batch_sequences = sequences[batch_idx:batch_end]
        
        logger.info(f"Processing batch {batch_idx//batch_size + 1}/{total_batches} "
                   f"(sequences {batch_idx+1}-{batch_end})")
        
        # For each sequence in the batch, check if it's nested in any longer sequence
        for acc_query, seq_query in batch_sequences:
            # Find all sequences that are longer than the query
            longer_sequences = [(acc, seq) for acc, seq in sequences 
                               if len(seq) > len(seq_query) and acc != acc_query]
            
            if not longer_sequences:
                continue
            
            # Use batch alignment to check if query is contained in any longer sequence
            containing_acc = check_sequence_containment_batch(seq_query, longer_sequences, 
                                                            min_coverage, min_identity)
            
            if containing_acc:
                nested_pairs.append((acc_query, containing_acc))
                logger.debug(f"Found nested sequence: {acc_query} is nested in {containing_acc}")
    
    elapsed = time.time() - start_time
    logger.info(f"Batch alignment-based nested sequence detection completed in {elapsed:.1f} seconds")
    logger.info(f"Found {len(nested_pairs)} nested sequence pairs")
    
    return nested_pairs


def check_sequence_containment_batch(query_seq: str, ref_sequences: List[Tuple[str, str]], 
                                   min_coverage: float = 0.95, min_identity: float = 0.95) -> Optional[str]:
    """
    Check if query sequence is contained within any of the reference sequences using batch alignment.
    
    Args:
        query_seq: Query sequence (shorter)
        ref_sequences: List of (accession, sequence) tuples for reference sequences (longer)
        min_coverage: Minimum coverage of query sequence required
        min_identity: Minimum sequence identity required
        
    Returns:
        str: Accession of the best containing sequence, or None if no containment found
    """
    import subprocess
    import tempfile
    
    query_len = len(query_seq)
    
    # Use minimap2 for batch alignment-based detection
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as query_file:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as ref_file:
            # Write query sequence
            query_file.write(f">query\n{query_seq}\n")
            query_file.flush()
            
            # Write all reference sequences
            for acc, seq in ref_sequences:
                ref_file.write(f">{acc}\n{seq}\n")
            ref_file.flush()
            
            # Run minimap2 alignment with short read preset for better sensitivity
            cmd = f"minimap2 -x sr -c -t 1 {ref_file.name} {query_file.name}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.debug(f"minimap2 failed for batch comparison: {result.stderr}")
                return None
            
            # Parse alignment results and find the best match
            best_match = None
            best_score = 0
            
            for line in result.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                
                ref_name = fields[5]  # Reference sequence name (accession)
                matches = int(fields[9])
                align_length = int(fields[10])
                
                coverage = align_length / query_len
                identity = matches / align_length if align_length > 0 else 0
                
                if coverage >= min_coverage and identity >= min_identity:
                    # Calculate a combined score (coverage * identity)
                    score = coverage * identity
                    
                    if score > best_score:
                        best_score = score
                        best_match = ref_name
                        logger.debug(f"Found containment: {ref_name} (coverage={coverage:.3f}, identity={identity:.3f})")
            
            return best_match


def check_sequence_containment(query_seq: str, ref_seq: str, min_coverage: float = 0.95, min_identity: float = 0.95) -> bool:
    """
    Check if query sequence is contained within reference sequence using alignment.
    
    Args:
        query_seq: Query sequence (shorter)
        ref_seq: Reference sequence (longer)
        min_coverage: Minimum coverage of query sequence required
        min_identity: Minimum sequence identity required
        
    Returns:
        bool: True if query is contained within reference
    """
    import subprocess
    import tempfile
    
    query_len = len(query_seq)
    
    # Use minimap2 for alignment-based detection
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as query_file:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as ref_file:
            # Write sequences to temporary files
            query_file.write(f">query\n{query_seq}\n")
            ref_file.write(f">reference\n{ref_seq}\n")
            query_file.flush()
            ref_file.flush()
            
            # Run minimap2 alignment with short read preset for better sensitivity
            cmd = f"minimap2 -c --cs -t 1 {ref_file.name} {query_file.name}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                logger.debug(f"minimap2 failed for sequence comparison: {result.stderr}")
                return False
            
            # Parse alignment results
            for line in result.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) < 10:
                    continue
                
                matches = int(fields[9])
                align_length = int(fields[10])
                
                coverage = align_length / query_len
                identity = matches / align_length if align_length > 0 else 0
                
                if coverage >= min_coverage and identity >= min_identity:
                    logger.debug(f"Found containment: coverage={coverage:.3f}, identity={identity:.3f}")
                    return True
    
    return False


def consolidate_accessions_by_hash(hash_groups: Dict[str, List[str]]) -> List[Dict]:
    """
    Consolidate accessions based on identical hash groups.
    
    Args:
        hash_groups (Dict[str, List[str]]): Groups of accessions with identical hashes
        
    Returns:
        List[Dict]: List of consolidation operations to perform
    """
    consolidations = []
    
    for hash_val, accessions in hash_groups.items():
        # Sort accessions to ensure consistent ordering
        accessions.sort()
        
        # Choose the accession with the lowest number as the primary
        primary_accession = accessions[0]
        secondary_accessions = accessions[1:]
        
        consolidation = {
            'type': 'hash_duplicate',
            'primary_accession': primary_accession,
            'secondary_accessions': secondary_accessions,
            'reason': f'Identical normal/reverse complement hash: {hash_val}'
        }
        
        consolidations.append(consolidation)
        logger.info(f"Will consolidate {len(secondary_accessions)} accessions under {primary_accession}")
    
    return consolidations


def consolidate_reverse_complement_pairs(rev_comp_pairs: List[Tuple[str, str]]) -> List[Dict]:
    """
    Consolidate reverse complement pairs under the accession with the lower number.
    
    Args:
        rev_comp_pairs (List[Tuple[str, str]]): List of (accession1, accession2) pairs
        
    Returns:
        List[Dict]: List of consolidation operations to perform
    """
    consolidations = []
    
    for acc1, acc2 in rev_comp_pairs:
        # Choose the accession with the lower number as primary
        if acc1 < acc2:
            primary_accession = acc1
            secondary_accession = acc2
        else:
            primary_accession = acc2
            secondary_accession = acc1
        
        consolidation = {
            'type': 'reverse_complement',
            'primary_accession': primary_accession,
            'secondary_accessions': [secondary_accession],
            'reason': f'Reverse complement of {primary_accession}'
        }
        
        consolidations.append(consolidation)
        logger.info(f"Will consolidate {secondary_accession} under {primary_accession} (reverse complement pair)")
    
    return consolidations


def consolidate_nested_sequences(nested_pairs: List[Tuple[str, str]]) -> List[Dict]:
    """
    Consolidate nested sequences under the containing sequence's accession.
    
    Args:
        nested_pairs (List[Tuple[str, str]]): List of (nested, containing) accession pairs
        
    Returns:
        List[Dict]: List of consolidation operations to perform
    """
    consolidations = []
    
    # Group nested sequences by their containing sequence
    nested_groups = defaultdict(list)
    for nested_acc, containing_acc in nested_pairs:
        nested_groups[containing_acc].append(nested_acc)
    
    for containing_acc, nested_accessions in nested_groups.items():
        consolidation = {
            'type': 'nested_sequence',
            'primary_accession': containing_acc,
            'secondary_accessions': nested_accessions,
            'reason': f'Nested within longer sequence {containing_acc}'
        }
        
        consolidations.append(consolidation)
        logger.info(f"Will consolidate {len(nested_accessions)} nested accessions under {containing_acc}")
    
    return consolidations


def apply_consolidations(consolidations: List[Dict], dry_run: bool = True) -> None:
    """
    Apply consolidation operations to the database.
    
    Args:
        consolidations (List[Dict]): List of consolidation operations
        dry_run (bool): If True, only log what would be done without making changes
    """
    session = StarbaseSession()
    
    try:
        for consolidation in consolidations:
            primary_acc = consolidation['primary_accession']
            secondary_accs = consolidation['secondary_accessions']
            reason = consolidation['reason']
            
            logger.info(f"Processing consolidation: {reason}")
            
            if dry_run:
                logger.info(f"DRY RUN: Would consolidate {secondary_accs} under {primary_acc}")
                continue
            
            # Get the primary accession record
            primary_record = session.query(Accessions).filter(
                Accessions.accession_tag == primary_acc
            ).first()
            
            if not primary_record:
                logger.error(f"Primary accession {primary_acc} not found")
                continue
            
            # Update all secondary accessions to point to the primary accession
            for secondary_acc in secondary_accs:
                # Get the secondary accession record
                secondary_record = session.query(Accessions).filter(
                    Accessions.accession_tag == secondary_acc
                ).first()
                
                if not secondary_record:
                    logger.warning(f"Secondary accession {secondary_acc} not found, skipping")
                    continue
                
                # Update ships table - find ships that reference the secondary accession
                ships_to_update = session.query(Ships).filter(
                    Ships.accession_id == secondary_record.id
                ).all()
                
                for ship in ships_to_update:
                    ship.accession_id = primary_record.id
                
                # Update joined_ships table - find joined_ships that reference the secondary accession
                joined_to_update = session.query(JoinedShips).filter(
                    JoinedShips.ship_id == secondary_record.id
                ).all()
                
                for joined in joined_to_update:
                    joined.ship_id = primary_record.id
                
                # Update captains table - find captains that reference the secondary accession
                captains_to_update = session.query(Captains).filter(
                    Captains.ship_id == secondary_record.id
                ).all()
                
                for captain in captains_to_update:
                    captain.ship_id = primary_record.id
                
                # Update starship_features table - find features that reference the secondary accession
                features_to_update = session.query(StarshipFeatures).filter(
                    StarshipFeatures.ship_id == secondary_record.id
                ).all()
                
                for feature in features_to_update:
                    feature.ship_id = primary_record.id
                
                # Update gff_entries table - find gff entries that reference the secondary accession
                gff_to_update = session.query(Gff).filter(
                    Gff.ship_id == secondary_record.id
                ).all()
                
                for gff in gff_to_update:
                    gff.ship_id = primary_record.id
                
                logger.info(f"Updated {len(ships_to_update)} ships, {len(joined_to_update)} joined_ships, {len(captains_to_update)} captains, {len(features_to_update)} features, and {len(gff_to_update)} gff entries for {secondary_acc}")
            
            # Commit the changes
            session.commit()
            logger.info(f"Successfully consolidated {len(secondary_accs)} accessions under {primary_acc}")
            
    except Exception as e:
        logger.error(f"Error applying consolidations: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def generate_consolidation_report(consolidations: List[Dict]) -> str:
    """
    Generate a report of all consolidation operations.
    
    Args:
        consolidations (List[Dict]): List of consolidation operations
        
    Returns:
        str: Formatted report
    """
    report = []
    report.append("=" * 80)
    report.append("DATABASE CLEANUP REPORT - ACCESSION CONSOLIDATION")
    report.append("=" * 80)
    report.append("")
    
    # Group by type
    hash_consolidations = [c for c in consolidations if c['type'] == 'hash_duplicate']
    rev_comp_consolidations = [c for c in consolidations if c['type'] == 'reverse_complement']
    nested_consolidations = [c for c in consolidations if c['type'] == 'nested_sequence']
    
    report.append(f"HASH DUPLICATE CONSOLIDATIONS: {len(hash_consolidations)}")
    report.append("-" * 50)
    for consolidation in hash_consolidations:
        report.append(f"Primary: {consolidation['primary_accession']}")
        report.append(f"Secondary: {', '.join(consolidation['secondary_accessions'])}")
        report.append(f"Reason: {consolidation['reason']}")
        report.append("")
    
    report.append(f"REVERSE COMPLEMENT CONSOLIDATIONS: {len(rev_comp_consolidations)}")
    report.append("-" * 50)
    for consolidation in rev_comp_consolidations:
        report.append(f"Primary: {consolidation['primary_accession']}")
        report.append(f"Secondary: {', '.join(consolidation['secondary_accessions'])}")
        report.append(f"Reason: {consolidation['reason']}")
        report.append("")
    
    report.append(f"NESTED SEQUENCE CONSOLIDATIONS: {len(nested_consolidations)}")
    report.append("-" * 50)
    for consolidation in nested_consolidations:
        report.append(f"Primary: {consolidation['primary_accession']}")
        report.append(f"Secondary: {', '.join(consolidation['secondary_accessions'])}")
        report.append(f"Reason: {consolidation['reason']}")
        report.append("")
    
    total_secondary = sum(len(c['secondary_accessions']) for c in consolidations)
    report.append(f"SUMMARY: {len(consolidations)} consolidation operations affecting {total_secondary} accessions")
    
    return "\n".join(report)

def add_new_accession(ship_id, new_accession):
    """
    Create a new accession for a ship entry and update related records.
    
    Args:
        ship_id (int): The ID of the ship that needs an accession
        new_accession (str): The new accession tag to assign
    """
    session = StarbaseSession()
    try:
        # Fetch the ship entry
        ship = session.query(Ships).filter_by(id=ship_id).one_or_none()
        if not ship:
            raise ValueError(f"No Ship found with id {ship_id}")

        # Check if ship already has an accession
        if ship.accession_id is not None:
            logger.warning(f"Ship {ship_id} already has accession_id {ship.accession_id}, skipping")
            return

        # Check if the accession tag already exists
        existing_accession = session.query(Accessions).filter_by(accession_tag=new_accession).first()
        if existing_accession:
            raise ValueError(f"Accession tag {new_accession} already exists with id {existing_accession.id}")

        # Try to find a joined_ships entry that might correspond to this ship
        # Look for joined_ships entries that reference this specific ship_id
        joined_ship = session.query(JoinedShips).filter_by(ship_id=ship_id).first()
        
        if joined_ship:
            ship_name = joined_ship.starshipID or f"SHIP_{ship_id}"
        else:
            # Look for orphaned joined_ships entries without ship_id
            orphaned_joined_ship = session.query(JoinedShips).filter(
                JoinedShips.ship_id.is_(None),
                JoinedShips.starshipID.isnot(None)
            ).first()
            
            if orphaned_joined_ship:
                ship_name = orphaned_joined_ship.starshipID
                # Link this orphaned entry to our ship
                orphaned_joined_ship.ship_id = ship_id
                joined_ship = orphaned_joined_ship
            else:
                ship_name = f"SHIP_{ship_id}"  # fallback name

        # Create new accession entry
        accession = Accessions(
            ship_name=ship_name,
            accession_tag=new_accession,
            version_tag="1",
        )
        session.add(accession)
        session.flush()  # assign accession.id before commit

        # Update ship to reference the new accession
        ship.accession_id = accession.id

        if joined_ship:
            logger.info(f"Linked joined_ships entry {joined_ship.id} (starshipID: {joined_ship.starshipID}) to ship {ship_id}")

        # Update related tables that reference ships.id
        # GFF and StarshipFeatures should already be properly linked to ship_id (ships.id)
        # We don't need to update them since they should already reference the correct ship_id
        
        # Log existing relationships for verification
        gffs = session.query(Gff).filter_by(ship_id=ship_id).all()
        if gffs:
            logger.info(f"Found {len(gffs)} existing GFF entries already linked to ship {ship_id}")

        features = session.query(StarshipFeatures).filter_by(ship_id=ship_id).all()
        if features:
            logger.info(f"Found {len(features)} existing StarshipFeatures entries already linked to ship {ship_id}")

        # Persist everything
        session.commit()
        logger.info(f"Successfully created accession {new_accession} for ship {ship_id}")

    except Exception as e:
        logger.error(f"Error adding new accession: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def main(dry_run: bool = True, output_report: str = None):
    """
    Main function to run the accession cleanup process.
    Step 1: Fetch all sequences
    Step 2: Analyze sequence hashes    
    Step 3: Find reverse complement pairs
    Step 4: Find nested sequences
    Step 5: Generate consolidation operations
    Step 6: Generate report
    Step 7: Apply consolidations (if not dry run)
    Step 8: assign new accessions if all the other checks pass (if not dry run)

    Args:
        dry_run (bool): If True, only analyze and report without making changes
        output_report (str): Path to save the consolidation report
    """
    logger.info("Starting database accession cleanup process")
    
    logger.info("Step 1: Fetching all sequences from database")
    sequences_df = fetch_all_sequences()
    
    logger.info("Step 2: Analyzing sequence hashes")
    hash_groups = analyze_sequence_hashes(sequences_df)
    
    logger.info("Step 3: Finding reverse complement pairs")
    rev_comp_pairs = find_reverse_complement_pairs(sequences_df)
    
    logger.info("Step 4: Finding nested sequences")
    nested_pairs = find_nested_sequences(sequences_df)
    
    logger.info("Step 5: Generating consolidation operations")
    hash_consolidations = consolidate_accessions_by_hash(hash_groups)
    rev_comp_consolidations = consolidate_reverse_complement_pairs(rev_comp_pairs)
    nested_consolidations = consolidate_nested_sequences(nested_pairs)    
    all_consolidations = hash_consolidations + rev_comp_consolidations + nested_consolidations

    logger.info("Step 6: Generating consolidation report")
    report = generate_consolidation_report(all_consolidations)
    
    if output_report:
        with open(output_report, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to {output_report}")
    else:
        print(report)
    
    if not dry_run and all_consolidations:
        logger.info("Step 7: Applying consolidations to database")
        apply_consolidations(all_consolidations, dry_run=False)
    elif dry_run:
        logger.info("Step 7: Skipping database changes (dry run mode)")
    
    if not dry_run:
        logger.info("Step 8: Assigning new accessions to ships that are still missing them")
        try:
            # Fetch all existing ships to get the latest accession numbers
            existing_ships = fetch_ships(curated=False, with_sequence=False)  # Get all ships for numbering
            
            # Find ships that don't have accessions
            session = StarbaseSession()
            ships_without_accessions = session.query(Ships).filter(
                Ships.accession_id.is_(None),
                Ships.sequence.isnot(None),
                Ships.sequence != ''
            ).all()
            session.close()
            
            logger.info(f"Found {len(ships_without_accessions)} ships without accessions")
            
            for ship in ships_without_accessions:
                # Generate a new accession number
                new_accession = generate_new_accession(existing_ships)
                
                # Add the new accession
                add_new_accession(ship.id, new_accession)
                
                # Update the existing_ships dataframe to include the new accession
                # This ensures subsequent accessions get the next number
                new_row = pd.DataFrame({
                    'accession_tag': [new_accession],
                    'accession_id': [ship.id]  # temporary placeholder
                })
                existing_ships = pd.concat([existing_ships, new_row], ignore_index=True)
                
                logger.info(f"Assigned accession {new_accession} to ship {ship.id}")
                
        except Exception as e:
            logger.error(f"Error in Step 8 (assigning new accessions): {str(e)}")
            raise

    logger.info("Optimized database accession cleanup process completed")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Optimized clean up database accession tags")
    parser.add_argument("--apply", action="store_true", 
                       help="Apply changes to database (default is dry run)")
    parser.add_argument("--report", type=str, 
                       help="Path to save consolidation report")
    
    args = parser.parse_args()
    
    main(dry_run=not args.apply, output_report=args.report)

#!/usr/bin/env python3
"""
Script to populate the rev_comp_md5 column in the ships table.

This script calculates the reverse complement MD5 hash for all sequences
in the ships table and updates the rev_comp_md5 column.
"""

import hashlib
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.database.models.schema import Ships
from src.utils.seq_utils import (revcomp, clean_sequence)

logger = get_logger(__name__)


def calculate_normal_md5(sequence: str) -> str:
    """
    Calculate MD5 hash of the normal sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: MD5 hash of the normal sequence
    """
    seq = clean_sequence(sequence)
    if seq is None:
        return None
    
    # Generate MD5 hash for normal sequence
    normal_hash = hashlib.md5(seq.encode()).hexdigest()
    
    return normal_hash


def calculate_rev_comp_md5(sequence: str) -> str:
    """
    Calculate MD5 hash of the reverse complement of a sequence.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: MD5 hash of the reverse complement
    """
    seq = clean_sequence(sequence)
    if seq is None:
        return None
    
    # Generate reverse complement
    rev_comp_seq = revcomp(seq)
    # Convert BioPython Seq object to string
    rev_comp_str = str(rev_comp_seq)
    
    # Generate MD5 hash
    rev_comp_hash = hashlib.md5(rev_comp_str.encode()).hexdigest()
    
    return rev_comp_hash


def validate_and_correct_md5(batch_size: int = 1000, dry_run: bool = False):
    """
    Validate and correct the existing md5 column for all ships in the database.
    
    Args:
        batch_size (int): Number of records to process in each batch
        dry_run (bool): If True, only show what would be updated without making changes
    """
    session = StarbaseSession()
    
    try:
        # Get total count of ships with sequences
        total_ships = session.query(Ships).filter(
            Ships.sequence.isnot(None),
            Ships.sequence != ''
        ).count()
        
        logger.info(f"Found {total_ships} ships with sequences to validate md5")
        
        if dry_run:
            logger.info("DRY RUN MODE: No changes will be made to the database")
        
        # Process in batches
        processed = 0
        corrected = 0
        errors = 0
        
        while processed < total_ships:
            # Get batch of ships
            ships_batch = session.query(Ships).filter(
                Ships.sequence.isnot(None),
                Ships.sequence != ''
            ).offset(processed).limit(batch_size).all()
            
            if not ships_batch:
                break
            
            logger.info(f"Processing batch: ships {processed+1}-{processed+len(ships_batch)} of {total_ships}")
            
            for ship in ships_batch:
                if ship.sequence:
                    # Calculate correct MD5
                    correct_md5 = calculate_normal_md5(ship.sequence)
                    
                    if correct_md5 is None:
                        logger.warning(f"Ship ID {ship.id}: Could not calculate MD5 (invalid sequence)")
                        errors += 1
                        continue
                    
                    # Check if the stored MD5 is different from the calculated one
                    if ship.md5 != correct_md5:
                        logger.info(f"Ship ID {ship.id}: MD5 mismatch - stored: {ship.md5}, calculated: {correct_md5}")
                        if not dry_run:
                            ship.md5 = correct_md5
                        corrected += 1
                    else:
                        logger.debug(f"Ship ID {ship.id}: MD5 is correct")
            
            processed += len(ships_batch)
            
            # Commit batch if not dry run
            if not dry_run:
                session.commit()
                logger.info(f"Committed batch. Total corrected: {corrected}, errors: {errors}")
            else:
                logger.info(f"Would correct {corrected} ships in this batch, {errors} errors")
        
        logger.info(f"Completed processing {processed} ships")
        logger.info(f"Total ships corrected: {corrected}")
        logger.info(f"Total errors: {errors}")
        
        if dry_run:
            logger.info("DRY RUN COMPLETED: No changes were made to the database")
        else:
            logger.info("Successfully validated and corrected md5 column")
            
    except Exception as e:
        logger.error(f"Error validating md5: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def populate_rev_comp_md5(batch_size: int = 1000, dry_run: bool = False):
    """
    Populate the rev_comp_md5 column for all ships in the database.
    
    Args:
        batch_size (int): Number of records to process in each batch
        dry_run (bool): If True, only show what would be updated without making changes
    """
    session = StarbaseSession()
    
    try:
        # Get total count of ships with sequences
        total_ships = session.query(Ships).filter(
            Ships.sequence.isnot(None),
            Ships.sequence != ''
        ).count()
        
        logger.info(f"Found {total_ships} ships with sequences to process")
        
        if dry_run:
            logger.info("DRY RUN MODE: No changes will be made to the database")
        
        # Process in batches
        processed = 0
        updated = 0
        
        while processed < total_ships:
            # Get batch of ships
            ships_batch = session.query(Ships).filter(
                Ships.sequence.isnot(None),
                Ships.sequence != ''
            ).offset(processed).limit(batch_size).all()
            
            if not ships_batch:
                break
            
            logger.info(f"Processing batch: ships {processed+1}-{processed+len(ships_batch)} of {total_ships}")
            
            for ship in ships_batch:
                if ship.sequence:
                    # Calculate reverse complement MD5
                    rev_comp_md5 = calculate_rev_comp_md5(ship.sequence)
                    
                    if rev_comp_md5 is None:
                        logger.warning(f"Ship ID {ship.id}: Could not calculate rev_comp_md5 (invalid sequence)")
                        continue
                    
                    # Check if the value is different from what's already stored
                    if ship.rev_comp_md5 != rev_comp_md5:
                        if not dry_run:
                            ship.rev_comp_md5 = rev_comp_md5
                        updated += 1
                        logger.debug(f"Updated ship ID {ship.id}: rev_comp_md5 = {rev_comp_md5}")
                    else:
                        logger.debug(f"Ship ID {ship.id} already has correct rev_comp_md5")
            
            processed += len(ships_batch)
            
            # Commit batch if not dry run
            if not dry_run:
                session.commit()
                logger.info(f"Committed batch. Total updated: {updated}")
            else:
                logger.info(f"Would update {updated} ships in this batch")
        
        logger.info(f"Completed processing {processed} ships")
        logger.info(f"Total ships updated: {updated}")
        
        if dry_run:
            logger.info("DRY RUN COMPLETED: No changes were made to the database")
        else:
            logger.info("Successfully populated rev_comp_md5 column")
            
    except Exception as e:
        logger.error(f"Error populating rev_comp_md5: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def verify_rev_comp_md5():
    """
    Verify that the rev_comp_md5 column has been populated correctly.
    """
    session = StarbaseSession()
    
    try:
        # Check how many ships have rev_comp_md5 populated
        total_ships = session.query(Ships).filter(
            Ships.sequence.isnot(None),
            Ships.sequence != ''
        ).count()
        
        populated_ships = session.query(Ships).filter(
            Ships.sequence.isnot(None),
            Ships.sequence != '',
            Ships.rev_comp_md5.isnot(None)
        ).count()
        
        logger.info(f"Total ships with sequences: {total_ships}")
        logger.info(f"Ships with rev_comp_md5 populated: {populated_ships}")
        logger.info(f"Percentage populated: {(populated_ships/total_ships)*100:.1f}%")
        
        # Sample verification - check a few random ships
        sample_ships = session.query(Ships).filter(
            Ships.sequence.isnot(None),
            Ships.sequence != '',
            Ships.rev_comp_md5.isnot(None)
        ).limit(5).all()
        
        logger.info("\nSample verification:")
        for ship in sample_ships:
            calculated_rev_comp_md5 = calculate_rev_comp_md5(ship.sequence)
            is_correct = ship.rev_comp_md5 == calculated_rev_comp_md5
            logger.info(f"Ship ID {ship.id}: {'✓' if is_correct else '✗'} "
                       f"stored={ship.rev_comp_md5}, calculated={calculated_rev_comp_md5}")
        
    except Exception as e:
        logger.error(f"Error verifying rev_comp_md5: {str(e)}")
        raise
    finally:
        session.close()


def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Validate and populate MD5 columns in ships table")
    parser.add_argument("--dry-run", action="store_true", 
                       help="Show what would be updated without making changes")
    parser.add_argument("--verify", action="store_true",
                       help="Verify that rev_comp_md5 column is populated correctly")
    parser.add_argument("--validate-md5", action="store_true",
                       help="Only validate and correct the existing md5 column")
    parser.add_argument("--populate-rev-comp", action="store_true",
                       help="Only populate the rev_comp_md5 column")
    parser.add_argument("--batch-size", type=int, default=1000,
                       help="Number of records to process in each batch (default: 1000)")
    
    args = parser.parse_args()
    
    if args.verify:
        logger.info("Verifying rev_comp_md5 column...")
        verify_rev_comp_md5()
    elif args.validate_md5:
        logger.info("Validating and correcting md5 column...")
        validate_and_correct_md5(batch_size=args.batch_size, dry_run=args.dry_run)
    elif args.populate_rev_comp:
        logger.info("Populating rev_comp_md5 column...")
        populate_rev_comp_md5(batch_size=args.batch_size, dry_run=args.dry_run)
    else:
        # Default: run both validation and population
        logger.info("Step 1: Validating and correcting md5 column...")
        validate_and_correct_md5(batch_size=args.batch_size, dry_run=args.dry_run)
        
        logger.info("Step 2: Populating rev_comp_md5 column...")
        populate_rev_comp_md5(batch_size=args.batch_size, dry_run=args.dry_run)


if __name__ == "__main__":
    main()

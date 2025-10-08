#!/usr/bin/env python3
"""
Script to fix ships table primary key issues.

This script:
1. Assigns proper sequential IDs to ships with NULL id (2,452 entries)
2. Links joined_ships entries to ships (747 entries)
3. Creates missing ships entries if needed (with empty sequence)

The ships table currently has no primary key constraint enforced, allowing NULL ids.
This script salvages the data by assigning proper IDs.

Usage:
    python fix_ships_primary_keys.py --dry-run    # Preview changes
    python fix_ships_primary_keys.py              # Apply changes
"""

import sys
import os
import argparse

# Add project root to path
PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from src.database.cleanup.utils.database_cleanup import fix_ships_primary_key_issues
from src.config.logging import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Fix ships table primary key issues (NULL ids and missing entries)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Analyze and report what would be done without making changes'
    )
    
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("SHIPS PRIMARY KEY FIX")
    logger.info("=" * 80)
    logger.info(f"Mode: {'DRY RUN' if args.dry_run else 'APPLYING CHANGES'}")
    logger.info("")
    
    try:
        # Run the fix
        report = fix_ships_primary_key_issues(dry_run=args.dry_run)
        
        # Display summary
        logger.info("")
        logger.info("=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        logger.info(f"NULL IDs fixed: {report['summary']['null_ids_fixed']}")
        logger.info(f"Ships created: {report['summary']['ships_created']}")
        logger.info(f"Joined_ships linked: {report['summary']['joined_ships_linked']}")
        logger.info(f"Next available ID: {report['summary']['next_available_id']}")
        logger.info("")
        logger.info(f"Recommendation: {report['summary']['recommendation']}")
        logger.info("")
        
        # Show some examples
        if report['null_ids_fixed']:
            logger.info("=" * 80)
            logger.info("Example NULL ID fixes (first 5):")
            for i, fix in enumerate(report['null_ids_fixed'][:5]):
                logger.info(f"{i+1}. Row {fix['rowid']}: assign id={fix['new_id']}, "
                          f"accession_id={fix['accession_id']}, "
                          f"has_sequence={fix['has_sequence']}")
        
        if report['ships_created']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Ships created: {len(report['ships_created'])}")
            for i, ship in enumerate(report['ships_created'][:5]):
                logger.info(f"{i+1}. Ship {ship['ship_id']}: "
                          f"accession_id={ship['accession_id']}, "
                          f"starshipID={ship['starshipID']}")
        
        if report['joined_ships_linked']:
            logger.info("")
            logger.info("=" * 80)
            logger.info(f"Joined_ships entries linked: {len(report['joined_ships_linked'])}")
            logger.info("Examples (first 5):")
            for i, link in enumerate(report['joined_ships_linked'][:5]):
                logger.info(f"{i+1}. JS {link['joined_ships_id']} ({link['starshipID']}): "
                          f"linked to ship {link['ship_id']}")
        
        # Final message
        logger.info("")
        logger.info("=" * 80)
        if args.dry_run:
            logger.info("DRY RUN COMPLETE - No changes were made")
            logger.info("Run without --dry-run to apply these changes")
        else:
            logger.info("CHANGES APPLIED SUCCESSFULLY")
            logger.info("")
            logger.info("Next steps:")
            logger.info("1. Update ship sequences from FASTA file")
            logger.info("2. Verify foreign key relationships are now correct")
        logger.info("=" * 80)
        
        return 0
        
    except Exception as e:
        logger.error(f"Error during ships primary key fix: {str(e)}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())


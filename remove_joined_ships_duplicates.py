#!/usr/bin/env python3
"""
Script to remove duplicate joined_ships entries with _ACC_XXXX or _SHIP_XXXX suffixes.

This script:
1. Identifies duplicate joined_ships entries by their starshipID suffixes
2. Retains the original entries (lowest ID)
3. Coalesces NULL values from originals with non-NULL values from duplicates
4. Deletes the duplicate entries

Usage:
    python remove_joined_ships_duplicates.py --dry-run     # See what would be done
    python remove_joined_ships_duplicates.py               # Apply changes
"""

import sys
import os
import argparse
import json
from datetime import datetime

# Add project root to path
PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from src.database.cleanup.utils.database_cleanup import remove_suffixed_joined_ships_duplicates
from src.config.logging import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Remove duplicate joined_ships entries with _ACC_XXXX or _SHIP_XXXX suffixes'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Analyze and report what would be done without making changes'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Path to save the detailed report (JSON format)'
    )
    
    args = parser.parse_args()
    
    logger.info("=" * 80)
    logger.info("JOINED_SHIPS DUPLICATE REMOVAL")
    logger.info("=" * 80)
    logger.info(f"Mode: {'DRY RUN' if args.dry_run else 'APPLYING CHANGES'}")
    logger.info("")
    
    try:
        # Run the duplicate removal
        report = remove_suffixed_joined_ships_duplicates(dry_run=args.dry_run)
        
        # Display summary
        logger.info("")
        logger.info("=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Duplicate groups found: {report['summary']['duplicate_groups_found']}")
        logger.info(f"Entries updated (coalesced): {report['summary']['entries_updated']}")
        logger.info(f"Entries deleted: {report['summary']['entries_deleted']}")
        logger.info(f"Total fields coalesced: {report['summary']['total_fields_coalesced']}")
        logger.info("")
        
        # Show some examples
        if report['duplicates_found']:
            logger.info("Example duplicate groups:")
            for i, group in enumerate(report['duplicates_found'][:5]):
                logger.info(f"\nGroup {i+1}:")
                logger.info(f"  Base starshipID: {group['base_starship_id']}")
                logger.info(f"  Base ID: {group['base_id']}")
                logger.info(f"  Duplicate count: {group['duplicate_count']}")
                logger.info(f"  Duplicate IDs: {group['duplicate_ids']}")
                if group['fields_coalesced']:
                    logger.info(f"  Fields coalesced: {', '.join(group['fields_coalesced'])}")
        
        # Show examples of updates
        if report['entries_updated']:
            logger.info("\n" + "=" * 80)
            logger.info("Example coalesced values:")
            for i, update in enumerate(report['entries_updated'][:3]):
                logger.info(f"\nUpdate {i+1}:")
                logger.info(f"  Entry ID: {update['joined_ships_id']}")
                logger.info(f"  StarshipID: {update['starshipID']}")
                logger.info(f"  Fields updated:")
                for field, data in update['fields_updated'].items():
                    logger.info(f"    {field}: NULL -> {data['new_value']} (from {data['source_starship_id']})")
        
        # Save report if requested
        if args.output:
            output_path = args.output
            with open(output_path, 'w') as f:
                json.dump(report, f, indent=2, default=str)
            logger.info(f"\nDetailed report saved to: {output_path}")
        
        # Final message
        logger.info("\n" + "=" * 80)
        if args.dry_run:
            logger.info("DRY RUN COMPLETE - No changes were made")
            logger.info("Run without --dry-run to apply these changes")
        else:
            logger.info("CHANGES APPLIED SUCCESSFULLY")
        logger.info("=" * 80)
        
        return 0
        
    except Exception as e:
        logger.error(f"Error during duplicate removal: {str(e)}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())


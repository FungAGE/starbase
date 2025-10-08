#!/usr/bin/env python3
"""
Standalone script to fill in missing family IDs in the joined_ships table.

This script uses a multi-step approach:
1. Inherit from existing classifications (if starshipID already has family elsewhere)
2. Fill based on shared navis (if navis has consensus family)
3. Fill based on shared haplotype (if haplotype has consensus family)
4. Classify using captain sequences (using HMMER on captain sequences)
"""

import sys
import os
import json
import argparse

# Add the project root to the Python path
_HERE = os.path.dirname(__file__)
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..', '..'))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from src.config.logging import get_logger
from src.database.cleanup.utils.database_cleanup import fill_missing_family_ids

logger = get_logger(__name__)


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Fill missing family IDs in joined_ships table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry run (default) - see what would be changed without making changes
  python fill_family_ids.py
  
  # Apply the changes
  python fill_family_ids.py --apply
  
  # Save report to file
  python fill_family_ids.py --report family_id_report.json
  
  # Apply changes and save report
  python fill_family_ids.py --apply --report family_id_report.json
        """
    )
    
    parser.add_argument(
        '--apply',
        action='store_true',
        help='Apply the changes (default is dry run only)'
    )
    
    parser.add_argument(
        '--report',
        type=str,
        help='Path to save the report as JSON file'
    )
    
    args = parser.parse_args()
    
    # Run the function
    dry_run = not args.apply
    
    if dry_run:
        logger.info("Running in DRY RUN mode - no changes will be made")
    else:
        logger.info("Running in APPLY mode - changes will be made to the database")
    
    logger.info("Starting family ID filling process...")
    
    try:
        report = fill_missing_family_ids(dry_run=dry_run)
        
        # Print summary
        print("\n" + "=" * 80)
        print("FAMILY ID FILLING REPORT")
        print("=" * 80)
        print(f"\nTotal missing family IDs: {report['summary']['total_missing']}")
        print(f"Inherited from existing: {report['summary']['inherited_from_existing']}")
        print(f"Filled from navis consensus: {report['summary']['filled_from_navis']}")
        print(f"Filled from haplotype consensus: {report['summary']['filled_from_haplotype']}")
        print(f"Filled from captain classification: {report['summary']['filled_from_captain']}")
        print(f"No classification found: {report['summary']['no_classification_found']}")
        print(f"Errors: {report['summary']['errors']}")
        print(f"\nTotal filled: {report['summary']['total_filled']}")
        
        if report['errors']:
            print(f"\nErrors (showing first 10):")
            for error in report['errors'][:10]:
                print(f"  - {error['starshipID']}: {error['error']}")
            if len(report['errors']) > 10:
                print(f"  ... and {len(report['errors']) - 10} more errors")
        
        # Show some examples of each method
        if report['inherited_from_existing']:
            print(f"\nInherited from existing (showing first 5):")
            for item in report['inherited_from_existing'][:5]:
                print(f"  - {item['starshipID']}: family_id={item['assigned_family_id']}")
        
        if report['filled_from_navis']:
            print(f"\nFilled from navis consensus (showing first 5):")
            for item in report['filled_from_navis'][:5]:
                print(f"  - {item['starshipID']}: family_id={item['assigned_family_id']} (navis: {item['navis_name']}, consensus: {item['consensus_ratio']:.1%})")
        
        if report['filled_from_haplotype']:
            print(f"\nFilled from haplotype consensus (showing first 5):")
            for item in report['filled_from_haplotype'][:5]:
                print(f"  - {item['starshipID']}: family_id={item['assigned_family_id']} (haplotype: {item['haplotype_name']}, consensus: {item['consensus_ratio']:.1%})")
        
        if report['filled_from_captain']:
            print(f"\nFilled from captain classification (showing first 5):")
            for item in report['filled_from_captain'][:5]:
                print(f"  - {item['starshipID']}: family_id={item['assigned_family_id']} (family: {item['family_name']})")
        
        print("\n" + "=" * 80)
        
        # Save report if requested
        if args.report:
            with open(args.report, 'w') as f:
                json.dump(report, f, indent=2, default=str)
            logger.info(f"Report saved to: {args.report}")
        
        if dry_run:
            print("\nNOTE: This was a DRY RUN. No changes were made to the database.")
            print("Use --apply to actually apply the changes.")
        else:
            print("\nChanges have been applied to the database.")
        
        # Return success
        return 0
        
    except Exception as e:
        logger.error(f"Error running family ID filling: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())


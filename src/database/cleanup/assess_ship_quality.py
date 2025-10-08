#!/usr/bin/env python3
"""
Standalone script to assess and tag ship quality.

This script demonstrates how to use the new ship_quality_tags system
to automatically assess data completeness and assign quality tags.

Usage:
    python assess_ship_quality.py --help
    python assess_ship_quality.py --assess-all --dry-run
    python assess_ship_quality.py --ship-id 123 --add-tag missing_captain
    python assess_ship_quality.py --summary
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from src.database.cleanup.utils.quality_tags import (
    add_quality_tag,
    remove_quality_tag,
    get_ship_quality_tags,
    find_ships_with_tags,
    get_quality_tag_summary,
    assess_ship_quality,
    bulk_assess_ship_quality,
    STANDARD_TAG_TYPES
)
from src.config.logging import get_logger

logger = get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Assess and manage ship quality tags",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get summary of all quality tags
  python assess_ship_quality.py --summary
  
  # Assess all ships (dry run)
  python assess_ship_quality.py --assess-all --dry-run
  
  # Actually assess all ships
  python assess_ship_quality.py --assess-all
  
  # Add a quality tag to a specific ship
  python assess_ship_quality.py --ship-id 123 --add-tag missing_captain
  
  # Remove a quality tag from a ship
  python assess_ship_quality.py --ship-id 123 --remove-tag missing_captain
  
  # Get all tags for a ship
  python assess_ship_quality.py --ship-id 123 --get-tags
  
  # Find ships with specific tags
  python assess_ship_quality.py --find-ships missing_captain incomplete
  
  # Assess a specific ship
  python assess_ship_quality.py --ship-id 123 --assess
        """
    )
    
    # Main actions (mutually exclusive)
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--summary", action="store_true",
                             help="Show summary of all quality tags")
    action_group.add_argument("--assess-all", action="store_true",
                             help="Assess quality for all ships")
    action_group.add_argument("--assess", action="store_true",
                             help="Assess quality for a specific ship (requires --ship-id)")
    action_group.add_argument("--add-tag", type=str,
                             help="Add a quality tag to a ship (requires --ship-id)")
    action_group.add_argument("--remove-tag", type=str,
                             help="Remove a quality tag from a ship (requires --ship-id)")
    action_group.add_argument("--get-tags", action="store_true",
                             help="Get all tags for a ship (requires --ship-id)")
    action_group.add_argument("--find-ships", nargs="+", metavar="TAG",
                             help="Find ships with specific quality tags")
    
    # Options
    parser.add_argument("--ship-id", type=int,
                       help="Specific ship ID to operate on")
    parser.add_argument("--dry-run", action="store_true",
                       help="Show what would be done without making changes")
    parser.add_argument("--created-by", default="manual",
                       help="Who is creating the tag (default: manual)")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.add_tag or args.remove_tag or args.get_tags or args.assess:
        if not args.ship_id:
            parser.error("--ship-id is required for this action")
    
    if args.add_tag and args.add_tag not in STANDARD_TAG_TYPES:
        logger.warning(f"Tag '{args.add_tag}' is not a standard tag type. Standard types: {STANDARD_TAG_TYPES}")
    
    try:
        if args.summary:
            show_summary()
        elif args.assess_all:
            assess_all_ships(args.dry_run)
        elif args.assess:
            assess_single_ship(args.ship_id, args.dry_run)
        elif args.add_tag:
            add_tag_to_ship(args.ship_id, args.add_tag, args.created_by, args.dry_run)
        elif args.remove_tag:
            remove_tag_from_ship(args.ship_id, args.remove_tag, args.dry_run)
        elif args.get_tags:
            get_tags_for_ship(args.ship_id)
        elif args.find_ships:
            find_ships_with_quality_tags(args.find_ships)
            
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        sys.exit(1)


def show_summary():
    """Show summary of all quality tags."""
    print("=== Quality Tag Summary ===")
    summary = get_quality_tag_summary()
    
    if not summary:
        print("No quality tags found in database.")
        return
    
    total_tags = sum(summary.values())
    print(f"Total tags: {total_tags}")
    print()
    
    for tag_type, count in sorted(summary.items()):
        print(f"  {tag_type:20} {count:6d}")


def assess_all_ships(dry_run: bool):
    """Assess quality for all ships."""
    print("=== Assessing All Ships ===")
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
    
    summary = bulk_assess_ship_quality(dry_run=dry_run)
    
    print(f"Total ships: {summary['total_ships']}")
    print(f"Ships assessed: {summary['ships_assessed']}")
    print(f"Complete ships: {summary['complete_ships']}")
    print(f"Incomplete ships: {summary['incomplete_ships']}")
    
    if summary['tags_added']:
        print("\nTags added:")
        for tag, count in summary['tags_added'].items():
            print(f"  {tag}: {count}")
    
    if summary['tags_removed']:
        print("\nTags removed:")
        for tag, count in summary['tags_removed'].items():
            print(f"  {tag}: {count}")


def assess_single_ship(ship_id: int, dry_run: bool):
    """Assess quality for a single ship."""
    print(f"=== Assessing Ship {ship_id} ===")
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
        return
    
    assessment = assess_ship_quality(ship_id)
    
    print(f"Ship ID: {assessment['ship_id']}")
    print(f"StarshipID: {assessment['starshipID']}")
    print(f"Is complete: {assessment['is_complete']}")
    
    if assessment['tags_added']:
        print(f"Tags added: {', '.join(assessment['tags_added'])}")
    
    if assessment['tags_removed']:
        print(f"Tags removed: {', '.join(assessment['tags_removed'])}")


def add_tag_to_ship(ship_id: int, tag_type: str, created_by: str, dry_run: bool):
    """Add a quality tag to a ship."""
    print(f"=== Adding Tag '{tag_type}' to Ship {ship_id} ===")
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
        return
    
    success = add_quality_tag(ship_id, tag_type, created_by=created_by)
    if success:
        print(f"Successfully added tag '{tag_type}' to ship {ship_id}")
    else:
        print(f"Tag '{tag_type}' already exists for ship {ship_id}")


def remove_tag_from_ship(ship_id: int, tag_type: str, dry_run: bool):
    """Remove a quality tag from a ship."""
    print(f"=== Removing Tag '{tag_type}' from Ship {ship_id} ===")
    if dry_run:
        print("DRY RUN MODE - No changes will be made")
        return
    
    success = remove_quality_tag(ship_id, tag_type)
    if success:
        print(f"Successfully removed tag '{tag_type}' from ship {ship_id}")
    else:
        print(f"Tag '{tag_type}' not found for ship {ship_id}")


def get_tags_for_ship(ship_id: int):
    """Get all quality tags for a ship."""
    print(f"=== Quality Tags for Ship {ship_id} ===")
    tags = get_ship_quality_tags(ship_id)
    
    if not tags:
        print("No quality tags found for this ship.")
        return
    
    for tag in tags:
        print(f"  {tag['tag_type']:20} {tag['created_by']:10} {tag['created_at']}")
        if tag['tag_value']:
            print(f"    Value: {tag['tag_value']}")


def find_ships_with_quality_tags(tag_types: list):
    """Find ships with specific quality tags."""
    print(f"=== Ships with Tags: {', '.join(tag_types)} ===")
    
    ship_ids = find_ships_with_tags(tag_types, match_all=False)
    
    if not ship_ids:
        print("No ships found with these tags.")
        return
    
    print(f"Found {len(ship_ids)} ships:")
    for ship_id in ship_ids[:20]:  # Show first 20
        print(f"  Ship ID: {ship_id}")
    
    if len(ship_ids) > 20:
        print(f"  ... and {len(ship_ids) - 20} more")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Runner script for comprehensive database cleanup.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from src.database.database_cleanup import (
    main,
    check_ships_accessions_joined_ships_relationships,
    analyze_table_relationships,
    fix_ships_accessions_joined_ships_relationships,
    identify_duplicate_joined_ships,
    cleanup_duplicate_joined_ships,
    analyze_ship_id_mislabeling,
    fix_ship_id_mislabeling,
    consolidate_duplicate_ships
)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run comprehensive database cleanup")
    parser.add_argument("--apply", action="store_true", 
                       help="Apply changes to database (default is dry run)")
    parser.add_argument("--apply-fixes", action="store_true",
                       help="Apply fixes to ships-accessions-joined_ships relationship issues")
    parser.add_argument("--report", type=str, 
                       help="Path to save cleanup report")
    parser.add_argument("--analyze-relationships", action="store_true",
                       help="Run only relationship analysis")
    parser.add_argument("--check-relationships", action="store_true",
                       help="Check for relationship issues")
    parser.add_argument("--skip-accession-cleanup", action="store_true",
                       help="Skip the slow accession cleanup step")
    parser.add_argument("--identify-duplicates", action="store_true",
                       help="Identify duplicate joined_ships entries")
    parser.add_argument("--cleanup-duplicates", action="store_true",
                       help="Clean up duplicate joined_ships entries")
    parser.add_argument("--analyze-mislabeling", action="store_true",
                       help="Analyze ship_id mislabeling in joined_ships table")
    parser.add_argument("--fix-mislabeling", action="store_true",
                       help="Fix ship_id mislabeling in joined_ships table")
    parser.add_argument("--consolidate-duplicate-ships", action="store_true",
                       help="Consolidate duplicate sequences in ships table while maintaining foreign key relationships")
    
    args = parser.parse_args()
    
    if args.analyze_relationships:
        # Run only the relationship analysis
        print("Running relationship analysis...")
        analysis = analyze_table_relationships()
        
        print("\n" + "=" * 80)
        print("TABLE RELATIONSHIP ANALYSIS")
        print("=" * 80)
        
        print("\nRELATIONSHIP COUNTS:")
        for key, value in analysis['relationship_counts'].items():
            print(f"  {key}: {value}")
        
        print("\nDATA DISTRIBUTION:")
        for key, value in analysis['data_distribution'].items():
            if 'pct' in key:
                print(f"  {key}: {value:.1f}%")
            else:
                print(f"  {key}: {value}")
        
        if analysis['potential_issues']:
            print("\nPOTENTIAL ISSUES:")
            for issue in analysis['potential_issues']:
                print(f"  - {issue}")
        
        if analysis['recommendations']:
            print("\nRECOMMENDATIONS:")
            for rec in analysis['recommendations']:
                print(f"  - {rec}")
    
    elif args.check_relationships:
        # Check for relationship issues
        print("Checking relationships...")
        issues = check_ships_accessions_joined_ships_relationships()
        
        print(f"\nFound {len(issues['orphaned_ships'])} orphaned ships")
        print(f"Found {len(issues['orphaned_accessions'])} orphaned accessions")
        print(f"Found {len(issues['orphaned_joined_ships'])} orphaned joined_ships")
        print(f"Found {len(issues['missing_sequence_links'])} missing sequence links")
        print(f"Found {len(issues['ships_missing_from_joined_ships'])} ships missing from joined_ships")
        print(f"Found {len(issues['inconsistent_joined_ships'])} inconsistent joined_ships")
        
        if issues['summary_stats']:
            stats = issues['summary_stats']
            print(f"\nTable counts - Ships: {stats.get('total_ships', 0)}, Accessions: {stats.get('total_accessions', 0)}, JoinedShips: {stats.get('total_joined_ships', 0)}")
    
    elif args.identify_duplicates:
        # Identify duplicate joined_ships entries
        print("Identifying duplicate joined_ships entries...")
        duplicates = identify_duplicate_joined_ships()
        
        print("\n" + "=" * 80)
        print("DUPLICATE JOINED_SHIPS ANALYSIS")
        print("=" * 80)
        
        print(f"\nShips with duplicate entries: {duplicates['summary']['total_duplicate_ship_ids']}")
        print(f"Duplicate starshipIDs: {duplicates['summary']['total_duplicate_starship_ids']}")
        print(f"Auto-created entries: {duplicates['summary']['total_auto_created']}")
        
        if duplicates['duplicate_ship_ids']:
            print("\nShips with duplicate entries:")
            for dup in duplicates['duplicate_ship_ids'][:10]:  # Show first 10
                print(f"  Ship {dup['ship_id']}: {dup['count']} entries")
                for entry in dup['entries']:
                    print(f"    - ID: {entry['id']}, StarshipID: {entry['starshipID']}, Source: {entry['source']}")
            if len(duplicates['duplicate_ship_ids']) > 10:
                print(f"  ... and {len(duplicates['duplicate_ship_ids']) - 10} more")
    
    elif args.cleanup_duplicates:
        # Clean up duplicate joined_ships entries
        print("Cleaning up duplicate joined_ships entries...")
        cleanup_report = cleanup_duplicate_joined_ships(dry_run=not args.apply)
        
        print(f"\nCleaned up duplicates for {cleanup_report['summary']['ships_with_duplicates']} ships")
        print(f"Total entries removed: {cleanup_report['summary']['total_entries_removed']}")
        
        if cleanup_report['entries_removed']:
            print("\nEntries removed:")
            for entry in cleanup_report['entries_removed'][:10]:  # Show first 10
                print(f"  - ID: {entry['entry_id']}, Ship: {entry['ship_id']}, StarshipID: {entry['starshipID']}, Source: {entry['source']}")
            if len(cleanup_report['entries_removed']) > 10:
                print(f"  ... and {len(cleanup_report['entries_removed']) - 10} more")
    
    elif args.analyze_mislabeling:
        # Analyze ship_id mislabeling
        print("Analyzing ship_id mislabeling in joined_ships table...")
        mislabeling_analysis = analyze_ship_id_mislabeling()
        
        print("\n" + "=" * 80)
        print("SHIP_ID MISLABELING ANALYSIS")
        print("=" * 80)
        
        indicators = mislabeling_analysis['mislabeling_indicators']
        print(f"\nTotal joined_ships entries: {indicators['total_joined_entries']}")
        print(f"Ship_id values exist in ships table: {indicators['ship_ids_exist_in_ships']}")
        print(f"Ship_id values exist in accessions table: {indicators['ship_ids_exist_in_accessions']}")
        print(f"Mislabeled entries: {indicators['mislabeled_entries']}")
        print(f"Ships missing from joined_ships: {indicators['ships_without_joined']}")
        print(f"Accessions incorrectly referenced: {indicators['accessions_incorrectly_referenced']}")
        
        if mislabeling_analysis['potential_corrections']:
            print(f"\nPotential corrections needed: {len(mislabeling_analysis['potential_corrections'])}")
            print("First few corrections:")
            for correction in mislabeling_analysis['potential_corrections'][:5]:
                print(f"  - {correction['action']}")
            if len(mislabeling_analysis['potential_corrections']) > 5:
                print(f"  ... and {len(mislabeling_analysis['potential_corrections']) - 5} more")
        
        print(f"\nRecommendation: {mislabeling_analysis['summary']['recommendation']}")
    
    elif args.fix_mislabeling:
        # Fix ship_id mislabeling
        print("Fixing ship_id mislabeling in joined_ships table...")
        fix_report = fix_ship_id_mislabeling(dry_run=not args.apply)
        
        print(f"\nFixed {fix_report['summary']['total_fixed']} mislabeled entries")
        print(f"Skipped {fix_report['summary']['total_skipped']} correctly labeled entries")
        print(f"Recommendation: {fix_report['summary']['recommendation']}")
        
        if fix_report['entries_fixed']:
            print("\nEntries fixed:")
            for entry in fix_report['entries_fixed'][:10]:  # Show first 10
                print(f"  - {entry['action']}")
            if len(fix_report['entries_fixed']) > 10:
                print(f"  ... and {len(fix_report['entries_fixed']) - 10} more")
    
    elif args.consolidate_duplicate_ships:
        # Consolidate duplicate sequences in ships table
        print("Consolidating duplicate sequences in ships table...")
        consolidation_report = consolidate_duplicate_ships(dry_run=not args.apply)
        
        print("\n" + "=" * 80)
        print("DUPLICATE SHIPS CONSOLIDATION REPORT")
        print("=" * 80)
        
        print(f"\nTotal duplicate groups consolidated: {consolidation_report['summary']['total_groups_consolidated']}")
        print(f"Total joined_ships entries updated: {consolidation_report['summary']['total_joined_ships_updated']}")
        print(f"Total duplicate ships removed: {consolidation_report['summary']['total_duplicate_ships_removed']}")
        print(f"Recommendation: {consolidation_report['summary']['recommendation']}")
        
        if consolidation_report['duplicate_groups_found']:
            print("\nDuplicate groups found (same canonical MD5, different ship IDs):")
            for group in consolidation_report['duplicate_groups_found'][:5]:
                print(f"  Canonical MD5: {group['canonical_md5']}")
                print(f"    Primary ship ID: {group['primary_ship_id']}")
                print(f"    Duplicate ship IDs: {', '.join(map(str, group['duplicate_ship_ids']))}")
                print(f"    Total ships in group: {group['count']}")
                print(f"    Original MD5s: {', '.join(group['original_md5s'])}")
            if len(consolidation_report['duplicate_groups_found']) > 5:
                print(f"  ... and {len(consolidation_report['duplicate_groups_found']) - 5} more")
        
        if consolidation_report['ships_consolidated']:
            print("\nShips consolidated:")
            for ship in consolidation_report['ships_consolidated'][:5]:
                print(f"  - {ship['action']}")
            if len(consolidation_report['ships_consolidated']) > 5:
                print(f"  ... and {len(consolidation_report['ships_consolidated']) - 5} more")
    
    else:
        # Run full cleanup process
        main(dry_run=not args.apply, output_report=args.report, apply_fixes=args.apply_fixes, skip_accession_cleanup=args.skip_accession_cleanup)

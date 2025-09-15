#!/usr/bin/env python3
"""
Runner script for comprehensive database cleanup.
"""

import sys
import os

import pandas as pd
from sqlalchemy import text, func
from typing import Dict
from datetime import datetime

# Add the project root to the Python path (for direct script execution)
_HERE = os.path.dirname(__file__)
_PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, '..', '..', '..'))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

try:
    from ...config.logging import get_logger
    from .utils.cleanup_accessions import main as run_accession_cleanup
    from .taxonomy_consolidator import TaxonomyConsolidator

    from .utils.database_cleanup import *

except ImportError:
    # Fallback for when running the script directly
    import sys
    import os
    sys.path.append(_PROJECT_ROOT)
    
    from src.config.logging import get_logger
    from src.database.cleanup.utils.cleanup_accessions import main as run_accession_cleanup
    from src.database.cleanup.taxonomy_consolidator import TaxonomyConsolidator

    from src.database.cleanup.utils.database_cleanup import *

logger = get_logger(__name__)


def main(dry_run: bool = True, output_report: str = None, apply_fixes: bool = False, skip_accession_cleanup: bool = False, no_record_issues: bool = False):
    """
    Main function to run comprehensive database cleanup.
    Step 1: Run accession cleanup (optional, can be skipped)
    - checks for matches to reverse complemented sequences
    - checks for nested sequences
    - aggregates sequences by accession
    Step 2: Check genome table
    Step 3: Check taxonomic table
    Step 4: Check genomic features
    Step 5: Check foreign keys
    Step 6: Check ships-accessions-joined_ships relationships
    Step 7: Check schema violations
    Step 8: Validate and fix ship_id relationships using sequence_reference table
    Step 9: Perform comprehensive foreign key validation
    Step 10: Generate comprehensive report
    Step 11: Record issues and fixes to cleanup_issues table (optional)
    Step 12: Apply fixes (if requested)

    Args:
        dry_run (bool): If True, only analyze and report without making changes
        output_report (str): Path to save the cleanup report
        apply_fixes (bool): If True, apply fixes to the identified issues
        skip_accession_cleanup (bool): If True, skip the slow accession cleanup step
        no_record_issues (bool): If True, skip recording issues and fixes to cleanup_issues table
    """
    logger.info("Starting comprehensive database cleanup process")
    
    all_issues = {}
    
    if skip_accession_cleanup:
        logger.info("Step 1: Skipping accession cleanup (--skip-accession-cleanup specified)")
    else:
        logger.info("Step 1: Running accession cleanup...")
        run_accession_cleanup(dry_run=dry_run)
    
    logger.info("Step 2: Checking genome table...")
    all_issues['genome'] = check_genome_table()
    
    logger.info("Step 3: Checking taxonomic table...")
    all_issues['taxonomy'] = check_taxonomic_table()
    
    logger.info("Step 4: Checking genomic features...")
    all_issues['features'] = check_genomic_features()
    
    logger.info("Step 5: Checking foreign key consistency...")
    all_issues['foreign_keys'] = check_foreign_keys()
    
    logger.info("Step 6: Checking ships-accessions-joined_ships relationships...")
    all_issues['ships_accessions_joined_ships'] = check_ships_accessions_joined_ships_relationships()
    
    logger.info("Step 7: Checking schema violations...")
    all_issues['schema'] = check_schema_violations()
    
    logger.info("Step 8: Validating ship_id relationships using sequence_reference table...")
    all_issues['ship_id_relationships'] = analyze_ship_id_relationships()

    logger.info("Step 9: Performing comprehensive foreign key validation...")
    all_issues['foreign_key_validation'] = validate_all_foreign_keys()

    logger.info("Step 10: Generating comprehensive report...")
    report = generate_cleanup_report(all_issues)
    
    if output_report:
        with open(output_report, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to {output_report}")
    else:
        print(report)

    # Record issues into tracking table (table creation happens even in dry-run)
    if not no_record_issues:
        try:
            logger.info("Step 11: Recording issues into cleanup_issues table...")
            create_cleanup_issues_table()
            record_summary = record_cleanup_issues(all_issues, source='pipeline', dry_run=dry_run)
            logger.info(f"Recorded issues summary: inserted={record_summary['inserted']} skipped={record_summary['skipped']}")
        except Exception as e:
            logger.error(f"Failed to record cleanup issues: {str(e)}")
    else:
        logger.info("Step 11: Skipping issue recording (--no-record-issues specified)")

    if apply_fixes:
        logger.info("Step 12: Applying fixes to identified issues...")

        # Apply general relationship fixes
        fixes_report = fix_ships_accessions_joined_ships_relationships(dry_run=False)

        # Apply ship_id relationship fixes
        ship_id_fixes = fix_ship_id_relationships(dry_run=False)

        # Apply joined_ships.ship_id foreign key fix
        joined_ships_fixes = fix_joined_ships_ship_id_foreign_key(dry_run=False)

        logger.info("Fixes applied successfully")

        # Record fix results into cleanup_issues table
        try:
            logger.info("Recording fix results into cleanup_issues table...")

            # Collect all fix results
            all_fixes = {
                'ships_accessions_joined_ships': fixes_report,
                'ship_id_relationships': ship_id_fixes,
                'joined_ships_ship_id': joined_ships_fixes
            }

            fixes_record_summary = record_cleanup_fixes(all_fixes, source='pipeline', dry_run=False)
            logger.info(f"Recorded fixes summary: fixes={fixes_record_summary['fixes_recorded']} skipped={fixes_record_summary['skipped']}")
        except Exception as e:
            logger.error(f"Failed to record cleanup fixes: {str(e)}")

        # Add fixes report to output
        if output_report:
            with open(output_report, 'a') as f:
                f.write("\n\n" + "=" * 80 + "\n")
                f.write("FIXES APPLIED\n")
                f.write("=" * 80 + "\n")
                f.write("GENERAL RELATIONSHIP FIXES:\n")
                f.write(f"Ships fixed: {len(fixes_report['ships_fixed'])}\n")
                f.write(f"Accessions fixed: {len(fixes_report['accessions_fixed'])}\n")
                f.write(f"JoinedShips fixed: {len(fixes_report['joined_ships_fixed'])}\n")
                f.write(f"New joined_ships entries created: {len(fixes_report['new_joined_ships_created'])}\n")
                f.write("\nSHIP_ID RELATIONSHIP FIXES:\n")
                f.write(f"Missing ship_ids fixed: {len(ship_id_fixes['missing_ship_ids_fixed'])}\n")
                f.write(f"Mismatched ship_ids fixed: {len(ship_id_fixes['mismatched_ship_ids_fixed'])}\n")
                f.write(f"Orphaned ships handled: {len(ship_id_fixes['orphaned_ships_handled'])}\n")
                f.write(f"Duplicate ship_ids resolved: {len(ship_id_fixes['duplicate_ship_ids_resolved'])}\n")
                f.write(f"Warnings: {len(ship_id_fixes['warnings'])}\n")
                f.write("\nJOINED_SHIPS SHIP_ID FOREIGN KEY FIXES:\n")
                f.write(f"Ship IDs fixed: {joined_ships_fixes['summary']['ship_ids_fixed']}\n")
                f.write(f"Invalid references found: {joined_ships_fixes['summary']['invalid_references']}\n")
                f.write(f"Accessions without ships: {joined_ships_fixes['summary']['no_ship_for_accession']}\n")
    
    logger.info("Comprehensive database cleanup process completed")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Comprehensive database cleanup for starbase")
    parser.add_argument("--apply", action="store_true", 
                       help="Apply changes to database (default is dry run)")
    parser.add_argument("--apply-fixes", action="store_true",
                       help="Apply fixes to ships-accessions-joined_ships relationship issues")
    parser.add_argument("--report", type=str, 
                       help="Path to save cleanup report")
    parser.add_argument("--analyze-relationships", action="store_true",
                       help="Run detailed analysis of table relationships only")
    parser.add_argument("--check-relationships", action="store_true",
                       help="Check for relationship issues")
    parser.add_argument("--skip-accession-cleanup", action="store_true",
                       help="Skip the slow accession cleanup step")
    parser.add_argument("--identify-duplicates", action="store_true",
                       help="Identify duplicate joined_ships entries")
    parser.add_argument("--cleanup-duplicates", action="store_true",
                       help="Clean up duplicate joined_ships entries")
    parser.add_argument("--fix-placeholder-starship-ids", action="store_true",
                       help="Fix placeholder starshipID entries like 'SHIP_XXX' with proper accession-based IDs")
    parser.add_argument("--analyze-mislabeling", action="store_true",
                       help="Analyze ship_id mislabeling in joined_ships table")
    parser.add_argument("--fix-mislabeling", action="store_true",
                       help="Fix ship_id mislabeling in joined_ships table")
    parser.add_argument("--consolidate-duplicate-ships", action="store_true",
                       help="Consolidate duplicate sequences in ships table while maintaining foreign key relationships")
    parser.add_argument("--validate-ship-relationships", action="store_true",
                       help="Run ship_id relationship validation only")
    parser.add_argument("--fix-ship-relationships", action="store_true",
                       help="Apply ship_id relationship fixes")
    parser.add_argument("--analyze-missing-genome-info", action="store_true",
                       help="Analyze missing genome information")
    parser.add_argument("--fix-missing-genome-taxonomy-from-joined", action="store_true",
                       help="Fill missing Genome.taxonomy_id from joined_ships.tax_id (mode)")
    parser.add_argument("--fix-missing-genome-taxonomy-via-ncbi", action="store_true",
                       help="Fill missing Genome.taxonomy_id by querying NCBI Taxonomy (requires --ome-map)")
    parser.add_argument("--ncbi-email", type=str,
                       help="Contact email for NCBI E-utilities")
    parser.add_argument("--ncbi-api-key", type=str,
                       help="Optional NCBI API key for higher rate limits")
    parser.add_argument("--fix-missing-genome-taxonomy-via-taxdump", action="store_true",
                       help="Fill missing Genome.taxonomy_id using NCBI taxdump names.dmp")
    parser.add_argument("--names-dmp", type=str,
                       help="Path to NCBI taxdump names.dmp file")
    parser.add_argument("--include-synonyms", action="store_true",
                       help="Include synonyms/equivalent/common names from names.dmp")
    parser.add_argument("--fix-missing-genome-info-via-ome-map", action="store_true",
                       help="Reconcile existing genome taxonomy_ids to match OME mapping TSV")
    parser.add_argument("--ome-map", type=str,
                       help="Path to OME taxonomy mapping TSV (columns: ome, genus, species_epithet, strain, json_tax, ...) ")
    parser.add_argument("--lookup-assemblies-from-contigs", action="store_true",
                       help="Query NCBI to map contigIDs to assembly accessions and record results to cleanup_issues")
    parser.add_argument("--lookup-limit", type=int, default=0,
                       help="Limit number of unique contigIDs to query (0 = no limit)")

    # Foreign key management options
    parser.add_argument("--enable-foreign-keys", action="store_true",
                       help="Enable foreign key constraints in the database")
    parser.add_argument("--disable-foreign-keys", action="store_true",
                       help="Disable foreign key constraints in the database")
    parser.add_argument("--check-foreign-key-enforcement", action="store_true",
                       help="Check if foreign key constraints are actually being enforced")
    parser.add_argument("--fix-joined-ships-ship-id", action="store_true",
                       help="Fix joined_ships.ship_id foreign key to point to ships.id instead of accessions.id")
    parser.add_argument("--validate-all-foreign-keys", action="store_true",
                       help="Perform comprehensive validation of all foreign key relationships")
    parser.add_argument("--no-record-issues", action="store_true",
                       help="Skip recording issues and fixes to cleanup_issues table")
    parser.add_argument("--show-issues-summary", action="store_true",
                       help="Show summary of issues in cleanup_issues table")
    parser.add_argument("--identify-taxonomy-duplicates", action="store_true",
                       help="Identify duplicate taxonomy entries based on taxID and consistent information")
    parser.add_argument("--fix-missing-tax-id-via-ome", action="store_true",
                       help="Fix missing tax_id in joined_ships using ome code consistency")

    # Consolidated taxonomy options
    parser.add_argument("--consolidate-taxonomy", action="store_true",
                       help="Run complete taxonomy consolidation pipeline")
    parser.add_argument("--clean-taxonomy-data", action="store_true", 
                       help="Clean taxonomy data (whitespace, duplicates, etc.)")
    parser.add_argument("--resolve-taxonomy-ids", action="store_true",
                       help="Resolve missing taxonomy IDs from multiple sources")
    parser.add_argument("--backfill-taxonomy-hierarchy", action="store_true",
                       help="Backfill missing taxonomy hierarchy fields")
    parser.add_argument("--validate-taxonomy-consistency", action="store_true",
                       help="Validate taxonomy consistency and flag issues")
    parser.add_argument("--consolidate-taxonomy-duplicates", action="store_true",
                       help="Consolidate duplicate taxonomy entries and update foreign keys")

    args = parser.parse_args()
    
    if args.analyze_relationships:
        # Run only the relationship analysis
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

    elif args.validate_ship_relationships:
        # Run only ship_id relationship validation
        print("🔍 Validating ship_id relationships...")
        analysis = analyze_ship_id_relationships()
        
        print("\n" + "=" * 80)
        print("SHIP_ID RELATIONSHIP ANALYSIS")
        print("=" * 80)
        
        # Sequence reference stats
        print(f"\n📊 SEQUENCE REFERENCE TABLE:")
        stats = analysis['sequence_reference_stats']
        print(f"   Total entries: {stats['total_entries']}")
        print(f"   Unique sequences: {stats['unique_sequences']}")
        print(f"   Duplicate sequences: {stats['duplicate_sequences']}")
        
        # Ships table analysis
        print(f"\n🚢 SHIPS TABLE:")
        ships = analysis['ships_table_analysis']
        print(f"   Total ships: {ships['total_ships']}")
        print(f"   With sequences: {ships['ships_with_sequences']}")
        print(f"   With MD5: {ships['ships_with_md5']}")
        print(f"   Without sequences: {ships['ships_without_sequences']}")
        
        # Joined ships analysis
        print(f"\n🔗 JOINED_SHIPS TABLE:")
        joined = analysis['joined_ships_analysis']
        print(f"   Total entries: {joined['total_joined_ships']}")
        print(f"   With ship_id: {joined['with_ship_id']}")
        print(f"   Without ship_id: {joined['without_ship_id']}")
        print(f"   Duplicate ship_id cases: {joined['duplicate_ship_ids']}")
        
        if joined['duplicate_details']:
            print(f"   Duplicate details:")
            for dup in joined['duplicate_details'][:3]:  # Show first 3
                print(f"     ship_id {dup['ship_id']}: {dup['count']} entries ({', '.join(dup['starshipIDs'][:3])})")
        
        # Relationship issues
        print(f"\n⚠️  RELATIONSHIP ISSUES:")
        issues = analysis['relationship_issues']
        if not issues:
            print("   ✅ No issues found!")
        else:
            for issue in issues:
                print(f"   {issue['type']}: {issue['count']} - {issue['description']}")
                if issue['examples']:
                    print(f"     Examples: {issue['examples']}")
        
        # Recommendations
        print(f"\n💡 RECOMMENDATIONS:")
        recommendations = analysis['recommendations']
        if not recommendations:
            print("   ✅ No fixes needed!")
        else:
            for i, rec in enumerate(recommendations, 1):
                print(f"   {i}. {rec}")
        
        if args.fix_ship_relationships:
            print(f"\n🔧 Applying ship_id relationship fixes...")
            fixes = fix_ship_id_relationships(dry_run=False)
            
            print(f"\n📈 FIXES APPLIED:")
            print(f"   Missing ship_ids fixed: {len(fixes['missing_ship_ids_fixed'])}")
            print(f"   Mismatched ship_ids fixed: {len(fixes['mismatched_ship_ids_fixed'])}")
            print(f"   Orphaned ships handled: {len(fixes['orphaned_ships_handled'])}")
            print(f"   Duplicate ship_ids resolved: {len(fixes['duplicate_ship_ids_resolved'])}")
            print(f"   Warnings: {len(fixes['warnings'])}")

    elif args.fix_ship_relationships:
        # Apply ship_id relationship fixes as standalone operation
        print("🔧 Applying ship_id relationship fixes...")
        fixes = fix_ship_id_relationships(dry_run=not args.apply)
        
        print(f"\n📈 FIXES APPLIED:")
        print(f"   Missing ship_ids fixed: {len(fixes['missing_ship_ids_fixed'])}")
        print(f"   Mismatched ship_ids fixed: {len(fixes['mismatched_ship_ids_fixed'])}")
        print(f"   Orphaned ships handled: {len(fixes['orphaned_ships_handled'])}")
        print(f"   Duplicate ship_ids resolved: {len(fixes['duplicate_ship_ids_resolved'])}")
        print(f"   Warnings: {len(fixes['warnings'])}")

        if fixes['missing_ship_ids_fixed']:
            print("\n✅ MISSING SHIP_IDs FIXED (first 10):")
            for fix in fixes['missing_ship_ids_fixed'][:10]:
                print(f"   {fix['starshipID']}: ship_id={fix['ship_id']}")
        
        if fixes['mismatched_ship_ids_fixed']:
            print("\n🔄 MISMATCHED SHIP_IDs FIXED (first 10):")
            for fix in fixes['mismatched_ship_ids_fixed'][:10]:
                print(f"   {fix['starshipID']}: {fix['old_ship_id']} -> {fix['new_ship_id']}")
        
        if fixes['orphaned_ships_handled']:
            print("\n🔗 ORPHANED SHIPS HANDLED (first 10):")
            for fix in fixes['orphaned_ships_handled'][:10]:
                print(f"   Ship {fix['ship_id']}: linked to joined_ships {fix['joined_id']}")
        
        if fixes['warnings']:
            print("\n⚠️  WARNINGS (first 5):")
            for warning in fixes['warnings'][:5]:
                print(f"   {warning}")

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
    
    elif args.fix_placeholder_starship_ids:
        # Fix placeholder starshipID entries
        print("Fixing placeholder starshipID entries...")
        fix_report = fix_placeholder_starship_ids(dry_run=not args.apply)
        
        print(f"\nFound {fix_report['summary']['total_placeholder_entries']} placeholder entries")
        print(f"Fixed: {fix_report['summary']['fixed_count']}")
        print(f"No accession found: {fix_report['summary']['no_accession_count']}")
        print(f"Warnings: {fix_report['summary']['warnings_count']}")
        
        if fix_report['placeholder_ids_fixed']:
            print("\nFixed placeholder entries:")
            for fix in fix_report['placeholder_ids_fixed'][:10]:  # Show first 10
                print(f"  - {fix['old_starshipID']} -> {fix['new_starshipID']} (Joined ID: {fix['joined_id']}, Ship: {fix['ship_id']})")
            if len(fix_report['placeholder_ids_fixed']) > 10:
                print(f"  ... and {len(fix_report['placeholder_ids_fixed']) - 10} more")
        
        if fix_report['no_accession_found']:
            print("\nEntries without accession information:")
            for entry in fix_report['no_accession_found'][:5]:  # Show first 5
                print(f"  - {entry['starshipID']} (Joined ID: {entry['joined_id']}, Ship: {entry['ship_id']})")
            if len(fix_report['no_accession_found']) > 5:
                print(f"  ... and {len(fix_report['no_accession_found']) - 5} more")
    
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
    
    elif args.analyze_missing_genome_info:
        # Analyze missing genome information
        print("Analyzing missing genome information...")
        analysis = analyze_missing_genome_info()
        
        print("\n" + "=" * 80)
        print("MISSING GENOME INFORMATION ANALYSIS")
        print("=" * 80)
        
        summary = analysis['summary']
        print(f"\n📊 SUMMARY:")
        print(f"   Total genomes: {summary['total_genomes']}")
        print(f"   Missing taxonomy_id: {summary['missing_taxonomy_id_count']}")
        print(f"   Missing ome: {summary['missing_ome_count']}")
        print(f"   Missing both: {summary['missing_both_count']}")
        print(f"   Potential fixes: {summary['potential_fixes_count']}")
        
        if analysis['potential_fixes']:
            print(f"\n🔧 POTENTIAL FIXES:")
            for fix in analysis['potential_fixes'][:10]:  # Show first 10
                print(f"   {fix['description']}")
            if len(analysis['potential_fixes']) > 10:
                print(f"   ... and {len(analysis['potential_fixes']) - 10} more")
        
        if analysis['missing_taxonomy_id']:
            print(f"\n❌ GENOMES MISSING TAXONOMY_ID:")
            for genome in analysis['missing_taxonomy_id'][:5]:  # Show first 5
                print(f"   Genome {genome['genome_id']}: ome='{genome['ome']}'")
            if len(analysis['missing_taxonomy_id']) > 5:
                print(f"   ... and {len(analysis['missing_taxonomy_id']) - 5} more")
        
        if analysis['missing_ome']:
            print(f"\n❌ GENOMES MISSING OME:")
            for genome in analysis['missing_ome'][:5]:  # Show first 5
                print(f"   Genome {genome['genome_id']}: taxonomy_id={genome['taxonomy_id']}")
            if len(analysis['missing_ome']) > 5:
                print(f"   ... and {len(analysis['missing_ome']) - 5} more")
        
        print(f"\n💡 RECOMMENDATION: {summary['recommendation']}")
    
    elif args.fix_missing_genome_taxonomy_from_joined:
        print("Fixing missing genome taxonomy_id from joined_ships (mode tax_id)...")
        fix_report = fix_missing_genome_taxonomy_from_joined_ships(dry_run=not args.apply)
        
        print("\n" + "=" * 80)
        print("GENOME TAXONOMY FROM JOINED_SHIPS FIX REPORT")
        print("=" * 80)
        
        summary = fix_report['summary']
        print(f"\n📈 FIXES:")
        print(f"   Total missing: {summary['total_missing']}")
        print(f"   Total updated: {summary['total_updated']}")
        print(f"   Total unresolved: {summary['total_unresolved']}")
        
        if fix_report['genomes_updated']:
            print("\n✅ UPDATED (first 10):")
            for item in fix_report['genomes_updated'][:10]:
                print(f"   Genome {item['genome_id']} ome='{item['ome']}' -> taxonomy_id={item['new_taxonomy_id']}")
        
        if fix_report['no_tax_id_found']:
            print("\n⚠️  UNRESOLVED (first 10):")
            for item in fix_report['no_tax_id_found'][:10]:
                print(f"   Genome {item['genome_id']} ome='{item['ome']}': {item['issue']}")

    elif args.fix_missing_genome_taxonomy_via_ncbi:
        if not args.ome_map:
            print("Error: --ome-map is required for --fix-missing-genome-taxonomy-via-ncbi")
        else:
            print("Fixing missing genome taxonomy_id via NCBI...")
            fix_report = fix_missing_genome_taxonomy_via_ncbi(
                map_csv=args.ome_map,
                email=args.ncbi_email,
                api_key=args.ncbi_api_key,
                dry_run=not args.apply
            )

            print("\n" + "=" * 80)
            print("GENOME TAXONOMY VIA NCBI FIX REPORT")
            print("=" * 80)
            
            summary = fix_report['summary']
            print(f"\n📈 FIXES:")
            print(f"   Total missing: {summary['total_missing']}")
            print(f"   Linked: {summary['linked']}")
            print(f"   Taxa created: {summary['taxa_created']}")
            print(f"   Skipped (no name): {summary['skipped_no_name']}")
            print(f"   Skipped (not found): {summary['skipped_not_found']}")

            if fix_report['genomes_linked']:
                print("\n✅ LINKED (first 10):")
                for item in fix_report['genomes_linked'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}' -> taxonomy_id={item['taxonomy_id']} ({item['name']})")

            if fix_report['taxa_created']:
                print("\n🆕 TAXA CREATED (first 10):")
                for item in fix_report['taxa_created'][:10]:
                    print(f"   Taxonomy {item['taxonomy_id']}: {item['name']} (taxID={item['taxID']})")

            if fix_report['skipped_no_name']:
                print("\n⚠️  SKIPPED (no name) (first 10):")
                for item in fix_report['skipped_no_name'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}'")

            if fix_report['skipped_not_found']:
                print("\n⚠️  SKIPPED (not found) (first 10):")
                for item in fix_report['skipped_not_found'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}', query='{item['query']}'")

    elif args.fix_missing_genome_taxonomy_via_taxdump:
        if not args.names_dmp:
            print("Error: --names-dmp is required for --fix-missing-genome-taxonomy-via-taxdump")
        else:
            print("Fixing missing genome taxonomy_id via taxdump names.dmp...")
            fix_report = fix_missing_genome_taxonomy_via_taxdump(
                names_dmp_path=args.names_dmp,
                include_synonyms=args.include_synonyms,
                dry_run=not args.apply
            )

            print("\n" + "=" * 80)
            print("GENOME TAXONOMY VIA TAXDUMP FIX REPORT")
            print("=" * 80)
            
            summary = fix_report['summary']
            print(f"\n📈 FIXES:")
            print(f"   Total missing: {summary['total_missing']}")
            print(f"   Linked: {summary['linked']}")
            print(f"   Taxa created: {summary['taxa_created']}")
            print(f"   Skipped (no match): {summary['skipped_no_match']}")

            if fix_report['genomes_linked']:
                print("\n✅ LINKED (first 10):")
                for item in fix_report['genomes_linked'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}' -> taxonomy_id={item['taxonomy_id']} ({item['name']})")

            if fix_report['taxa_created']:
                print("\n🆕 TAXA CREATED (first 10):")
                for item in fix_report['taxa_created'][:10]:
                    print(f"   Taxonomy {item['taxonomy_id']}: {item['name']} (taxID={item['taxID']})")

    elif args.fix_missing_genome_info_via_ome_map:
        if not args.ome_map:
            print("Error: --ome-map is required for --fix-missing-genome-info-via-ome-map")
        else:
            print("Reconciling genome info according to OME mapping TSV...")
            fix_report = reconcile_genome_taxonomy_via_map(
                map_path=args.ome_map,
                dry_run=not args.apply
            )

            print("\n" + "=" * 80)
            print("GENOME INFO RECONCILIATION REPORT")
            print("=" * 80)
            
            summary = fix_report['summary']
            print(f"\n📈 CHANGES:")
            print(f"   Checked: {summary['checked']}")
            print(f"   Updated: {summary['updated']}")
            print(f"   Taxa created: {summary['taxa_created']}")
            print(f"   Assembly accessions updated: {summary['assembly_accessions_updated']}")
            print(f"   Skipped (no map): {summary['skipped_no_map']}")

            if fix_report['genomes_updated']:
                print("\n✅ UPDATED (first 10):")
                for item in fix_report['genomes_updated'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}' -> taxonomy_id {item['old_taxonomy_id']} => {item['new_taxonomy_id']} ({item['name']})")

            if fix_report['assembly_accessions_updated']:
                print("\n🔄 ASSEMBLY ACCESSIONS UPDATED (first 10):")
                for item in fix_report['assembly_accessions_updated'][:10]:
                    print(f"   Genome {item['genome_id']} ome='{item['ome']}' -> assembly_accession '{item['old_assembly_accession']}' => '{item['new_assembly_accession']}'")

    elif args.lookup_assemblies_from_contigs:
        print("Looking up assembly accessions from contigIDs via NCBI...")
        result = lookup_assembly_accessions_from_contigs(
            email=args.ncbi_email,
            api_key=args.ncbi_api_key,
            limit=args.lookup_limit,
            dry_run=not args.apply
        )

        print("\n" + "=" * 80)
        print("ASSEMBLY LOOKUP FROM CONTIGIDS")
        print("=" * 80)
        print(f"\nQueried: {result['queried']} of {result['total_contigs']} unique contigs")
        print(f"Resolved assemblies: {result['resolved']}")
        print(f"Skipped empty: {result['skipped_empty']}")
        if result.get('examples'):
            print("\nExamples:")
            for ex in result['examples']:
                print(f"  contigID='{ex['contigID']}', term='{ex['normalized_term']}', assembly='{ex['assemblyaccession']}', organism='{ex['organism']}'")

    elif args.enable_foreign_keys:
        print("Enabling foreign key constraints...")
        result = enable_foreign_keys(dry_run=not args.apply)

        print("\n" + "=" * 80)
        print("FOREIGN KEY CONSTRAINTS ENABLED")
        print("=" * 80)
        print(f"\nStatus: {result['status']}")
        for detail in result['details']:
            print(f"  {detail}")

    elif args.disable_foreign_keys:
        print("Disabling foreign key constraints...")
        result = disable_foreign_keys(dry_run=not args.apply)

        print("\n" + "=" * 80)
        print("FOREIGN KEY CONSTRAINTS DISABLED")
        print("=" * 80)
        print(f"\nStatus: {result['status']}")
        for detail in result['details']:
            print(f"  {detail}")

    elif args.check_foreign_key_enforcement:
        print("Checking foreign key enforcement...")
        result = check_foreign_key_enforcement()

        print("\n" + "=" * 80)
        print("FOREIGN KEY ENFORCEMENT CHECK")
        print("=" * 80)
        print(f"\nForeign keys enabled: {result['foreign_keys_enabled']}")
        print(f"Tables with constraints: {len(result['tables_with_constraints'])}")
        print(f"Tables missing constraints: {len(result['missing_constraints'])}")

        if result['tables_with_constraints']:
            print("\nTables with foreign key constraints:")
            for table_info in result['tables_with_constraints']:
                print(f"  {table_info['table']}: {table_info['constraints']} constraints")

        if result['missing_constraints']:
            print("\nTables missing expected foreign key constraints:")
            for missing in result['missing_constraints']:
                print(f"  {missing['table']}: {len(missing['expected_constraints'])} expected constraints")

        if result['recommendations']:
            print("\nRecommendations:")
            for rec in result['recommendations']:
                print(f"  - {rec}")

    elif args.fix_joined_ships_ship_id:
        print("Fixing joined_ships.ship_id foreign key...")
        result = fix_joined_ships_ship_id_foreign_key(dry_run=not args.apply)

        print("\n" + "=" * 80)
        print("JOINED_SHIPS SHIP_ID FOREIGN KEY FIX")
        print("=" * 80)
        print(f"\n⚠️  CAUTION: This function maps accession IDs to ship IDs")
        print(f"Ship IDs fixed: {result['summary']['ship_ids_fixed']}")
        print(f"Invalid references found: {result['summary']['invalid_references']}")
        print(f"Accessions without ships: {result['summary']['no_ship_for_accession']}")

        if result['ship_ids_fixed']:
            print(f"\n✅ SHIP IDs FIXED (first 5):")
            for fix in result['ship_ids_fixed'][:5]:
                print(f"  {fix['starshipID']}: {fix['old_ship_id']} -> {fix['new_ship_id']} (accession: {fix['accession_tag']})")

        if result['no_ship_for_accession']:
            print(f"\n⚠️  ACCESSIONS WITHOUT SHIPS (first 5):")
            for item in result['no_ship_for_accession'][:5]:
                print(f"  {item['starshipID']}: accession {item['accession_tag']} has no ship record")

        if result['invalid_references_found']:
            print(f"\n❌ INVALID REFERENCES FOUND (first 5):")
            for invalid in result['invalid_references_found'][:5]:
                print(f"  {invalid['starshipID']}: invalid ship_id {invalid['invalid_ship_id']}")

        print(f"\n💡 Recommendation: {result['summary']['recommendation']}")
        print(f"⚠️  Warning: {result['summary']['warning']}")

    elif args.validate_all_foreign_keys:
        print("Performing comprehensive foreign key validation...")
        result = validate_all_foreign_keys()

        print("\n" + "=" * 80)
        print("COMPREHENSIVE FOREIGN KEY VALIDATION")
        print("=" * 80)
        print(f"\nForeign keys enabled: {result['foreign_keys_enabled']}")
        print(f"Relationships checked: {result['total_relationships_checked']}")
        print(f"Total violations found: {result['violations_found']}")
        print(f"Relationships with violations: {result['summary']['relationships_with_violations']}")

        if result['relationships']:
            print("\nRelationship details:")
            for rel_name, rel_info in result['relationships'].items():
                if isinstance(rel_info, dict) and 'orphaned_records' in rel_info:
                    if rel_info['orphaned_records'] > 0:
                        print(f"  ❌ {rel_name}: {rel_info['orphaned_records']} orphaned records")
                        if rel_info['examples']:
                            print(f"     Examples: {rel_info['examples'][:3]}")
                    else:
                        print(f"  ✅ {rel_name}: OK")
                elif 'error' in rel_info:
                    print(f"  ⚠️  {rel_name}: Error - {rel_info['error']}")

        if result['recommendations']:
            print("\nRecommendations:")
            for rec in result['recommendations']:
                print(f"  - {rec}")

    elif args.show_issues_summary:
        print("Showing cleanup_issues table summary...")
        result = get_cleanup_issues_summary()

        print("\n" + "=" * 80)
        print("CLEANUP ISSUES SUMMARY")
        print("=" * 80)

        print(f"\n📊 OVERVIEW:")
        print(f"   Total issues: {result['total_issues']}")
        print(f"   Open issues: {result['open_issues']}")
        print(f"   Fixed issues: {result['fixed_issues']}")
        print(f"   Issues by category: {result['issues_by_category']}")

        if result['recent_issues']:
            print(f"\n🕒 RECENT ISSUES (last 10):")
            for issue in result['recent_issues']:
                status_icon = "✅" if issue['status'] == 'FIXED' else "🔓"
                print(f"   {status_icon} {issue['category']}/{issue['issue_type']}: {issue['table_name']}.{issue['record_id']} ({issue['created_at']})")

        print(f"\n📈 STATUS BREAKDOWN:")
        for status, count in result['status_breakdown'].items():
            print(f"   {status}: {count}")

    elif args.identify_taxonomy_duplicates:
        # Identify taxonomy duplicates
        print("Identifying taxonomy duplicates...")
        duplicates = identify_taxonomy_duplicates()
        
        print("\n" + "=" * 80)
        print("TAXONOMY DUPLICATES ANALYSIS")
        print("=" * 80)
        
        print(f"\nTaxID duplicate groups: {duplicates['summary']['taxid_duplicate_groups']}")
        print(f"Total taxID duplicates: {duplicates['summary']['total_taxid_duplicates']}")
        print(f"Consistent duplicate groups: {duplicates['summary']['consistent_duplicate_groups']}")
        print(f"Total consistent duplicates: {duplicates['summary']['total_consistent_duplicates']}")
        print(f"Total entries to remove: {duplicates['summary']['total_entries_to_remove']}")
        
        if duplicates['taxid_duplicates']:
            print("\nTaxID DUPLICATE GROUPS:")
            for group in duplicates['taxid_duplicates'][:10]:  # Show first 10
                status = "✅ CONSISTENT" if group['consistent'] else "❌ INCONSISTENT"
                print(f"  TaxID {group['taxID']} ({status}): {group['count']} entries")
                print(f"    Primary: ID {group['primary_id']}")
                print(f"    Duplicates: {group['duplicate_ids']}")
                if not group['consistent']:
                    print(f"    Conflicts: {len(group['inconsistent_fields'])} field(s)")
                for entry in group['entries']:
                    print(f"      - ID {entry['id']}: {entry['name']} ({entry['genus']} {entry['species']} sect. {entry['section'] or 'N/A'})")
                print()
            if len(duplicates['taxid_duplicates']) > 10:
                print(f"  ... and {len(duplicates['taxid_duplicates']) - 10} more groups")
        
        if duplicates['consistent_duplicates']:
            print(f"\nCONSISTENT SEMANTIC DUPLICATES:")
            for group in duplicates['consistent_duplicates'][:5]:  # Show first 5
                print(f"  Group: {group['count']} entries with different taxIDs")
                print(f"    Primary: ID {group['primary_id']}")
                print(f"    Duplicates: {group['duplicate_ids']}")
                print(f"    TaxIDs: {group['taxids']}")
                for entry in group['entries']:
                    print(f"      - ID {entry['id']}: {entry['name']} (taxID: {entry['taxID']})")
                print()
            if len(duplicates['consistent_duplicates']) > 5:
                print(f"  ... and {len(duplicates['consistent_duplicates']) - 5} more groups")

    elif args.fix_missing_tax_id_via_ome:
        # Fix missing tax_id using ome code consistency
        print("Fixing missing tax_id in joined_ships using ome code consistency...")
        fix_report = fix_missing_tax_id_via_ome_consistency(
            dry_run=not args.apply,
            ome_map_path=args.ome_map
        )
        
        print("\n" + "=" * 80)
        print("OME TAX_ID CONSISTENCY FIX")
        print("=" * 80)
        
        print(f"\nOME groups found: {fix_report['summary']['ome_groups_found']}")
        print(f"OME groups analyzed: {fix_report['summary']['ome_groups_analyzed']}")
        print(f"Tax_ids filled: {fix_report['summary']['tax_ids_filled']}")
        print(f"Inconsistent groups: {fix_report['summary']['inconsistent_groups']}")
        print(f"Genomes created: {fix_report['summary']['genomes_created']}")
        print(f"OME map used: {fix_report['summary']['ome_map_used']}")
        print(f"Warnings: {fix_report['summary']['warnings']}")
        
        if fix_report['tax_ids_filled']:
            print(f"\nTAX_IDS FILLED (first 10):")
            for fix in fix_report['tax_ids_filled'][:10]:
                print(f"  {fix['starshipID']}: ome='{fix['ome_code']}' → tax_id={fix['assigned_tax_id']}")
            if len(fix_report['tax_ids_filled']) > 10:
                print(f"  ... and {len(fix_report['tax_ids_filled']) - 10} more")
        
        if fix_report['ome_groups_analyzed']:
            print(f"\nOME GROUPS ANALYZED (first 10):")
            for group in fix_report['ome_groups_analyzed'][:10]:
                print(f"  {group['ome_code']}: {group['total_entries']} entries, {group['status']}")
                if group['status'] == 'consistent':
                    print(f"    → Filled {group['filled']} missing tax_ids with {group['consensus_tax_id']}")
                elif group['status'] == 'created_from_ome_map':
                    print(f"    → Created genome {group['genome_id']}, filled {group['filled']} missing tax_ids with {group['consensus_tax_id']}")
            if len(fix_report['ome_groups_analyzed']) > 10:
                print(f"  ... and {len(fix_report['ome_groups_analyzed']) - 10} more")
        
        if fix_report['genomes_created']:
            print(f"\nGENOMES CREATED FROM OME MAP (first 10):")
            for genome in fix_report['genomes_created'][:10]:
                print(f"  {genome['ome_code']}: genome_id={genome['genome_id']}, taxonomy_id={genome['taxonomy_id']} ({genome['name']})")
            if len(fix_report['genomes_created']) > 10:
                print(f"  ... and {len(fix_report['genomes_created']) - 10} more")
        
        if fix_report['inconsistent_ome_groups']:
            print(f"\nINCONSISTENT OME GROUPS (manual review needed):")
            for group in fix_report['inconsistent_ome_groups'][:5]:
                print(f"  {group['ome_code']}: {group['total_entries']} entries, tax_ids={group['tax_ids']}")
            if len(fix_report['inconsistent_ome_groups']) > 5:
                print(f"  ... and {len(fix_report['inconsistent_ome_groups']) - 5} more")

    elif args.consolidate_taxonomy:
        # Run complete taxonomy consolidation pipeline
        print("Running comprehensive taxonomy consolidation...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.run_full_consolidation(
                clean_data=True,
                resolve_ids=True,
                backfill_hierarchy=True,
                validate_consistency=True,
                dry_run=not args.apply,
                names_dmp_path=args.names_dmp,
                ncbi_email=args.ncbi_email,
                ncbi_api_key=args.ncbi_api_key,
                include_synonyms=args.include_synonyms
            )
        
        print("\n" + "=" * 80)
        print("TAXONOMY CONSOLIDATION SUMMARY")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nTotal operations: {summary['total_operations']}")
        print(f"Successful operations: {summary['successful_operations']}")
        print(f"Issues found: {summary['issues_found']}")
        
        if summary['recommendations']:
            print("\nRECOMMENDATIONS:")
            for rec in summary['recommendations']:
                print(f"  • {rec}")
        
        # Show detailed results for each phase
        if report['cleaning']:
            cleaning = report['cleaning']['summary']
            print(f"\nCLEANING PHASE:")
            print(f"  Whitespace fixes: {cleaning['whitespace_fixes']}")
            print(f"  Species cleanups: {cleaning['species_cleanup']}")
            print(f"  Strain cleanups: {cleaning['strain_cleanup']}")
            print(f"  Name cleanups: {cleaning['name_cleanup']}")
            print(f"  Hierarchy issues: {cleaning['hierarchy_issues']}")
        
        if report['id_resolution']:
            id_res = report['id_resolution']['summary']
            print(f"\nID RESOLUTION PHASE:")
            print(f"  Total missing taxIDs: {id_res['total_missing_taxid']}")
            print(f"  Resolved: {id_res['resolved']}")
            print(f"  Unresolved: {id_res['unresolved']}")
            print(f"  From existing: {id_res['resolved_from_existing']}")
            print(f"  From taxdump: {id_res['resolved_from_taxdump']}")
            print(f"  From NCBI: {id_res['resolved_from_ncbi']}")
        
        if report['hierarchy_backfill']:
            backfill = report['hierarchy_backfill']['summary']
            print(f"\nHIERARCHY BACKFILL PHASE:")
            print(f"  Fields backfilled: {backfill['fields_backfilled']}")
            print(f"  Consistency issues: {backfill['consistency_issues']}")
        
        if report['consistency_validation']:
            validation = report['consistency_validation']['summary']
            print(f"\nCONSISTENCY VALIDATION PHASE:")
            print(f"  Family inconsistencies: {validation['family_inconsistencies']}")
            print(f"  Hierarchy violations: {validation['hierarchy_violations']}")
            print(f"  Orphaned entries: {validation['orphaned_entries']}")
            print(f"  Duplicate entries: {validation['duplicate_entries']}")
            print(f"  Data quality issues: {validation['data_quality_issues']}")

    elif args.clean_taxonomy_data:
        # Clean taxonomy data only
        print("Cleaning taxonomy data...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.clean_taxonomy_data(dry_run=not args.apply)
        
        print("\n" + "=" * 80)
        print("TAXONOMY DATA CLEANING")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nTotal processed: {summary['total_processed']}")
        print(f"Whitespace fixes: {summary['whitespace_fixes']}")
        print(f"Species cleanups: {summary['species_cleanup']}")
        print(f"Strain cleanups: {summary['strain_cleanup']}")
        print(f"Name cleanups: {summary['name_cleanup']}")
        print(f"Hierarchy issues: {summary['hierarchy_issues']}")
        
        if report['whitespace_fixes']:
            print(f"\nWHITESPACE FIXES (first 10):")
            for fix in report['whitespace_fixes'][:10]:
                print(f"  ID {fix['taxonomy_id']} {fix['field']}: '{fix['old_value']}' → '{fix['new_value']}'")
            if len(report['whitespace_fixes']) > 10:
                print(f"  ... and {len(report['whitespace_fixes']) - 10} more")
        
        if report['species_cleanup']:
            print(f"\nSPECIES CLEANUP (first 10):")
            for cleanup in report['species_cleanup'][:10]:
                print(f"  ID {cleanup['taxonomy_id']}: '{cleanup['old_species']}' → '{cleanup['new_species']}' (genus: {cleanup['genus']})")
            if len(report['species_cleanup']) > 10:
                print(f"  ... and {len(report['species_cleanup']) - 10} more")

    elif args.resolve_taxonomy_ids:
        # Resolve taxonomy IDs only
        print("Resolving taxonomy IDs...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.resolve_taxonomy_ids(
                dry_run=not args.apply,
                names_dmp_path=args.names_dmp,
                ncbi_email=args.ncbi_email,
                ncbi_api_key=args.ncbi_api_key,
                include_synonyms=args.include_synonyms
            )
        
        print("\n" + "=" * 80)
        print("TAXONOMY ID RESOLUTION")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nTotal missing taxIDs: {summary['total_missing_taxid']}")
        print(f"Resolved: {summary['resolved']}")
        print(f"Unresolved: {summary['unresolved']}")
        print(f"From existing: {summary['resolved_from_existing']}")
        print(f"From taxdump: {summary['resolved_from_taxdump']}")
        print(f"From NCBI: {summary['resolved_from_ncbi']}")
        
        if report['resolved_from_existing']:
            print(f"\nRESOLVED FROM EXISTING (first 10):")
            for resolved in report['resolved_from_existing'][:10]:
                print(f"  ID {resolved['taxonomy_id']}: {resolved['name']} → taxID {resolved['taxID']}")
            if len(report['resolved_from_existing']) > 10:
                print(f"  ... and {len(report['resolved_from_existing']) - 10} more")
        
        if report['unresolved']:
            print(f"\nUNRESOLVED (first 10):")
            for unresolved in report['unresolved'][:10]:
                print(f"  ID {unresolved['taxonomy_id']}: {unresolved['name']} ({unresolved['reason']})")
            if len(report['unresolved']) > 10:
                print(f"  ... and {len(report['unresolved']) - 10} more")

    elif args.backfill_taxonomy_hierarchy:
        # Backfill taxonomy hierarchy only
        print("Backfilling taxonomy hierarchy...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.backfill_taxonomy_hierarchy(dry_run=not args.apply)
        
        print("\n" + "=" * 80)
        print("TAXONOMY HIERARCHY BACKFILL")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nTotal processed: {summary['total_processed']}")
        print(f"Fields backfilled: {summary['fields_backfilled']}")
        print(f"Consistency issues: {summary['consistency_issues']}")
        
        if report['fields_backfilled']:
            print(f"\nFIELDS BACKFILLED (first 10):")
            for backfill in report['fields_backfilled'][:10]:
                print(f"  ID {backfill['taxonomy_id']} {backfill['field']}: '{backfill['inferred_value']}'")
            if len(report['fields_backfilled']) > 10:
                print(f"  ... and {len(report['fields_backfilled']) - 10} more")

    elif args.validate_taxonomy_consistency:
        # Validate taxonomy consistency only
        print("Validating taxonomy consistency...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.validate_taxonomy_consistency()
        
        print("\n" + "=" * 80)
        print("TAXONOMY CONSISTENCY VALIDATION")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nFamily inconsistencies: {summary['family_inconsistencies']}")
        print(f"Hierarchy violations: {summary['hierarchy_violations']}")
        print(f"Orphaned entries: {summary['orphaned_entries']}")
        print(f"Duplicate entries: {summary['duplicate_entries']}")
        print(f"Data quality issues: {summary['data_quality_issues']}")
        
        if report['family_inconsistencies']:
            print(f"\nFAMILY INCONSISTENCIES (first 10):")
            for issue in report['family_inconsistencies'][:10]:
                print(f"  Family '{issue['family']}' {issue['field']}: {issue['values']} ({issue['count']} entries)")
            if len(report['family_inconsistencies']) > 10:
                print(f"  ... and {len(report['family_inconsistencies']) - 10} more")
        
        if report['duplicate_entries']:
            print(f"\nDUPLICATE ENTRIES (first 10):")
            for dup in report['duplicate_entries'][:10]:
                print(f"  {dup['genus']} {dup['species']} {dup['strain']} (name: {dup['name']}): {dup['count']} entries")
                for entry in dup['entries'][:3]:
                    print(f"    ID {entry['id']}: {entry['name']} (taxID: {entry['taxID']})")
                if len(dup['entries']) > 3:
                    print(f"    ... and {len(dup['entries']) - 3} more")
            if len(report['duplicate_entries']) > 10:
                print(f"  ... and {len(report['duplicate_entries']) - 10} more")

    elif args.consolidate_taxonomy_duplicates:
        # Consolidate duplicate taxonomy entries
        print("Consolidating duplicate taxonomy entries...")
        with TaxonomyConsolidator() as consolidator:
            report = consolidator.consolidate_duplicate_entries(dry_run=not args.apply)
        
        print("\n" + "=" * 80)
        print("TAXONOMY DUPLICATE CONSOLIDATION")
        print("=" * 80)
        
        summary = report['summary']
        print(f"\nDuplicates found: {summary['duplicates_found']}")
        print(f"Consolidations performed: {summary['consolidations_performed']}")
        print(f"Foreign key updates: {summary['foreign_key_updates']}")
        
        if report['consolidations_performed']:
            print(f"\nCONSOLIDATIONS PERFORMED (first 10):")
            for consolidation in report['consolidations_performed'][:10]:
                print(f"  {consolidation['genus']} {consolidation['species']} {consolidation['strain']}")
                print(f"    Keep: ID {consolidation['keep_entry']['id']} - {consolidation['keep_entry']['name']} (taxID: {consolidation['keep_entry']['taxID']})")
                for remove_entry in consolidation['remove_entries']:
                    print(f"    Remove: ID {remove_entry['id']} - {remove_entry['name']} (taxID: {remove_entry['taxID']})")
            if len(report['consolidations_performed']) > 10:
                print(f"  ... and {len(report['consolidations_performed']) - 10} more")
        
        if report['foreign_key_updates']:
            print(f"\nFOREIGN KEY UPDATES (first 10):")
            for fk_update in report['foreign_key_updates'][:10]:
                print(f"  {fk_update['table']}: {fk_update['records_updated']} records updated")
                print(f"    {fk_update['old_taxonomy_id']} → {fk_update['new_taxonomy_id']}")
            if len(report['foreign_key_updates']) > 10:
                print(f"  ... and {len(report['foreign_key_updates']) - 10} more")

    else:
        # Run full cleanup process
        main(dry_run=not args.apply, output_report=args.report, apply_fixes=args.apply_fixes, skip_accession_cleanup=args.skip_accession_cleanup, no_record_issues=args.no_record_issues) 
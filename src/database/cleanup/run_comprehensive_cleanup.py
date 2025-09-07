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

    from .utils.database_cleanup import *

except ImportError:
    # Fallback for when running the script directly
    import sys
    import os
    sys.path.append(_PROJECT_ROOT)
    
    from src.config.logging import get_logger
    from src.database.cleanup.utils.cleanup_accessions import main as run_accession_cleanup

    from src.database.cleanup.utils.database_cleanup import *

logger = get_logger(__name__)


def main(dry_run: bool = True, output_report: str = None, apply_fixes: bool = False, skip_accession_cleanup: bool = False):
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
    Step 9: Generate comprehensive report
    Step 10: Apply fixes (if requested)
    
    Args:
        dry_run (bool): If True, only analyze and report without making changes
        output_report (str): Path to save the cleanup report
        apply_fixes (bool): If True, apply fixes to the identified issues
        skip_accession_cleanup (bool): If True, skip the slow accession cleanup step
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
    
    logger.info("Step 9: Generating comprehensive report...")
    report = generate_cleanup_report(all_issues)
    
    if output_report:
        with open(output_report, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to {output_report}")
    else:
        print(report)

    # Record issues into tracking table (table creation happens even in dry-run)
    try:
        logger.info("Recording issues into cleanup_issues table...")
        create_cleanup_issues_table()
        record_summary = record_cleanup_issues(all_issues, source='pipeline', dry_run=dry_run)
        logger.info(f"Recorded issues summary: inserted={record_summary['inserted']} skipped={record_summary['skipped']}")
    except Exception as e:
        logger.error(f"Failed to record cleanup issues: {str(e)}")
    
    if apply_fixes:
        logger.info("Step 10: Applying fixes to identified issues...")
        
        # Apply general relationship fixes
        fixes_report = fix_ships_accessions_joined_ships_relationships(dry_run=False)
        
        # Apply ship_id relationship fixes
        ship_id_fixes = fix_ship_id_relationships(dry_run=False)
        
        logger.info("Fixes applied successfully")
        
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
                       help="Fill missing Genome.taxonomy_id by querying NCBI Taxonomy (requires --ome-name-map)")
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
        if not args.ome_name_map:
            print("Error: --ome-name-map is required for --fix-missing-genome-taxonomy-via-ncbi")
        else:
            print("Fixing missing genome taxonomy_id via NCBI...")
            fix_report = fix_missing_genome_taxonomy_via_ncbi(
                map_csv=args.ome_name_map,
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
            print("Error: --ome-map is required for --fix-genome-info-via-ome-map")
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

    else:
        # Run full cleanup process
        main(dry_run=not args.apply, output_report=args.report, apply_fixes=args.apply_fixes, skip_accession_cleanup=args.skip_accession_cleanup) 
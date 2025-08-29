#!/usr/bin/env python3
"""
Comprehensive database cleaning script for Starbase.

This script performs all database cleaning tasks outlined in the plan:
1. Accession tag correction
2. Genome table validation and updates
3. Taxonomic table consistency checks
4. Genomic features validation
5. Foreign key consistency updates
6. General schema validation
"""

import pandas as pd
from sqlalchemy import text
from typing import Dict
from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.database.models.schema import (
    Accessions, Ships, JoinedShips, Genome, Taxonomy, 
    StarshipFeatures, Captains
)
from src.database.cleanup_accessions import main as run_accession_cleanup
import time

logger = get_logger(__name__)


def check_genome_table() -> Dict:
    """
    Check genome table for missing information and inconsistencies.
    
    Returns:
        Dict: Report of genome table issues found
    """
    logger.info("Checking genome table...")
    session = StarbaseSession()
    
    issues = {
        'missing_genome_info': [],
        'inconsistent_taxonomy': [],
        'orphaned_genomes': [],
        'duplicate_genomes': []
    }
    
    try:
        # Check for genomes with missing essential information
        missing_info = session.query(Genome).filter(
            (Genome.ome.is_(None)) | 
            (Genome.ome == '') |
            (Genome.taxonomy_id.is_(None))
        ).all()
        
        for genome in missing_info:
            issues['missing_genome_info'].append({
                'id': genome.id,
                'ome': genome.ome,
                'taxonomy_id': genome.taxonomy_id,
                'issue': 'Missing ome or taxonomy_id'
            })
        
        # Check for orphaned genomes (no associated accessions)
        # * Note: This relationship might not be direct, so we'll check differently
        orphaned = session.query(Genome).outerjoin(JoinedShips, Genome.id == JoinedShips.genome_id).filter(
            JoinedShips.id.is_(None)
        ).all()
        
        for genome in orphaned:
            issues['orphaned_genomes'].append({
                'id': genome.id,
                'ome': genome.ome,
                'issue': 'No associated accessions'
            })
        
        # Check for duplicate genome entries
        duplicates = session.query(Genome.ome, Genome.taxonomy_id).group_by(
            Genome.ome, Genome.taxonomy_id
        ).having(text('COUNT(*) > 1')).all()
        
        for ome, tax_id in duplicates:
            issues['duplicate_genomes'].append({
                'ome': ome,
                'taxonomy_id': tax_id,
                'issue': 'Duplicate genome entry'
            })
        
        logger.info(f"Found {len(issues['missing_genome_info'])} genomes with missing info")
        logger.info(f"Found {len(issues['orphaned_genomes'])} orphaned genomes")
        logger.info(f"Found {len(issues['duplicate_genomes'])} duplicate genome entries")
        
    except Exception as e:
        logger.error(f"Error checking genome table: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def check_taxonomic_table() -> Dict:
    """
    Check taxonomic table for internal consistency.
    
    Returns:
        Dict: Report of taxonomic table issues found
    """
    logger.info("Checking taxonomic table...")
    session = StarbaseSession()
    
    issues = {
        'missing_taxonomic_levels': [],
        'inconsistent_hierarchy': [],
        'orphaned_taxonomy': [],
        'duplicate_taxonomy': []
    }
    
    try:
        # Check for taxonomy entries with missing essential information
        missing_levels = session.query(Taxonomy).filter(
            (Taxonomy.name.is_(None)) | 
            (Taxonomy.name == '') |
            (Taxonomy.superkingdom.is_(None)) |
            (Taxonomy.kingdom.is_(None))
        ).all()
        
        for tax in missing_levels:
            issues['missing_taxonomic_levels'].append({
                'id': tax.id,
                'name': tax.name,
                'superkingdom': tax.superkingdom,
                'kingdom': tax.kingdom,
                'issue': 'Missing essential taxonomic information'
            })
        
        # Check for orphaned taxonomy (no associated genomes)
        orphaned = session.query(Taxonomy).outerjoin(Genome, Taxonomy.id == Genome.taxonomy_id).filter(
            Genome.id.is_(None)
        ).all()
        
        for tax in orphaned:
            issues['orphaned_taxonomy'].append({
                'id': tax.id,
                'name': tax.name,
                'issue': 'No associated genomes'
            })
        
        # Check for duplicate taxonomy entries
        duplicates = session.query(Taxonomy.name, Taxonomy.superkingdom, Taxonomy.kingdom).group_by(
            Taxonomy.name, Taxonomy.superkingdom, Taxonomy.kingdom
        ).having(text('COUNT(*) > 1')).all()
        
        for name, superkingdom, kingdom in duplicates:
            issues['duplicate_taxonomy'].append({
                'name': name,
                'superkingdom': superkingdom,
                'kingdom': kingdom,
                'issue': 'Duplicate taxonomy entry'
            })
        
        logger.info(f"Found {len(issues['missing_taxonomic_levels'])} taxonomy entries with missing levels")
        logger.info(f"Found {len(issues['orphaned_taxonomy'])} orphaned taxonomy entries")
        logger.info(f"Found {len(issues['duplicate_taxonomy'])} duplicate taxonomy entries")
        
    except Exception as e:
        logger.error(f"Error checking taxonomic table: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def check_genomic_features() -> Dict:
    """
    Check genomic features for coordinate consistency and genome linkage.
    
    Returns:
        Dict: Report of genomic features issues found
    """
    logger.info("Checking genomic features...")
    session = StarbaseSession()
    
    issues = {
        'invalid_coordinates': [],
        'orphaned_features': [],
        'inconsistent_strands': [],
        'overlapping_features': []
    }
    
    try:
        # Check for features with invalid coordinates
        invalid_coords = session.query(StarshipFeatures).filter(
            (StarshipFeatures.elementBegin < 0) |
            (StarshipFeatures.elementEnd < 0) |
            (StarshipFeatures.elementBegin >= StarshipFeatures.elementEnd)
        ).all()
        
        for feature in invalid_coords:
            issues['invalid_coordinates'].append({
                'id': feature.id,
                'contigID': feature.contigID,
                'elementBegin': feature.elementBegin,
                'elementEnd': feature.elementEnd,
                'issue': 'Invalid coordinate range'
            })
        
        # Check for orphaned features (no associated ships)
        orphaned = session.query(StarshipFeatures).outerjoin(Accessions, StarshipFeatures.ship_id == Accessions.id).filter(
            Accessions.id.is_(None)
        ).all()
        
        for feature in orphaned:
            issues['orphaned_features'].append({
                'id': feature.id,
                'contigID': feature.contigID,
                'issue': 'No associated ship'
            })
        
        # Check for features with inconsistent strand information
        inconsistent_strands = session.query(StarshipFeatures).filter(
            (StarshipFeatures.strand.notin_(['+', '-', None])) &
            (StarshipFeatures.strand != '')
        ).all()
        
        for feature in inconsistent_strands:
            issues['inconsistent_strands'].append({
                'id': feature.id,
                'contigID': feature.contigID,
                'strand': feature.strand,
                'issue': 'Invalid strand value'
            })
        
        logger.info(f"Found {len(issues['invalid_coordinates'])} features with invalid coordinates")
        logger.info(f"Found {len(issues['orphaned_features'])} orphaned features")
        logger.info(f"Found {len(issues['inconsistent_strands'])} features with inconsistent strands")
        
    except Exception as e:
        logger.error(f"Error checking genomic features: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def check_foreign_keys() -> Dict:
    """
    Check foreign key consistency across all tables.
    
    Returns:
        Dict: Report of foreign key issues found
    """
    logger.info("Checking foreign key consistency...")
    session = StarbaseSession()
    
    issues = {
        'broken_ship_accession_links': [],
        'broken_joined_ship_links': [],
        'broken_captain_links': [],
        'broken_genome_taxonomy_links': []
    }
    
    try:
        # Check for ships with broken accession links
        broken_ship_links = session.query(Ships).outerjoin(Accessions, Ships.accession_id == Accessions.id).filter(
            Accessions.id.is_(None)
        ).all()
        
        for ship in broken_ship_links:
            issues['broken_ship_accession_links'].append({
                'ship_id': ship.id,
                'accession_id': ship.accession_id,
                'issue': 'Ship references non-existent accession'
            })
        
        # Check for joined_ships with broken links
        broken_joined_links = session.query(JoinedShips).outerjoin(Accessions, JoinedShips.ship_id == Accessions.id).filter(
            Accessions.id.is_(None)
        ).all()
        
        for joined in broken_joined_links:
            issues['broken_joined_ship_links'].append({
                'joined_id': joined.id,
                'ship_id': joined.ship_id,
                'issue': 'Joined ship references non-existent accession'
            })
        
        # Check for captains with broken ship links
        broken_captain_links = session.query(Captains).outerjoin(Accessions, Captains.ship_id == Accessions.id).filter(
            Accessions.id.is_(None)
        ).all()
        
        for captain in broken_captain_links:
            issues['broken_captain_links'].append({
                'captain_id': captain.id,
                'ship_id': captain.ship_id,
                'issue': 'Captain references non-existent ship'
            })
        
        # Check for genomes with broken taxonomy links
        broken_genome_links = session.query(Genome).outerjoin(Taxonomy, Genome.taxonomy_id == Taxonomy.id).filter(
            Taxonomy.id.is_(None)
        ).all()
        
        for genome in broken_genome_links:
            issues['broken_genome_taxonomy_links'].append({
                'genome_id': genome.id,
                'taxonomy_id': genome.taxonomy_id,
                'issue': 'Genome references non-existent taxonomy'
            })
        
        logger.info(f"Found {len(issues['broken_ship_accession_links'])} broken ship-accession links")
        logger.info(f"Found {len(issues['broken_joined_ship_links'])} broken joined ship links")
        logger.info(f"Found {len(issues['broken_captain_links'])} broken captain links")
        logger.info(f"Found {len(issues['broken_genome_taxonomy_links'])} broken genome-taxonomy links")
        
    except Exception as e:
        logger.error(f"Error checking foreign keys: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def check_schema_violations() -> Dict:
    """
    Check for database entries that violate schema constraints.
    
    Returns:
        Dict: Report of schema violations found
    """
    logger.info("Checking schema violations...")
    session = StarbaseSession()
    
    issues = {
        'null_violations': [],
        'unique_violations': [],
        'type_violations': []
    }
    
    try:
        # Check for null violations in required fields
        null_violations = []
        
        # Check accessions table
        null_accessions = session.query(Accessions).filter(
            (Accessions.accession_tag.is_(None)) | 
            (Accessions.accession_tag == '')
        ).all()
        
        for acc in null_accessions:
            null_violations.append({
                'table': 'accessions',
                'id': acc.id,
                'field': 'accession_tag',
                'issue': 'Required field is null or empty'
            })
        
        # Check ships table
        null_ships = session.query(Ships).filter(
            (Ships.sequence.is_(None)) | 
            (Ships.sequence == '')
        ).all()
        
        for ship in null_ships:
            null_violations.append({
                'table': 'ships',
                'id': ship.id,
                'field': 'sequence',
                'issue': 'Required field is null or empty'
            })
        
        issues['null_violations'] = null_violations
        
        # Check for potential unique violations
        # This would require more complex analysis based on your specific unique constraints
        
        logger.info(f"Found {len(issues['null_violations'])} null/empty violations")
        
    except Exception as e:
        logger.error(f"Error checking schema violations: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def generate_cleanup_report(all_issues: Dict) -> str:
    """
    Generate a comprehensive report of all database issues found.
    
    Args:
        all_issues (Dict): Dictionary containing all issues found
        
    Returns:
        str: Formatted report
    """
    report = []
    report.append("=" * 80)
    report.append("COMPREHENSIVE DATABASE CLEANUP REPORT")
    report.append("=" * 80)
    report.append("")
    
    # Genome table issues
    genome_issues = all_issues.get('genome', {})
    report.append("GENOME TABLE ISSUES")
    report.append("-" * 50)
    report.append(f"Missing genome info: {len(genome_issues.get('missing_genome_info', []))}")
    report.append(f"Orphaned genomes: {len(genome_issues.get('orphaned_genomes', []))}")
    report.append(f"Duplicate genomes: {len(genome_issues.get('duplicate_genomes', []))}")
    report.append("")
    
    # Taxonomy table issues
    taxonomy_issues = all_issues.get('taxonomy', {})
    report.append("TAXONOMY TABLE ISSUES")
    report.append("-" * 50)
    report.append(f"Missing taxonomic levels: {len(taxonomy_issues.get('missing_taxonomic_levels', []))}")
    report.append(f"Orphaned taxonomy: {len(taxonomy_issues.get('orphaned_taxonomy', []))}")
    report.append(f"Duplicate taxonomy: {len(taxonomy_issues.get('duplicate_taxonomy', []))}")
    report.append("")
    
    # Genomic features issues
    features_issues = all_issues.get('features', {})
    report.append("GENOMIC FEATURES ISSUES")
    report.append("-" * 50)
    report.append(f"Invalid coordinates: {len(features_issues.get('invalid_coordinates', []))}")
    report.append(f"Orphaned features: {len(features_issues.get('orphaned_features', []))}")
    report.append(f"Inconsistent strands: {len(features_issues.get('inconsistent_strands', []))}")
    report.append("")
    
    # Foreign key issues
    fk_issues = all_issues.get('foreign_keys', {})
    report.append("FOREIGN KEY ISSUES")
    report.append("-" * 50)
    report.append(f"Broken ship-accession links: {len(fk_issues.get('broken_ship_accession_links', []))}")
    report.append(f"Broken joined ship links: {len(fk_issues.get('broken_joined_ship_links', []))}")
    report.append(f"Broken captain links: {len(fk_issues.get('broken_captain_links', []))}")
    report.append(f"Broken genome-taxonomy links: {len(fk_issues.get('broken_genome_taxonomy_links', []))}")
    report.append("")
    
    # Schema violations
    schema_issues = all_issues.get('schema', {})
    report.append("SCHEMA VIOLATIONS")
    report.append("-" * 50)
    report.append(f"Null/empty violations: {len(schema_issues.get('null_violations', []))}")
    report.append("")
    
    # Summary
    total_issues = sum(
        len(issues) for category in all_issues.values() 
        for issues in category.values() if isinstance(issues, list)
    )
    report.append(f"SUMMARY: {total_issues} total issues found across all categories")
    
    return "\n".join(report)


def main(dry_run: bool = True, output_report: str = None):
    """
    Main function to run comprehensive database cleanup.
    Step 1: Run accession cleanup
    - checks for matches to reverse complemented sequences
    - checks for nested sequences
    - aggregates sequences by accession
    Step 2: Check genome table
    Step 3: Check taxonomic table
    Step 4: Check genomic features
    Step 5: Check foreign keys
    Step 6: Check schema violations
    Step 7: Generate comprehensive report
    
    Args:
        dry_run (bool): If True, only analyze and report without making changes
        output_report (str): Path to save the cleanup report
    """
    logger.info("Starting comprehensive database cleanup process")
    
    all_issues = {}
    
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
    
    logger.info("Step 6: Checking schema violations...")
    all_issues['schema'] = check_schema_violations()
    
    logger.info("Step 7: Generating comprehensive report...")
    report = generate_cleanup_report(all_issues)
    
    if output_report:
        with open(output_report, 'w') as f:
            f.write(report)
        logger.info(f"Report saved to {output_report}")
    else:
        print(report)
    
    logger.info("Comprehensive database cleanup process completed")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Comprehensive database cleanup for Starbase")
    parser.add_argument("--apply", action="store_true", 
                       help="Apply changes to database (default is dry run)")
    parser.add_argument("--report", type=str, 
                       help="Path to save cleanup report")
    
    args = parser.parse_args()
    
    main(dry_run=not args.apply, output_report=args.report)
 
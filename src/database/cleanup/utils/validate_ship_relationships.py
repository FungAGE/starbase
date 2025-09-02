#!/usr/bin/env python3
"""
Validate and correct ship_id relationships using sequence_reference table.

This tool uses the sequence_reference table (populated from FASTA files) as the source
of truth to validate and correct ship_id relationships in the main database tables.
"""

import argparse
import logging
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# Import handling for both module and direct execution
try:
    from ..config.database import StarbaseSession
    from ..config.logging import get_logger
    from ..utils.seq_utils import clean_sequence, revcomp
    from ..utils.classification_utils import generate_md5_hash
except ImportError:
    # Fallback for when running the script directly
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    from src.config.database import StarbaseSession
    from src.config.logging import get_logger
    from src.utils.seq_utils import clean_sequence, revcomp
    from src.utils.classification_utils import generate_md5_hash

logger = get_logger(__name__)


def analyze_ship_id_relationships() -> Dict:
    """
    Analyze the current state of ship_id relationships across all tables.
    
    Returns:
        Dict: Comprehensive analysis of relationship issues
    """
    logger.info("Analyzing ship_id relationships...")
    session = StarbaseSession()
    
    analysis = {
        'sequence_reference_stats': {},
        'ships_table_analysis': {},
        'joined_ships_analysis': {},
        'relationship_issues': {},
        'recommendations': []
    }
    
    try:
        # 1. Analyze sequence_reference table
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM sequence_reference")
            total_refs = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(DISTINCT canonical_md5) FROM sequence_reference")
            unique_md5s = cursor.fetchone()[0]
            cursor.close()
        except Exception as e:
            logger.error(f"Error analyzing sequence_reference table: {str(e)}")
            total_refs = 0
            unique_md5s = 0
        
        analysis['sequence_reference_stats'] = {
            'total_entries': total_refs,
            'unique_sequences': unique_md5s,
            'duplicate_sequences': total_refs - unique_md5s
        }
        
        # 2. Analyze ships table
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM ships")
            total_ships = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM ships WHERE sequence IS NOT NULL AND sequence != ''")
            ships_with_sequences = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM ships WHERE md5 IS NOT NULL AND md5 != ''")
            ships_with_md5 = cursor.fetchone()[0]
            cursor.close()
        except Exception as e:
            logger.error(f"Error analyzing ships table: {str(e)}")
            total_ships = 0
            ships_with_sequences = 0
            ships_with_md5 = 0
        
        analysis['ships_table_analysis'] = {
            'total_ships': total_ships,
            'ships_with_sequences': ships_with_sequences,
            'ships_with_md5': ships_with_md5,
            'ships_without_sequences': total_ships - ships_with_sequences
        }
        
        # 3. Analyze joined_ships table
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM joined_ships")
            total_joined = cursor.fetchone()[0]
            cursor.execute("SELECT COUNT(*) FROM joined_ships WHERE ship_id IS NOT NULL")
            joined_with_ship_id = cursor.fetchone()[0]
            joined_without_ship_id = total_joined - joined_with_ship_id
            
            # Check for duplicate ship_ids (different starshipIDs with same ship_id)
            cursor.execute("""
                SELECT ship_id, COUNT(*) as count, GROUP_CONCAT(starshipID) as starshipIDs
                FROM joined_ships 
                WHERE ship_id IS NOT NULL 
                GROUP BY ship_id 
                HAVING COUNT(*) > 1
            """)
            duplicate_ship_ids = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error analyzing joined_ships table: {str(e)}")
            total_joined = 0
            joined_with_ship_id = 0
            joined_without_ship_id = 0
            duplicate_ship_ids = []
        
        analysis['joined_ships_analysis'] = {
            'total_joined_ships': total_joined,
            'with_ship_id': joined_with_ship_id,
            'without_ship_id': joined_without_ship_id,
            'duplicate_ship_ids': len(duplicate_ship_ids),
            'duplicate_details': [
                {
                    'ship_id': row[0],
                    'count': row[1],
                    'starshipIDs': row[2].split(',') if row[2] else []
                }
                for row in duplicate_ship_ids
            ]
        }
        
        # 4. Identify relationship issues
        issues = []
        
        # Issue 1: Joined_ships entries that should have ship_id but don't
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT js.id, js.starshipID, sr.starshipID as ref_starshipID
                FROM joined_ships js
                LEFT JOIN sequence_reference sr ON js.starshipID = sr.starshipID
                WHERE js.ship_id IS NULL AND sr.starshipID IS NOT NULL
            """)
            missing_ship_ids = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying missing ship_ids: {str(e)}")
            missing_ship_ids = []
        
        if missing_ship_ids:
            issues.append({
                'type': 'missing_ship_id',
                'count': len(missing_ship_ids),
                'description': 'Joined_ships entries that should have ship_id but don\'t',
                'examples': [{'joined_id': row[0], 'starshipID': row[1]} for row in missing_ship_ids[:5]]
            })
        
        # Issue 2: Joined_ships with ship_id that doesn't match sequence_reference
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT js.id, js.starshipID, js.ship_id, s.md5, sr.canonical_md5
                FROM joined_ships js
                JOIN ships s ON js.ship_id = s.id
                JOIN sequence_reference sr ON js.starshipID = sr.starshipID
                WHERE s.md5 != sr.canonical_md5 AND s.rev_comp_md5 != sr.canonical_md5
            """)
            mismatched_ship_ids = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying mismatched ship_ids: {str(e)}")
            mismatched_ship_ids = []
        
        if mismatched_ship_ids:
            issues.append({
                'type': 'mismatched_ship_id',
                'count': len(mismatched_ship_ids),
                'description': 'Joined_ships with ship_id that doesn\'t match sequence_reference',
                'examples': [{'joined_id': row[0], 'starshipID': row[1], 'ship_id': row[2]} for row in mismatched_ship_ids[:5]]
            })
        
        # Issue 3: Ships that exist but aren't referenced by any joined_ships
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT s.id, s.md5
                FROM ships s
                LEFT JOIN joined_ships js ON s.id = js.ship_id
                WHERE js.id IS NULL
            """)
            orphaned_ships = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying orphaned ships: {str(e)}")
            orphaned_ships = []
        
        if orphaned_ships:
            issues.append({
                'type': 'orphaned_ships',
                'count': len(orphaned_ships),
                'description': 'Ships that exist but aren\'t referenced by any joined_ships',
                'examples': [{'ship_id': row[0], 'md5': row[1]} for row in orphaned_ships[:5]]
            })
        
        analysis['relationship_issues'] = issues
        
        # 5. Generate recommendations
        recommendations = []
        
        if len(missing_ship_ids) > 0:
            recommendations.append(f"Fix {len(missing_ship_ids)} joined_ships entries missing ship_id")
        
        if len(mismatched_ship_ids) > 0:
            recommendations.append(f"Fix {len(mismatched_ship_ids)} joined_ships entries with mismatched ship_id")
        
        if len(orphaned_ships) > 0:
            recommendations.append(f"Handle {len(orphaned_ships)} orphaned ships")
        
        if len(duplicate_ship_ids) > 0:
            recommendations.append(f"Resolve {len(duplicate_ship_ids)} cases of duplicate ship_id assignments")
        
        analysis['recommendations'] = recommendations
        
    except Exception as e:
        logger.error(f"Error analyzing relationships: {str(e)}")
        raise
    finally:
        session.close()
    
    return analysis


def fix_ship_id_relationships(dry_run: bool = True) -> Dict:
    """
    Fix ship_id relationships using sequence_reference table as source of truth.
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of fixes applied or would be applied
    """
    logger.info("Fixing ship_id relationships...")
    session = StarbaseSession()
    
    fixes = {
        'missing_ship_ids_fixed': [],
        'mismatched_ship_ids_fixed': [],
        'orphaned_ships_handled': [],
        'duplicate_ship_ids_resolved': [],
        'warnings': []
    }
    
    try:
        # Fix 1: Add missing ship_id to joined_ships entries
        missing_ship_ids = session.execute("""
            SELECT js.id, js.starshipID, sr.canonical_md5
            FROM joined_ships js
            JOIN sequence_reference sr ON js.starshipID = sr.starshipID
            WHERE js.ship_id IS NULL
        """).fetchall()
        
        for row in missing_ship_ids:
            joined_id, starshipID, canonical_md5 = row
            
            # Validate canonical_md5 is not None or empty
            if not canonical_md5 or canonical_md5.strip() == '':
                fixes['warnings'].append({
                    'joined_id': joined_id,
                    'starshipID': starshipID,
                    'issue': f'Empty or None canonical_md5 for starshipID {starshipID}'
                })
                continue
            
            # Find the corresponding ship using raw sqlite3
            try:
                raw_conn = session.connection().connection
                cursor = raw_conn.cursor()
                cursor.execute("""
                    SELECT id FROM ships 
                    WHERE md5 = ? OR rev_comp_md5 = ?
                """, (canonical_md5, canonical_md5))
                ship = cursor.fetchone()
                cursor.close()
            except Exception as e:
                logger.error(f"Error finding ship for {starshipID}: {str(e)}")
                fixes['warnings'].append({
                    'joined_id': joined_id,
                    'starshipID': starshipID,
                    'issue': f'Error finding ship: {str(e)}'
                })
                continue
            
            if ship:
                if not dry_run:
                    try:
                        raw_conn = session.connection().connection
                        cursor = raw_conn.cursor()
                        cursor.execute("""
                            UPDATE joined_ships SET ship_id = ? WHERE id = ?
                        """, (ship[0], joined_id))
                        cursor.close()
                        
                        fixes['missing_ship_ids_fixed'].append({
                            'joined_id': joined_id,
                            'starshipID': starshipID,
                            'ship_id': ship[0],
                            'action': f'Added missing ship_id {ship[0]}'
                        })
                    except Exception as e:
                        logger.error(f"Error updating joined_ships for {starshipID}: {str(e)}")
                        fixes['warnings'].append({
                            'joined_id': joined_id,
                            'starshipID': starshipID,
                            'issue': f'Error updating: {str(e)}'
                        })
                else:
                    fixes['missing_ship_ids_fixed'].append({
                        'joined_id': joined_id,
                        'starshipID': starshipID,
                        'ship_id': ship[0] if ship else 'UNKNOWN',
                        'action': f'Would add missing ship_id {ship[0] if ship else "UNKNOWN"}'
                    })
            else:
                fixes['warnings'].append({
                    'joined_id': joined_id,
                    'starshipID': starshipID,
                    'issue': f'No ship found for canonical_md5 {canonical_md5}'
                })
        
        # Fix 2: Fix mismatched ship_id assignments
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT js.id, js.starshipID, js.ship_id, s.md5, sr.canonical_md5
                FROM joined_ships js
                JOIN ships s ON js.ship_id = s.id
                JOIN sequence_reference sr ON js.starshipID = sr.starshipID
                WHERE s.md5 != sr.canonical_md5 AND s.rev_comp_md5 != sr.canonical_md5
            """)
            mismatched_ship_ids = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying mismatched ship_ids: {str(e)}")
            mismatched_ship_ids = []
        
        for row in mismatched_ship_ids:
            joined_id, starshipID, current_ship_id, current_md5, correct_canonical_md5 = row
            
            # Validate correct_canonical_md5 is not None or empty
            if not correct_canonical_md5 or correct_canonical_md5.strip() == '':
                fixes['warnings'].append({
                    'joined_id': joined_id,
                    'starshipID': starshipID,
                    'issue': f'Empty or None canonical_md5 for starshipID {starshipID}'
                })
                continue
            
            # Find the correct ship using raw sqlite3
            try:
                raw_conn = session.connection().connection
                cursor = raw_conn.cursor()
                cursor.execute("""
                    SELECT id FROM ships 
                    WHERE md5 = ? OR rev_comp_md5 = ?
                """, (correct_canonical_md5, correct_canonical_md5))
                correct_ship = cursor.fetchone()
                cursor.close()
            except Exception as e:
                logger.error(f"Error finding correct ship for {starshipID}: {str(e)}")
                fixes['warnings'].append({
                    'joined_id': joined_id,
                    'starshipID': starshipID,
                    'issue': f'Error finding correct ship: {str(e)}'
                })
                continue
            
            if correct_ship and correct_ship[0] != current_ship_id:
                if not dry_run:
                    try:
                        raw_conn = session.connection().connection
                        cursor = raw_conn.cursor()
                        cursor.execute("""
                            UPDATE joined_ships SET ship_id = ? WHERE id = ?
                        """, (correct_ship[0], joined_id))
                        cursor.close()
                        
                        fixes['mismatched_ship_ids_fixed'].append({
                            'joined_id': joined_id,
                            'starshipID': starshipID,
                            'old_ship_id': current_ship_id,
                            'new_ship_id': correct_ship[0],
                            'action': f'Fixed ship_id from {current_ship_id} to {correct_ship[0]}'
                        })
                    except Exception as e:
                        logger.error(f"Error updating mismatched ship_id for {starshipID}: {str(e)}")
                        fixes['warnings'].append({
                            'joined_id': joined_id,
                            'starshipID': starshipID,
                            'issue': f'Error updating mismatched ship_id: {str(e)}'
                        })
                else:
                    fixes['mismatched_ship_ids_fixed'].append({
                        'joined_id': joined_id,
                        'starshipID': starshipID,
                        'old_ship_id': current_ship_id,
                        'new_ship_id': correct_ship[0],
                        'action': f'Would fix ship_id from {current_ship_id} to {correct_ship[0]}'
                    })
        
        # Fix 3: Handle orphaned ships (ships not referenced by any joined_ships)
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT s.id, s.md5
                FROM ships s
                LEFT JOIN joined_ships js ON s.id = js.ship_id
                WHERE js.id IS NULL
            """)
            orphaned_ships = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying orphaned ships: {str(e)}")
            orphaned_ships = []
        
        for row in orphaned_ships:
            ship_id, md5 = row
            
            # Check if this ship should be referenced by a joined_ships entry
            try:
                raw_conn = session.connection().connection
                cursor = raw_conn.cursor()
                cursor.execute("""
                    SELECT starshipID FROM sequence_reference 
                    WHERE canonical_md5 = ?
                """, (md5,))
                matching_ref = cursor.fetchone()
                cursor.close()
            except Exception as e:
                logger.error(f"Error finding matching reference for ship {ship_id}: {str(e)}")
                fixes['warnings'].append({
                    'ship_id': ship_id,
                    'md5': md5,
                    'issue': f'Error finding matching reference: {str(e)}'
                })
                continue
            
            if matching_ref:
                # This ship should be referenced - find the joined_ships entry
                try:
                    raw_conn = session.connection().connection
                    cursor = raw_conn.cursor()
                    cursor.execute("""
                        SELECT id FROM joined_ships 
                        WHERE starshipID = ? AND ship_id IS NULL
                    """, (matching_ref[0],))
                    joined_entry = cursor.fetchone()
                    cursor.close()
                except Exception as e:
                    logger.error(f"Error finding joined_ships entry for {matching_ref[0]}: {str(e)}")
                    fixes['warnings'].append({
                        'ship_id': ship_id,
                        'md5': md5,
                        'issue': f'Error finding joined_ships entry: {str(e)}'
                    })
                    continue
                
                if joined_entry:
                    if not dry_run:
                        try:
                            raw_conn = session.connection().connection
                            cursor = raw_conn.cursor()
                            cursor.execute("""
                                UPDATE joined_ships SET ship_id = ? WHERE id = ?
                            """, (ship_id, joined_entry[0]))
                            cursor.close()
                            
                            fixes['orphaned_ships_handled'].append({
                                'ship_id': ship_id,
                                'joined_id': joined_entry[0],
                                'starshipID': matching_ref[0],
                                'action': f'Linked orphaned ship {ship_id} to joined_ships {joined_entry[0]}'
                            })
                        except Exception as e:
                            logger.error(f"Error linking orphaned ship {ship_id}: {str(e)}")
                            fixes['warnings'].append({
                                'ship_id': ship_id,
                                'md5': md5,
                                'issue': f'Error linking orphaned ship: {str(e)}'
                            })
                    else:
                        fixes['orphaned_ships_handled'].append({
                            'ship_id': ship_id,
                            'joined_id': joined_entry[0] if joined_entry else 'UNKNOWN',
                            'starshipID': matching_ref[0],
                            'action': f'Would link orphaned ship {ship_id} to joined_ships {joined_entry[0] if joined_entry else "UNKNOWN"}'
                        })
                else:
                    fixes['warnings'].append({
                        'ship_id': ship_id,
                        'md5': md5,
                        'issue': f'Orphaned ship with no matching joined_ships entry for {matching_ref[0]}'
                    })
            else:
                fixes['warnings'].append({
                    'ship_id': ship_id,
                    'md5': md5,
                    'issue': 'Orphaned ship with no matching sequence_reference entry'
                })
        
        # Fix 4: Resolve duplicate ship_id assignments
        try:
            raw_conn = session.connection().connection
            cursor = raw_conn.cursor()
            cursor.execute("""
                SELECT ship_id, COUNT(*) as count, GROUP_CONCAT(starshipID) as starshipIDs
                FROM joined_ships 
                WHERE ship_id IS NOT NULL 
                GROUP BY ship_id 
                HAVING COUNT(*) > 1
            """)
            duplicate_ship_ids = cursor.fetchall()
            cursor.close()
        except Exception as e:
            logger.error(f"Error querying duplicate ship_ids: {str(e)}")
            duplicate_ship_ids = []
        
        for row in duplicate_ship_ids:
            ship_id, count, starshipIDs_str = row
            starshipIDs = starshipIDs_str.split(',') if starshipIDs_str else []
            
            # For each starshipID, verify if it should actually reference this ship_id
            for starshipID in starshipIDs:
                # Check if this starshipID should reference this ship_id
                try:
                    raw_conn = session.connection().connection
                    cursor = raw_conn.cursor()
                    cursor.execute("""
                        SELECT COUNT(*) FROM sequence_reference sr
                        JOIN ships s ON sr.canonical_md5 = s.md5 OR sr.canonical_md5 = s.rev_comp_md5
                        WHERE sr.starshipID = ? AND s.id = ?
                    """, (starshipID, ship_id))
                    should_reference = cursor.fetchone()[0]
                    cursor.close()
                except Exception as e:
                    logger.error(f"Error checking if {starshipID} should reference ship_id {ship_id}: {str(e)}")
                    continue
                
                if not should_reference:
                    # This is a wrong assignment - find the correct ship_id
                    try:
                        raw_conn = session.connection().connection
                        cursor = raw_conn.cursor()
                        cursor.execute("""
                            SELECT s.id FROM sequence_reference sr
                            JOIN ships s ON sr.canonical_md5 = s.md5 OR sr.canonical_md5 = s.rev_comp_md5
                            WHERE sr.starshipID = ?
                        """, (starshipID,))
                        correct_ship = cursor.fetchone()
                        cursor.close()
                    except Exception as e:
                        logger.error(f"Error finding correct ship for {starshipID}: {str(e)}")
                        continue
                    
                    if correct_ship:
                        if not dry_run:
                            try:
                                raw_conn = session.connection().connection
                                cursor = raw_conn.cursor()
                                cursor.execute("""
                                    UPDATE joined_ships SET ship_id = ? WHERE starshipID = ?
                                """, (correct_ship[0], starshipID))
                                cursor.close()
                                
                                fixes['duplicate_ship_ids_resolved'].append({
                                    'starshipID': starshipID,
                                    'old_ship_id': ship_id,
                                    'new_ship_id': correct_ship[0],
                                    'action': f'Fixed duplicate ship_id assignment'
                                })
                            except Exception as e:
                                logger.error(f"Error fixing duplicate ship_id for {starshipID}: {str(e)}")
                                fixes['warnings'].append({
                                    'starshipID': starshipID,
                                    'issue': f'Error fixing duplicate ship_id: {str(e)}'
                                })
                        else:
                            fixes['duplicate_ship_ids_resolved'].append({
                                'starshipID': starshipID,
                                'old_ship_id': ship_id,
                                'new_ship_id': correct_ship[0],
                                'action': f'Would fix duplicate ship_id assignment'
                            })
        
        if not dry_run:
            session.commit()
            logger.info("Ship_id relationships fixed successfully")
        else:
            logger.info("Dry run completed - no changes made")
            
    except Exception as e:
        logger.error(f"Error fixing relationships: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes


def print_analysis_report(analysis: Dict) -> None:
    """Print a formatted analysis report."""
    print("\n" + "="*80)
    print("SHIP_ID RELATIONSHIP ANALYSIS REPORT")
    print("="*80)
    
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


def print_fixes_report(fixes: Dict) -> None:
    """Print a formatted fixes report."""
    print("\n" + "="*80)
    print("SHIP_ID RELATIONSHIP FIXES REPORT")
    print("="*80)
    
    total_fixes = 0
    
    # Missing ship_ids
    missing = fixes['missing_ship_ids_fixed']
    if missing:
        print(f"\n✅ MISSING SHIP_ID FIXES ({len(missing)}):")
        for fix in missing[:5]:  # Show first 5
            print(f"   {fix['starshipID']}: {fix['action']}")
        if len(missing) > 5:
            print(f"   ... and {len(missing) - 5} more")
        total_fixes += len(missing)
    
    # Mismatched ship_ids
    mismatched = fixes['mismatched_ship_ids_fixed']
    if mismatched:
        print(f"\n🔄 MISMATCHED SHIP_ID FIXES ({len(mismatched)}):")
        for fix in mismatched[:5]:  # Show first 5
            print(f"   {fix['starshipID']}: {fix['action']}")
        if len(mismatched) > 5:
            print(f"   ... and {len(mismatched) - 5} more")
        total_fixes += len(mismatched)
    
    # Orphaned ships
    orphaned = fixes['orphaned_ships_handled']
    if orphaned:
        print(f"\n🔗 ORPHANED SHIPS HANDLED ({len(orphaned)}):")
        for fix in orphaned[:5]:  # Show first 5
            print(f"   Ship {fix['ship_id']}: {fix['action']}")
        if len(orphaned) > 5:
            print(f"   ... and {len(orphaned) - 5} more")
        total_fixes += len(orphaned)
    
    # Duplicate ship_ids
    duplicates = fixes['duplicate_ship_ids_resolved']
    if duplicates:
        print(f"\n🔄 DUPLICATE SHIP_ID RESOLVED ({len(duplicates)}):")
        for fix in duplicates[:5]:  # Show first 5
            print(f"   {fix['starshipID']}: {fix['action']}")
        if len(duplicates) > 5:
            print(f"   ... and {len(duplicates) - 5} more")
        total_fixes += len(duplicates)
    
    # Warnings
    warnings = fixes['warnings']
    if warnings:
        print(f"\n⚠️  WARNINGS ({len(warnings)}):")
        for warning in warnings[:5]:  # Show first 5
            print(f"   {warning}")
        if len(warnings) > 5:
            print(f"   ... and {len(warnings) - 5} more")
    
    print(f"\n📈 SUMMARY: {total_fixes} fixes applied")


def main():
    """Main function to run the ship_id relationship validation and correction."""
    parser = argparse.ArgumentParser(
        description="Validate and correct ship_id relationships using sequence_reference table"
    )
    
    parser.add_argument(
        '--analyze',
        action='store_true',
        help='Analyze current relationship state'
    )
    
    parser.add_argument(
        '--fix',
        action='store_true',
        help='Fix relationship issues'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        default=True,
        help='Run in dry-run mode (default: True)'
    )
    
    parser.add_argument(
        '--apply',
        action='store_true',
        help='Apply fixes (overrides --dry-run)'
    )
    
    args = parser.parse_args()
    
    if not args.analyze and not args.fix:
        print("Please specify --analyze and/or --fix")
        return
    
    if args.apply:
        args.dry_run = False
    
    try:
        if args.analyze:
            print("🔍 Analyzing ship_id relationships...")
            analysis = analyze_ship_id_relationships()
            print_analysis_report(analysis)
        
        if args.fix:
            print(f"\n🔧 Fixing ship_id relationships (dry_run: {args.dry_run})...")
            fixes = fix_ship_id_relationships(dry_run=args.dry_run)
            print_fixes_report(fixes)
            
            if args.dry_run and (fixes['missing_ship_ids_fixed'] or fixes['mismatched_ship_ids_fixed'] or 
                                fixes['orphaned_ships_handled'] or fixes['duplicate_ship_ids_resolved']):
                print(f"\n💡 To apply these fixes, run with --apply flag")
        
        print(f"\n✅ Ship_id relationship validation completed successfully")
        
    except Exception as e:
        logger.error(f"Error in ship_id relationship validation: {str(e)}")
        raise


if __name__ == "__main__":
    main()

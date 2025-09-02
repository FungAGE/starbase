import pandas as pd
from sqlalchemy import text, func
from typing import Dict
from datetime import datetime
try:
    from ..config.database import StarbaseSession
    from ..config.logging import get_logger
    from .models.schema import (
        Accessions, Ships, JoinedShips, Genome, Taxonomy, 
        StarshipFeatures, Captains
    )
    from .cleanup_accessions import main as run_accession_cleanup
except ImportError:
    # Fallback for when running the script directly
    import sys
    import os
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    
    from src.config.database import StarbaseSession
    from src.config.logging import get_logger
    from src.database.models.schema import (
        Accessions, Ships, JoinedShips, Genome, Taxonomy, 
        StarshipFeatures, Captains
    )
    from src.database.cleanup_accessions import main as run_accession_cleanup

logger = get_logger(__name__)


def check_ships_accessions_joined_ships_relationships() -> Dict:
    """
    Check the relationships and consistency between ships, accessions, and joined_ships tables.
    
    Based on the expected relationships:
    - joined_ships should contain the most number of records
    - not all entries in joined_ships will be present in ships (only if it has a sequence)
    - not all entries in joined_ships will have an accession_id (only after classification)
    - ships should link to accessions when they have an accession assigned
    
    Returns:
        Dict: Report of relationship issues found
    """
    logger.info("Checking ships-accessions-joined_ships relationships...")
    session = StarbaseSession()
    
    issues = {
        'orphaned_ships': [],
        'orphaned_accessions': [],
        'inconsistent_joined_ships': [],
        'missing_sequence_links': [],
        'orphaned_joined_ships': [],
        'ships_missing_from_joined_ships': [],
        'summary_stats': {}
    }
    
    try:
        # Get basic counts for summary
        total_ships = session.query(Ships).count()
        total_accessions = session.query(Accessions).count()
        total_joined_ships = session.query(JoinedShips).count()
        
        issues['summary_stats'] = {
            'total_ships': total_ships,
            'total_accessions': total_accessions,
            'total_joined_ships': total_joined_ships
        }
        
        logger.info(f"Table counts - Ships: {total_ships}, Accessions: {total_accessions}, JoinedShips: {total_joined_ships}")
        
        # Check 1: Ships that reference non-existent accessions
        orphaned_ships = session.query(Ships).outerjoin(Accessions, Ships.accession_id == Accessions.id).filter(
            (Ships.accession_id.isnot(None)) & (Accessions.id.is_(None))
        ).all()
        
        for ship in orphaned_ships:
            issues['orphaned_ships'].append({
                'ship_id': ship.id,
                'accession_id': ship.accession_id,
                'md5': ship.md5,
                'issue': 'Ship references non-existent accession'
            })
        
        # Check 2: Accessions that have no associated ships or joined_ships
        orphaned_accessions = session.query(Accessions).outerjoin(Ships, Accessions.id == Ships.accession_id).outerjoin(
            JoinedShips, Accessions.id == JoinedShips.ship_id
        ).filter(
            (Ships.id.is_(None)) & (JoinedShips.id.is_(None))
        ).all()
        
        for acc in orphaned_accessions:
            issues['orphaned_accessions'].append({
                'accession_id': acc.id,
                'accession_tag': acc.accession_tag,
                'ship_name': acc.ship_name,
                'issue': 'Accession has no associated ships or joined_ships'
            })
        
        # Check 3: JoinedShips entries that reference non-existent accessions
        orphaned_joined_ships = session.query(JoinedShips).outerjoin(Accessions, JoinedShips.ship_id == Accessions.id).filter(
            (JoinedShips.ship_id.isnot(None)) & (Accessions.id.is_(None))
        ).all()
        
        for joined in orphaned_joined_ships:
            issues['orphaned_joined_ships'].append({
                'joined_id': joined.id,
                'starshipID': joined.starshipID,
                'ship_id': joined.ship_id,
                'issue': 'JoinedShips references non-existent accession'
            })
        
        # Check 4: JoinedShips entries that should have sequences but don't link to ships
        # This is more complex - we need to identify which joined_ships entries should have sequences
        # For now, we'll check for joined_ships with ship_id but no corresponding ship entry
        missing_sequence_links = session.query(JoinedShips).outerjoin(Ships, JoinedShips.ship_id == Ships.accession_id).filter(
            (JoinedShips.ship_id.isnot(None)) & (Ships.id.is_(None))
        ).all()
        
        for joined in missing_sequence_links:
            issues['missing_sequence_links'].append({
                'joined_id': joined.id,
                'starshipID': joined.starshipID,
                'ship_id': joined.ship_id,
                'issue': 'JoinedShips has accession_id but no corresponding ship entry (missing sequence)'
            })
        
        # Check 5: Ships that are missing from joined_ships
        # This ensures all ships have corresponding entries in the main connecting table
        ships_missing_from_joined = session.query(Ships).outerjoin(JoinedShips, Ships.id == JoinedShips.ship_id).filter(
            JoinedShips.id.is_(None)
        ).all()
        
        # Debug: Log the count to see what's happening
        logger.info(f"Found {len(ships_missing_from_joined)} ships missing from joined_ships")
        
        for ship in ships_missing_from_joined:
            issues['ships_missing_from_joined_ships'].append({
                'ship_id': ship.id,
                'md5': ship.md5,
                'accession_id': ship.accession_id,
                'issue': 'Ship exists but has no corresponding entry in joined_ships table'
            })
        
        # Check 6: Inconsistent joined_ships entries
        # Check for joined_ships with missing essential information
        inconsistent_joined = session.query(JoinedShips).filter(
            (JoinedShips.starshipID.is_(None)) | 
            (JoinedShips.starshipID == '') |
            (JoinedShips.evidence.is_(None)) |
            (JoinedShips.evidence == '')
        ).all()
        
        for joined in inconsistent_joined:
            issues['inconsistent_joined_ships'].append({
                'joined_id': joined.id,
                'starshipID': joined.starshipID,
                'evidence': joined.evidence,
                'issue': 'Missing essential information (starshipID or evidence)'
            })
        
        # Log summary of issues found
        logger.info(f"Found {len(issues['orphaned_ships'])} ships with broken accession links")
        logger.info(f"Found {len(issues['orphaned_accessions'])} orphaned accessions")
        logger.info(f"Found {len(issues['orphaned_joined_ships'])} joined_ships with broken accession links")
        logger.info(f"Found {len(issues['missing_sequence_links'])} joined_ships missing sequence links")
        logger.info(f"Found {len(issues['ships_missing_from_joined_ships'])} ships missing from joined_ships")
        logger.info(f"Found {len(issues['inconsistent_joined_ships'])} inconsistent joined_ships entries")
        
    except Exception as e:
        logger.error(f"Error checking ships-accessions-joined_ships relationships: {str(e)}")
        raise
    finally:
        session.close()
    
    return issues


def fix_ships_accessions_joined_ships_relationships(dry_run: bool = True) -> Dict:
    """
    Fix the relationship issues found between ships, accessions, and joined_ships tables.
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed without making changes
        
    Returns:
        Dict: Report of fixes applied or would be applied
    """
    logger.info("Fixing ships-accessions-joined_ships relationships...")
    session = StarbaseSession()
    
    fixes_applied = {
        'ships_fixed': [],
        'accessions_fixed': [],
        'joined_ships_fixed': [],
        'new_joined_ships_created': [],
        'warnings': []
    }
    
    try:
        # Fix 1: Remove broken accession_id references from ships
        orphaned_ships = session.query(Ships).outerjoin(Accessions, Ships.accession_id == Accessions.id).filter(
            (Ships.accession_id.isnot(None)) & (Accessions.id.is_(None))
        ).all()
        
        for ship in orphaned_ships:
            if not dry_run:
                old_accession_id = ship.accession_id
                ship.accession_id = None
                session.add(ship)
                fixes_applied['ships_fixed'].append({
                    'ship_id': ship.id,
                    'action': f'Removed broken accession_id reference: {old_accession_id} -> None'
                })
            else:
                fixes_applied['ships_fixed'].append({
                    'ship_id': ship.id,
                    'action': f'Would remove broken accession_id reference: {ship.accession_id} -> None'
                })
        
        # Fix 2: Remove orphaned accessions that have no associated data
        orphaned_accessions = session.query(Accessions).outerjoin(Ships, Accessions.id == Ships.accession_id).outerjoin(
            JoinedShips, Accessions.id == JoinedShips.ship_id
        ).filter(
            (Ships.id.is_(None)) & (JoinedShips.id.is_(None))
        ).all()
        
        for acc in orphaned_accessions:
            if not dry_run:
                session.delete(acc)
                fixes_applied['accessions_fixed'].append({
                    'accession_id': acc.id,
                    'accession_tag': acc.accession_tag,
                    'action': 'Deleted orphaned accession'
                })
            else:
                fixes_applied['accessions_fixed'].append({
                    'accession_id': acc.id,
                    'accession_tag': acc.accession_tag,
                    'action': 'Would delete orphaned accession'
                })
        
        # Fix 3: Remove broken ship_id references from joined_ships
        orphaned_joined_ships = session.query(JoinedShips).outerjoin(Accessions, JoinedShips.ship_id == Accessions.id).filter(
            (JoinedShips.ship_id.isnot(None)) & (Accessions.id.is_(None))
        ).all()
        
        for joined in orphaned_joined_ships:
            if not dry_run:
                old_ship_id = joined.ship_id
                joined.ship_id = None
                session.add(joined)
                fixes_applied['joined_ships_fixed'].append({
                    'joined_id': joined.id,
                    'starshipID': joined.starshipID,
                    'action': f'Removed broken ship_id reference: {old_ship_id} -> None'
                })
            else:
                fixes_applied['joined_ships_fixed'].append({
                    'joined_id': joined.id,
                    'starshipID': joined.starshipID,
                    'action': f'Would remove broken ship_id reference: {joined.ship_id} -> None'
                })
        
        # Fix 4: Create missing joined_ships entries for ships that don't have them
        # Use a more explicit query to avoid duplicates
        ships_missing_from_joined = session.query(Ships).outerjoin(
            JoinedShips, Ships.id == JoinedShips.ship_id
        ).filter(
            JoinedShips.id.is_(None)
        ).all()
        
        # Debug: Log the count to see what's happening
        logger.info(f"Creating {len(ships_missing_from_joined)} new joined_ships entries")
        
        for ship in ships_missing_from_joined:
            if not dry_run:
                # Double-check that this ship doesn't already have a joined_ships entry
                existing_entry = session.query(JoinedShips).filter(JoinedShips.ship_id == ship.id).first()
                if existing_entry:
                    logger.warning(f"Ship {ship.id} already has joined_ships entry {existing_entry.id}, skipping creation")
                    continue
                
                # Create a new joined_ships entry for this ship
                new_joined_ship = JoinedShips(
                    starshipID=f"SHIP_{ship.id}",
                    evidence="AUTO_CREATED",
                    source="database_cleanup",
                    curated_status="auto_created",
                    ship_id=ship.id,
                    created_at=datetime.now(),
                    updated_at=datetime.now()
                )
                
                session.add(new_joined_ship)
                fixes_applied['new_joined_ships_created'].append({
                    'ship_id': ship.id,
                    'md5': ship.md5,
                    'action': f'Created new joined_ships entry with ID: SHIP_{ship.id}'
                })
            else:
                fixes_applied['new_joined_ships_created'].append({
                    'ship_id': ship.id,
                    'md5': ship.md5,
                    'action': f'Would create new joined_ships entry with ID: SHIP_{ship.id}'
                })
        
        # Fix 5: Clean up inconsistent joined_ships entries
        inconsistent_joined = session.query(JoinedShips).filter(
            (JoinedShips.starshipID.is_(None)) | 
            (JoinedShips.starshipID == '') |
            (JoinedShips.evidence.is_(None)) |
            (JoinedShips.evidence == '')
        ).all()
        
        for joined in inconsistent_joined:
            if not dry_run:
                # Set default values for missing required fields
                if not joined.starshipID or joined.starshipID == '':
                    joined.starshipID = f"UNKNOWN_{joined.id}"
                if not joined.evidence or joined.evidence == '':
                    joined.evidence = "UNKNOWN"
                
                session.add(joined)
                fixes_applied['joined_ships_fixed'].append({
                    'joined_id': joined.id,
                    'starshipID': joined.starshipID,
                    'action': 'Fixed missing required fields with default values'
                })
            else:
                fixes_applied['joined_ships_fixed'].append({
                    'joined_id': joined.id,
                    'starshipID': joined.starshipID,
                    'action': 'Would fix missing required fields with default values'
                })
        
        # Commit changes if not dry run
        if not dry_run:
            session.commit()
            logger.info("Applied fixes to database")
        else:
            logger.info("Dry run - no changes applied to database")
        
        # Log summary of fixes
        logger.info(f"Fixed {len(fixes_applied['ships_fixed'])} ship issues")
        logger.info(f"Fixed {len(fixes_applied['accessions_fixed'])} accession issues")
        logger.info(f"Fixed {len(fixes_applied['joined_ships_fixed'])} joined_ships issues")
        logger.info(f"Created {len(fixes_applied['new_joined_ships_created'])} new joined_ships entries")
        
    except Exception as e:
        logger.error(f"Error fixing ships-accessions-joined_ships relationships: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes_applied


def analyze_table_relationships() -> Dict:
    """
    Perform detailed analysis of the relationships between ships, accessions, and joined_ships tables.
    
    This function provides insights into:
    - How many records in each table have relationships with other tables
    - The distribution of relationship types
    - Potential data quality issues
    
    Returns:
        Dict: Detailed analysis of table relationships
    """
    logger.info("Analyzing table relationships in detail...")
    session = StarbaseSession()
    
    analysis = {
        'relationship_counts': {},
        'data_distribution': {},
        'potential_issues': {},
        'recommendations': []
    }
    
    try:
        # Count total records in each table
        total_ships = session.query(Ships).count()
        total_accessions = session.query(Accessions).count()
        total_joined_ships = session.query(JoinedShips).count()
        
        # Count ships with accession_id
        ships_with_accession = session.query(Ships).filter(Ships.accession_id.isnot(None)).count()
        ships_without_accession = total_ships - ships_with_accession
        
        # Count joined_ships with ship_id (accession_id)
        joined_with_ship = session.query(JoinedShips).filter(JoinedShips.ship_id.isnot(None)).count()
        joined_without_ship = total_joined_ships - joined_with_ship
        
        # Count accessions that are referenced by ships
        accessions_referenced_by_ships = session.query(Accessions).join(Ships, Accessions.id == Ships.accession_id).count()
        
        # Count accessions that are referenced by joined_ships
        accessions_referenced_by_joined = session.query(Accessions).join(JoinedShips, Accessions.id == JoinedShips.ship_id).count()
        
        # Count accessions that are referenced by both
        accessions_referenced_by_both = session.query(Accessions).join(Ships, Accessions.id == Ships.accession_id).join(
            JoinedShips, Accessions.id == JoinedShips.ship_id
        ).count()
        
        analysis['relationship_counts'] = {
            'total_ships': total_ships,
            'total_accessions': total_accessions,
            'total_joined_ships': total_joined_ships,
            'ships_with_accession': ships_with_accession,
            'ships_without_accession': ships_without_accession,
            'joined_with_ship': joined_with_ship,
            'joined_without_ship': joined_without_ship,
            'accessions_referenced_by_ships': accessions_referenced_by_ships,
            'accessions_referenced_by_joined': accessions_referenced_by_joined,
            'accessions_referenced_by_both': accessions_referenced_by_both
        }
        
        # Calculate percentages and ratios
        analysis['data_distribution'] = {
            'ships_with_accession_pct': (ships_with_accession / total_ships * 100) if total_ships > 0 else 0,
            'joined_with_ship_pct': (joined_with_ship / total_joined_ships * 100) if total_joined_ships > 0 else 0,
            'accessions_utilization_pct': (max(accessions_referenced_by_ships, accessions_referenced_by_joined) / total_accessions * 100) if total_accessions > 0 else 0
        }
        
        # Identify potential issues
        potential_issues = []
        
        if ships_without_accession > 0:
            potential_issues.append(f"{ships_without_accession} ships have no accession assigned")
        
        if joined_without_ship > 0:
            potential_issues.append(f"{joined_without_ship} joined_ships entries have no ship/accession link")
        
        if total_accessions > 0 and accessions_referenced_by_ships == 0 and accessions_referenced_by_joined == 0:
            potential_issues.append("All accessions are orphaned (no references from ships or joined_ships)")
        
        analysis['potential_issues'] = potential_issues
        
        # Generate recommendations
        recommendations = []
        
        if ships_without_accession > 0:
            recommendations.append("Consider assigning accessions to ships that have sequences but no accession")
        
        if joined_without_ship > 0:
            recommendations.append("Review joined_ships entries without ship links - these may need accession assignment after classification")
        
        if total_accessions > 0 and accessions_referenced_by_ships == 0 and accessions_referenced_by_joined == 0:
            recommendations.append("Investigate why all accessions are orphaned - this may indicate a data import issue")
        
        analysis['recommendations'] = recommendations
        
        # Log analysis summary
        logger.info(f"Analysis complete:")
        logger.info(f"  - {ships_with_accession}/{total_ships} ships have accessions ({analysis['data_distribution']['ships_with_accession_pct']:.1f}%)")
        logger.info(f"  - {joined_with_ship}/{total_joined_ships} joined_ships have ship links ({analysis['data_distribution']['joined_with_ship_pct']:.1f}%)")
        logger.info(f"  - {len(potential_issues)} potential issues identified")
        logger.info(f"  - {len(recommendations)} recommendations generated")
        
    except Exception as e:
        logger.error(f"Error analyzing table relationships: {str(e)}")
        raise
    finally:
        session.close()
    
    return analysis


def consolidate_duplicate_ships(dry_run: bool = True) -> Dict:
    """
    Consolidate duplicate sequences in the ships table while maintaining foreign key relationships.
    
    This function will:
    1. Find ships with identical sequences (same MD5 or reverse complement MD5)
    2. Keep one ship entry and update all references to point to it
    3. Delete duplicate ship entries
    4. Maintain referential integrity in joined_ships table
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of consolidations made or would be made
    """
    logger.info("Consolidating duplicate sequences in ships table...")
    session = StarbaseSession()
    
    consolidations = {
        'duplicate_groups_found': [],
        'ships_consolidated': [],
        'joined_ships_updated': [],
        'duplicate_ships_removed': [],
        'summary': {}
    }
    
    try:
        # Import the sequence utilities
        try:
            from ..utils.seq_utils import clean_sequence, revcomp
            from ..utils.classification_utils import generate_md5_hash
        except ImportError:
            # Fallback for when running the script directly
            import sys
            import os
            sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
            from src.utils.seq_utils import clean_sequence, revcomp
            from src.utils.classification_utils import generate_md5_hash
        
        # Get all ships with sequences
        all_ships = session.query(Ships).filter(
            (Ships.sequence.isnot(None)) & 
            (Ships.sequence != '')
        ).all()
        
        logger.info(f"Analyzing {len(all_ships)} ships with sequences for duplicates")
        
        # Create a mapping of canonical MD5 -> list of ships
        # Canonical MD5 is the lexicographically smaller of original and reverse complement
        canonical_md5_groups = {}
        
        for ship in all_ships:
            if ship.sequence is None or ship.sequence == '':
                continue
                
            # Clean the sequence
            clean_seq = clean_sequence(ship.sequence)
            if clean_seq is None:
                logger.warning(f"Could not clean sequence for ship {ship.id}, skipping")
                continue
            
            # Generate both MD5 hashes
            md5_hash = generate_md5_hash(clean_seq)
            md5_hash_revcomp = generate_md5_hash(revcomp(clean_seq))
            
            if md5_hash is None or md5_hash_revcomp is None:
                logger.warning(f"Could not generate MD5 for ship {ship.id}, skipping")
                continue
            
            # Use the lexicographically smaller MD5 as the canonical key
            canonical_md5 = min(md5_hash, md5_hash_revcomp)
            
            if canonical_md5 not in canonical_md5_groups:
                canonical_md5_groups[canonical_md5] = []
            canonical_md5_groups[canonical_md5].append(ship)
        
        logger.info(f"Found {len(canonical_md5_groups)} unique sequence groups (considering reverse complements)")
        
        total_consolidated = 0
        total_joined_ships_updated = 0
        total_duplicate_ships_removed = 0
        
        # Process groups with multiple ships (duplicates)
        for canonical_md5, ships in canonical_md5_groups.items():
            if len(ships) == 1:
                continue  # No duplicates in this group
            
            # Sort ships by ID to ensure consistent primary selection
            ships.sort(key=lambda s: s.id)
            primary_ship = ships[0]
            duplicate_ships = ships[1:]
            
            consolidations['duplicate_groups_found'].append({
                'canonical_md5': canonical_md5,
                'primary_ship_id': primary_ship.id,
                'duplicate_ship_ids': [s.id for s in duplicate_ships],
                'count': len(ships),
                'original_md5s': [s.md5 for s in ships]
            })
            
            if not dry_run:
                # Update all joined_ships entries that reference duplicate ships
                for duplicate_ship in duplicate_ships:
                    # Find all joined_ships entries that reference this duplicate ship
                    joined_entries = session.query(JoinedShips).filter(
                        JoinedShips.ship_id == duplicate_ship.id
                    ).all()
                    
                    for entry in joined_entries:
                        old_ship_id = entry.ship_id
                        entry.ship_id = primary_ship.id
                        session.add(entry)
                        total_joined_ships_updated += 1
                        
                        consolidations['joined_ships_updated'].append({
                            'entry_id': entry.id,
                            'starshipID': entry.starshipID,
                            'old_ship_id': old_ship_id,
                            'new_ship_id': primary_ship.id,
                            'action': f'Updated reference from duplicate ship {old_ship_id} to primary ship {primary_ship.id}'
                        })
                    
                    # Delete the duplicate ship
                    session.delete(duplicate_ship)
                    total_duplicate_ships_removed += 1
                    
                    consolidations['duplicate_ships_removed'].append({
                        'ship_id': duplicate_ship.id,
                        'md5': duplicate_ship.md5,
                        'action': f'Deleted duplicate ship, consolidated into primary ship {primary_ship.id}'
                    })
                
                total_consolidated += 1
                
                consolidations['ships_consolidated'].append({
                    'primary_ship_id': primary_ship.id,
                    'canonical_md5': canonical_md5,
                    'duplicates_consolidated': len(duplicate_ships),
                    'action': f'Consolidated {len(duplicate_ships)} duplicate ships into primary ship {primary_ship.id}'
                })
            else:
                # Dry run - just report what would be done
                consolidations['ships_consolidated'].append({
                    'primary_ship_id': primary_ship.id,
                    'canonical_md5': canonical_md5,
                    'duplicates_consolidated': len(duplicate_ships),
                    'action': f'Would consolidate {len(duplicate_ships)} duplicate ships into primary ship {primary_ship.id}'
                })
                
                # Count what would be updated
                for duplicate_ship in duplicate_ships:
                    joined_entries = session.query(JoinedShips).filter(
                        JoinedShips.ship_id == duplicate_ship.id
                    ).all()
                    
                    for entry in joined_entries:
                        consolidations['joined_ships_updated'].append({
                            'entry_id': entry.id,
                            'starshipID': entry.starshipID,
                            'old_ship_id': duplicate_ship.id,
                            'new_ship_id': primary_ship.id,
                            'action': f'Would update reference from duplicate ship {duplicate_ship.id} to primary ship {primary_ship.id}'
                        })
                        total_joined_ships_updated += 1
                    
                    consolidations['duplicate_ships_removed'].append({
                        'ship_id': duplicate_ship.id,
                        'md5': duplicate_ship.md5,
                        'action': f'Would delete duplicate ship, consolidate into primary ship {primary_ship.id}'
                    })
                    total_duplicate_ships_removed += 1
                
                total_consolidated += 1
        
        if not dry_run:
            session.commit()
            logger.info(f"Consolidated {total_consolidated} duplicate groups, updated {total_joined_ships_updated} joined_ships entries, removed {total_duplicate_ships_removed} duplicate ships")
        else:
            logger.info(f"Would consolidate {total_consolidated} duplicate groups, would update {total_joined_ships_updated} joined_ships entries, would remove {total_duplicate_ships_removed} duplicate ships")
        
        consolidations['summary'] = {
            'total_groups_consolidated': total_consolidated,
            'total_joined_ships_updated': total_joined_ships_updated,
            'total_duplicate_ships_removed': total_duplicate_ships_removed,
            'unique_sequence_groups': len(canonical_md5_groups),
            'recommendation': 'Run duplicate detection again after consolidation to verify cleanup'
        }
        
    except Exception as e:
        logger.error(f"Error consolidating duplicate ships: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return consolidations


def analyze_ship_id_mislabeling() -> Dict:
    """
    Analyze the current joined_ships table to detect if ship_id is incorrectly 
    pointing to accessions.id instead of ships.id.
    
    This function will help identify the scope of the mislabeling issue.
    
    Returns:
        Dict: Analysis of the mislabeling issue
    """
    logger.info("Analyzing ship_id mislabeling in joined_ships table...")
    session = StarbaseSession()
    
    analysis = {
        'mislabeling_indicators': [],
        'potential_corrections': [],
        'summary': {}
    }
    
    try:
        # Get all joined_ships entries with ship_id populated
        joined_entries = session.query(JoinedShips).filter(
            JoinedShips.ship_id.isnot(None)
        ).all()
        
        total_entries = len(joined_entries)
        logger.info(f"Analyzing {total_entries} joined_ships entries with ship_id")
        
        # Check 1: How many ship_id values actually exist in ships table
        ship_ids_in_joined = [e.ship_id for e in joined_entries]
        ships_exist = session.query(Ships).filter(Ships.id.in_(ship_ids_in_joined)).count()
        
        # Check 2: How many ship_id values exist in accessions table
        accessions_exist = session.query(Accessions).filter(Accessions.id.in_(ship_ids_in_joined)).count()
        
        # Check 3: Find entries where ship_id points to accessions but should point to ships
        mislabeled_entries = []
        for entry in joined_entries:
            # Check if this ship_id exists in accessions but not in ships
            accession_exists = session.query(Accessions).filter(Accessions.id == entry.ship_id).first()
            ship_exists = session.query(Ships).filter(Ships.id == entry.ship_id).first()
            
            if accession_exists and not ship_exists:
                # This entry is mislabeled - ship_id points to accessions.id
                mislabeled_entries.append({
                    'joined_id': entry.id,
                    'ship_id': entry.ship_id,
                    'starshipID': entry.starshipID,
                    'source': entry.source,
                    'issue': 'ship_id points to accessions.id instead of ships.id'
                })
        
        # Check 4: Find ships that should have joined_ships entries but don't
        ships_without_joined = session.query(Ships).outerjoin(
            JoinedShips, Ships.id == JoinedShips.ship_id
        ).filter(JoinedShips.id.is_(None)).all()
        
        # Check 5: Find accessions that are incorrectly referenced as ship_id
        accessions_incorrectly_referenced = session.query(Accessions).join(
            JoinedShips, Accessions.id == JoinedShips.ship_id
        ).all()
        
        analysis['mislabeling_indicators'] = {
            'total_joined_entries': total_entries,
            'ship_ids_exist_in_ships': ships_exist,
            'ship_ids_exist_in_accessions': accessions_exist,
            'mislabeled_entries': len(mislabeled_entries),
            'ships_without_joined': len(ships_without_joined),
            'accessions_incorrectly_referenced': len(accessions_incorrectly_referenced)
        }
        
        # Generate potential corrections
        for entry in mislabeled_entries[:100]:  # Limit to first 100 for analysis
            # Find the correct ship_id by looking for ships with this accession
            correct_ship = session.query(Ships).filter(
                Ships.accession_id == entry['ship_id']
            ).first()
            
            if correct_ship:
                analysis['potential_corrections'].append({
                    'joined_id': entry['joined_id'],
                    'current_ship_id': entry['ship_id'],
                    'correct_ship_id': correct_ship.id,
                    'starshipID': entry['starshipID'],
                    'action': f'Change ship_id from {entry["ship_id"]} (accession) to {correct_ship.id} (ship)'
                })
        
        analysis['summary'] = {
            'mislabeling_detected': len(mislabeled_entries) > 0,
            'total_corrections_needed': len(analysis['potential_corrections']),
            'recommendation': 'Fix ship_id references to point to ships.id instead of accessions.id'
        }
        
        logger.info(f"Analysis complete:")
        logger.info(f"  - {ships_exist}/{total_entries} ship_id values exist in ships table")
        logger.info(f"  - {accessions_exist}/{total_entries} ship_id values exist in accessions table")
        logger.info(f"  - {len(mislabeled_entries)} entries have mislabeled ship_id")
        logger.info(f"  - {len(ships_without_joined)} ships missing from joined_ships")
        
    except Exception as e:
        logger.error(f"Error analyzing ship_id mislabeling: {str(e)}")
        raise
    finally:
        session.close()
    
    return analysis


def fix_ship_id_mislabeling(dry_run: bool = True) -> Dict:
    """
    Fix the mislabeling of ship_id in joined_ships table by correcting references
    from accessions.id to the correct ships.id.
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of fixes applied or would be applied
    """
    logger.info("Fixing ship_id mislabeling in joined_ships table...")
    session = StarbaseSession()
    
    fixes_applied = {
        'entries_fixed': [],
        'entries_skipped': [],
        'summary': {}
    }
    
    try:
        # Find all mislabeled entries
        joined_entries = session.query(JoinedShips).filter(
            JoinedShips.ship_id.isnot(None)
        ).all()
        
        total_fixed = 0
        total_skipped = 0
        
        for entry in joined_entries:
            # Check if this ship_id exists in accessions but not in ships
            accession_exists = session.query(Accessions).filter(Accessions.id == entry.ship_id).first()
            ship_exists = session.query(Ships).filter(Ships.id == entry.ship_id).first()
            
            if accession_exists and not ship_exists:
                # This entry is mislabeled - find the correct ship_id
                correct_ship = session.query(Ships).filter(
                    Ships.accession_id == entry.ship_id
                ).first()
                
                if correct_ship:
                    if not dry_run:
                        # Fix the mislabeled ship_id
                        old_ship_id = entry.ship_id
                        entry.ship_id = correct_ship.id
                        session.add(entry)
                        
                        fixes_applied['entries_fixed'].append({
                            'joined_id': entry.id,
                            'starshipID': entry.starshipID,
                            'old_ship_id': old_ship_id,
                            'new_ship_id': correct_ship.id,
                            'action': f'Fixed ship_id from {old_ship_id} (accession) to {correct_ship.id} (ship)'
                        })
                        total_fixed += 1
                    else:
                        fixes_applied['entries_fixed'].append({
                            'joined_id': entry.id,
                            'starshipID': entry.starshipID,
                            'old_ship_id': entry.ship_id,
                            'new_ship_id': correct_ship.id if correct_ship else 'UNKNOWN',
                            'action': f'Would fix ship_id from {entry.ship_id} (accession) to {correct_ship.id if correct_ship else "UNKNOWN"} (ship)'
                        })
                        total_fixed += 1
                else:
                    # No corresponding ship found - this is a problem
                    fixes_applied['entries_skipped'].append({
                        'joined_id': entry.id,
                        'starshipID': entry.starshipID,
                        'ship_id': entry.ship_id,
                        'issue': 'No corresponding ship found for this accession'
                    })
                    total_skipped += 1
            else:
                # This entry is correctly labeled
                total_skipped += 1
        
        if not dry_run:
            session.commit()
            logger.info(f"Fixed {total_fixed} mislabeled ship_id references")
        else:
            logger.info(f"Would fix {total_fixed} mislabeled ship_id references")
        
        logger.info(f"Skipped {total_skipped} correctly labeled entries")
        
        fixes_applied['summary'] = {
            'total_fixed': total_fixed,
            'total_skipped': total_skipped,
            'recommendation': 'Run duplicate detection again after fixing mislabeling'
        }
        
    except Exception as e:
        logger.error(f"Error fixing ship_id mislabeling: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes_applied


def identify_duplicate_joined_ships() -> Dict:
    """
    Identify duplicate joined_ships entries that may have been created by previous cleanup runs.
    
    Returns:
        Dict: Report of duplicate entries found
    """
    logger.info("Identifying duplicate joined_ships entries...")
    session = StarbaseSession()
    
    duplicates = {
        'duplicate_ship_ids': [],
        'duplicate_starship_ids': [],
        'auto_created_duplicates': [],
        'summary': {}
    }
    
    try:
        # Find ships with multiple joined_ships entries
        duplicate_ship_entries = session.query(
            JoinedShips.ship_id, 
            func.count(JoinedShips.id).label('count')
        ).filter(
            JoinedShips.ship_id.isnot(None)
        ).group_by(JoinedShips.ship_id).having(
            func.count(JoinedShips.id) > 1
        ).all()
        
        for ship_id, count in duplicate_ship_entries:
            # Get all entries for this ship_id
            entries = session.query(JoinedShips).filter(JoinedShips.ship_id == ship_id).all()
            
            duplicates['duplicate_ship_ids'].append({
                'ship_id': ship_id,
                'count': count,
                'entry_ids': [e.id for e in entries],
                'entries': [
                    {
                        'id': e.id,
                        'starshipID': e.starshipID,
                        'source': e.source,
                        'created_at': e.created_at
                    } for e in entries
                ]
            })
        
        # Find duplicate starshipIDs
        duplicate_starship_entries = session.query(
            JoinedShips.starshipID, 
            func.count(JoinedShips.id).label('count')
        ).filter(
            JoinedShips.starshipID.isnot(None)
        ).group_by(JoinedShips.starshipID).having(
            func.count(JoinedShips.id) > 1
        ).all()
        
        for starship_id, count in duplicate_starship_entries:
            entries = session.query(JoinedShips).filter(JoinedShips.starshipID == starship_id).all()
            duplicates['duplicate_starship_ids'].append({
                'starshipID': starship_id,
                'count': count,
                'entry_ids': [e.id for e in entries]
            })
        
        # Find auto-created duplicates (entries with same source and similar patterns)
        auto_created_duplicates = session.query(
            JoinedShips.source,
            func.count(JoinedShips.id).label('count')
        ).filter(
            JoinedShips.source == 'database_cleanup'
        ).group_by(JoinedShips.source).all()
        
        for source, count in auto_created_duplicates:
            duplicates['auto_created_duplicates'].append({
                'source': source,
                'count': count
            })
        
        # Summary
        duplicates['summary'] = {
            'total_duplicate_ship_ids': len(duplicates['duplicate_ship_ids']),
            'total_duplicate_starship_ids': len(duplicates['duplicate_starship_ids']),
            'total_auto_created': sum(d['count'] for d in duplicates['auto_created_duplicates'])
        }
        
        logger.info(f"Found {duplicates['summary']['total_duplicate_ship_ids']} ships with duplicate joined_ships entries")
        logger.info(f"Found {duplicates['summary']['total_duplicate_starship_ids']} duplicate starshipIDs")
        logger.info(f"Found {duplicates['summary']['total_auto_created']} auto-created entries")
        
    except Exception as e:
        logger.error(f"Error identifying duplicate joined_ships: {str(e)}")
        raise
    finally:
        session.close()
    
    return duplicates


def cleanup_duplicate_joined_ships(dry_run: bool = True) -> Dict:
    """
    Clean up duplicate joined_ships entries, keeping only the most appropriate one for each ship.
    
    Strategy:
    1. For each ship with duplicates, keep the entry with the most complete information
    2. Prefer entries that are NOT from 'database_cleanup' source
    3. Prefer entries with more complete data
    
    Args:
        dry_run (bool): If True, only analyze and report what would be cleaned up
        
    Returns:
        Dict: Report of cleanup actions taken
    """
    logger.info("Cleaning up duplicate joined_ships entries...")
    session = StarbaseSession()
    
    cleanup_report = {
        'duplicates_cleaned': [],
        'entries_removed': [],
        'summary': {}
    }
    
    try:
        # Find ships with duplicate entries
        duplicate_ship_entries = session.query(
            JoinedShips.ship_id, 
            func.count(JoinedShips.id).label('count')
        ).filter(
            JoinedShips.ship_id.isnot(None)
        ).group_by(JoinedShips.ship_id).having(
            func.count(JoinedShips.id) > 1
        ).all()
        
        total_entries_to_remove = 0
        
        for ship_id, count in duplicate_ship_entries:
            # Get all entries for this ship_id
            entries = session.query(JoinedShips).filter(JoinedShips.ship_id == ship_id).order_by(
                JoinedShips.id
            ).all()
            
            if len(entries) <= 1:
                continue
            
            # Score entries to determine which one to keep
            # Higher score = better entry to keep
            scored_entries = []
            for entry in entries:
                score = 0
                
                # Prefer entries NOT from database_cleanup
                if entry.source != 'database_cleanup':
                    score += 100
                
                # Prefer entries with more complete data
                if entry.evidence and entry.evidence != '':
                    score += 10
                if entry.starshipID and entry.starshipID != '':
                    score += 10
                if entry.curated_status and entry.curated_status != '':
                    score += 5
                
                # Prefer older entries (more likely to be original)
                if entry.created_at:
                    score += 1
                
                scored_entries.append((entry, score))
            
            # Sort by score (highest first) and keep the best one
            scored_entries.sort(key=lambda x: x[1], reverse=True)
            entry_to_keep = scored_entries[0][0]
            entries_to_remove = [e for e, _ in scored_entries[1:]]
            
            if not dry_run:
                # Remove duplicate entries
                for entry in entries_to_remove:
                    session.delete(entry)
                    cleanup_report['entries_removed'].append({
                        'entry_id': entry.id,
                        'ship_id': entry.ship_id,
                        'starshipID': entry.starshipID,
                        'source': entry.source,
                        'action': 'Removed duplicate entry'
                    })
                
                total_entries_to_remove += len(entries_to_remove)
            else:
                # Just report what would be done
                for entry in entries_to_remove:
                    cleanup_report['entries_removed'].append({
                        'entry_id': entry.id,
                        'ship_id': entry.ship_id,
                        'starshipID': entry.starshipID,
                        'source': entry.source,
                        'action': 'Would remove duplicate entry'
                    })
                
                total_entries_to_remove += len(entries_to_remove)
            
            cleanup_report['duplicates_cleaned'].append({
                'ship_id': ship_id,
                'entries_found': len(entries),
                'entry_kept': entry_to_keep.id,
                'entries_removed': len(entries_to_remove)
            })
        
        if not dry_run:
            session.commit()
            logger.info(f"Removed {total_entries_to_remove} duplicate joined_ships entries")
        else:
            logger.info(f"Would remove {total_entries_to_remove} duplicate joined_ships entries")
        
        cleanup_report['summary'] = {
            'ships_with_duplicates': len(cleanup_report['duplicates_cleaned']),
            'total_entries_removed': total_entries_to_remove
        }
        
    except Exception as e:
        logger.error(f"Error cleaning up duplicate joined_ships: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return cleanup_report


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
    
    # Ships-Accessions-JoinedShips relationship issues
    relationship_issues = all_issues.get('ships_accessions_joined_ships', {})
    report.append("SHIPS-ACCESSIONS-JOINED_SHIPS RELATIONSHIP ISSUES")
    report.append("-" * 60)
    if 'summary_stats' in relationship_issues:
        stats = relationship_issues['summary_stats']
        report.append(f"Table counts - Ships: {stats.get('total_ships', 0)}, Accessions: {stats.get('total_accessions', 0)}, JoinedShips: {stats.get('total_joined_ships', 0)}")
        report.append("")
    report.append(f"Orphaned ships: {len(relationship_issues.get('orphaned_ships', []))}")
    report.append(f"Orphaned accessions: {len(relationship_issues.get('orphaned_accessions', []))}")
    report.append(f"Orphaned joined_ships: {len(relationship_issues.get('orphaned_joined_ships', []))}")
    report.append(f"Missing sequence links: {len(relationship_issues.get('missing_sequence_links', []))}")
    report.append(f"Ships missing from joined_ships: {len(relationship_issues.get('ships_missing_from_joined_ships', []))}")
    report.append(f"Inconsistent joined_ships: {len(relationship_issues.get('inconsistent_joined_ships', []))}")
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
    
    # Ship_id relationship issues
    ship_id_issues = all_issues.get('ship_id_relationships', {})
    report.append("SHIP_ID RELATIONSHIP ISSUES")
    report.append("-" * 50)
    if 'sequence_reference_stats' in ship_id_issues:
        stats = ship_id_issues['sequence_reference_stats']
        report.append(f"Sequence reference entries: {stats.get('total_entries', 0)}")
        report.append(f"Unique sequences: {stats.get('unique_sequences', 0)}")
        report.append("")
    if 'joined_ships_analysis' in ship_id_issues:
        joined_stats = ship_id_issues['joined_ships_analysis']
        report.append(f"Joined_ships with ship_id: {joined_stats.get('with_ship_id', 0)}")
        report.append(f"Joined_ships without ship_id: {joined_stats.get('without_ship_id', 0)}")
        report.append(f"Duplicate ship_id cases: {joined_stats.get('duplicate_ship_ids', 0)}")
        report.append("")
    if 'relationship_issues' in ship_id_issues:
        for issue in ship_id_issues['relationship_issues']:
            report.append(f"{issue['type']}: {issue['count']} - {issue['description']}")
    if 'recommendations' in ship_id_issues:
        for rec in ship_id_issues['recommendations']:
            report.append(f"  - {rec}")
    report.append("")
    
    # Summary
    total_issues = sum(
        len(issues) for category in all_issues.values() 
        for issues in category.values() if isinstance(issues, list)
    )
    report.append(f"SUMMARY: {total_issues} total issues found across all categories")
    
    return "\n".join(report)
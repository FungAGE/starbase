import pandas as pd
from sqlalchemy import text, func
from typing import Dict, List
import json
from datetime import datetime
import csv
import time
import urllib.parse
import urllib.request
import re
import os
try:
    from ...config.database import StarbaseSession
    from ...config.logging import get_logger
    from ...database.models.schema import (
        Accessions, Ships, JoinedShips, Genome, Taxonomy, 
        StarshipFeatures, Captains, FamilyNames, Navis, Haplotype
    )
    from .cleanup_accessions import main as run_accession_cleanup
except ImportError:
    # Fallback for when running the script directly
    import sys
    import os
    PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    if PROJECT_ROOT not in sys.path:
        sys.path.append(PROJECT_ROOT)
    
    from src.config.database import StarbaseSession
    from src.config.logging import get_logger
    from src.database.models.schema import (
        Accessions, Ships, JoinedShips, Genome, Taxonomy, 
        StarshipFeatures, Captains, FamilyNames, Navis, Haplotype
    )
    from src.database.cleanup.utils.cleanup_accessions import main as run_accession_cleanup

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
        
        # Check 3: JoinedShips entries that reference non-existent ships
        orphaned_joined_ships = session.query(JoinedShips).outerjoin(Ships, JoinedShips.ship_id == Ships.id).filter(
            (JoinedShips.ship_id.isnot(None)) & (Ships.id.is_(None))
        ).all()
        
        for joined in orphaned_joined_ships:
            issues['orphaned_joined_ships'].append({
                'joined_id': joined.id,
                'starshipID': joined.starshipID,
                'ship_id': joined.ship_id,
                'issue': 'JoinedShips references non-existent ship'
            })
        
        # Check 4: JoinedShips entries that reference ships with no sequences
        # Check for joined_ships with ship_id but the referenced ship has no sequence
        missing_sequence_links = session.query(JoinedShips).join(Ships, JoinedShips.ship_id == Ships.id).filter(
            (JoinedShips.ship_id.isnot(None)) &
            ((Ships.sequence.is_(None)) | (Ships.sequence == ''))
        ).all()
        
        for joined in missing_sequence_links:
            issues['missing_sequence_links'].append({
                'joined_id': joined.id,
                'starshipID': joined.starshipID,
                'ship_id': joined.ship_id,
                'issue': 'JoinedShips has ship_id but the referenced ship has no sequence'
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
        orphaned_joined_ships = session.query(JoinedShips).outerjoin(Ships, JoinedShips.ship_id == Ships.id).filter(
            (JoinedShips.ship_id.isnot(None)) & (Ships.id.is_(None))
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
                
                # Determine proper starshipID based on available information
                starship_id = f"SHIP_{ship.id}"  # fallback
                evidence = "AUTO_CREATED"
                
                # If ship has an accession, use the accession information
                if ship.accession_id:
                    accession = session.query(Accessions).filter(Accessions.id == ship.accession_id).first()
                    if accession and accession.ship_name:
                        starship_id = accession.ship_name
                        evidence = "AUTO_CREATED_FROM_ACCESSION"
                    elif accession and accession.accession_tag:
                        starship_id = accession.accession_tag
                        evidence = "AUTO_CREATED_FROM_ACCESSION_TAG"
                
                # Create a new joined_ships entry for this ship
                new_joined_ship = JoinedShips(
                    starshipID=starship_id,
                    evidence=evidence,
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
                    'starshipID': starship_id,
                    'evidence': evidence,
                    'action': f'Created new joined_ships entry with starshipID: {starship_id}'
                })
            else:
                # For dry run, also calculate what the starshipID would be
                starship_id = f"SHIP_{ship.id}"  # fallback
                evidence = "AUTO_CREATED"
                
                # If ship has an accession, use the accession information
                if ship.accession_id:
                    accession = session.query(Accessions).filter(Accessions.id == ship.accession_id).first()
                    if accession and accession.ship_name:
                        starship_id = accession.ship_name
                        evidence = "AUTO_CREATED_FROM_ACCESSION"
                    elif accession and accession.accession_tag:
                        starship_id = accession.accession_tag
                        evidence = "AUTO_CREATED_FROM_ACCESSION_TAG"
                
                fixes_applied['new_joined_ships_created'].append({
                    'ship_id': ship.id,
                    'md5': ship.md5,
                    'starshipID': starship_id,
                    'evidence': evidence,
                    'action': f'Would create new joined_ships entry with starshipID: {starship_id}'
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


def fix_placeholder_starship_ids(dry_run: bool = True) -> Dict:
    """
    Fix joined_ships entries that have placeholder starshipID patterns like 'SHIP_XXX'
    and replace them with proper accession-based starshipIDs.
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of fixes applied or would be applied
    """
    logger.info("Fixing placeholder starshipID entries...")
    session = StarbaseSession()
    
    fixes_applied = {
        'placeholder_ids_fixed': [],
        'no_accession_found': [],
        'warnings': [],
        'summary': {}
    }
    
    try:
        # Find joined_ships entries with placeholder starshipID patterns
        placeholder_entries = session.query(JoinedShips).filter(
            JoinedShips.starshipID.like('SHIP_%')
        ).all()
        
        logger.info(f"Found {len(placeholder_entries)} entries with placeholder starshipID patterns")
        
        for entry in placeholder_entries:
            # Get the associated ship
            ship = session.query(Ships).filter(Ships.id == entry.ship_id).first()
            if not ship:
                fixes_applied['warnings'].append({
                    'joined_id': entry.id,
                    'starshipID': entry.starshipID,
                    'issue': 'No associated ship found'
                })
                continue
            
            # Check if the ship has an accession
            if not ship.accession_id:
                fixes_applied['no_accession_found'].append({
                    'joined_id': entry.id,
                    'starshipID': entry.starshipID,
                    'ship_id': ship.id,
                    'issue': 'Ship has no accession_id'
                })
                continue
            
            # Get the accession information
            accession = session.query(Accessions).filter(Accessions.id == ship.accession_id).first()
            if not accession:
                fixes_applied['warnings'].append({
                    'joined_id': entry.id,
                    'starshipID': entry.starshipID,
                    'ship_id': ship.id,
                    'accession_id': ship.accession_id,
                    'issue': 'Accession not found'
                })
                continue
            
            # Determine the proper starshipID
            new_starship_id = None
            new_evidence = None
            
            if accession.ship_name:
                new_starship_id = accession.ship_name
                new_evidence = "FIXED_FROM_ACCESSION"
            elif accession.accession_tag:
                new_starship_id = accession.accession_tag
                new_evidence = "FIXED_FROM_ACCESSION_TAG"
            
            if new_starship_id:
                if not dry_run:
                    old_starship_id = entry.starshipID
                    entry.starshipID = new_starship_id
                    entry.evidence = new_evidence
                    entry.updated_at = datetime.now()
                    session.add(entry)
                    
                    fixes_applied['placeholder_ids_fixed'].append({
                        'joined_id': entry.id,
                        'old_starshipID': old_starship_id,
                        'new_starshipID': new_starship_id,
                        'evidence': new_evidence,
                        'ship_id': ship.id,
                        'accession_tag': accession.accession_tag,
                        'action': f'Fixed placeholder {old_starship_id} -> {new_starship_id}'
                    })
                else:
                    fixes_applied['placeholder_ids_fixed'].append({
                        'joined_id': entry.id,
                        'old_starshipID': entry.starshipID,
                        'new_starshipID': new_starship_id,
                        'evidence': new_evidence,
                        'ship_id': ship.id,
                        'accession_tag': accession.accession_tag,
                        'action': f'Would fix placeholder {entry.starshipID} -> {new_starship_id}'
                    })
            else:
                fixes_applied['no_accession_found'].append({
                    'joined_id': entry.id,
                    'starshipID': entry.starshipID,
                    'ship_id': ship.id,
                    'accession_id': ship.accession_id,
                    'issue': 'Accession has no ship_name or accession_tag'
                })
        
        if not dry_run:
            session.commit()
            logger.info(f"Fixed {len(fixes_applied['placeholder_ids_fixed'])} placeholder starshipID entries")
        else:
            logger.info(f"Would fix {len(fixes_applied['placeholder_ids_fixed'])} placeholder starshipID entries")
        
        fixes_applied['summary'] = {
            'total_placeholder_entries': len(placeholder_entries),
            'fixed_count': len(fixes_applied['placeholder_ids_fixed']),
            'no_accession_count': len(fixes_applied['no_accession_found']),
            'warnings_count': len(fixes_applied['warnings']),
            'recommendation': 'Run again after fixing accession assignments for remaining entries'
        }
        
        logger.info(f"Found {len(placeholder_entries)} placeholder entries")
        logger.info(f"Fixed {len(fixes_applied['placeholder_ids_fixed'])} entries")
        logger.info(f"Could not fix {len(fixes_applied['no_accession_found'])} entries (no accession)")
        logger.info(f"Warnings: {len(fixes_applied['warnings'])}")
        
    except Exception as e:
        logger.error(f"Error fixing placeholder starshipID entries: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes_applied


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


def enable_foreign_keys(dry_run: bool = True) -> Dict:
    """
    Enable foreign key constraints in SQLite database.

    Args:
        dry_run (bool): If True, only report what would be done

    Returns:
        Dict: Report of the operation
    """
    logger.info("Enabling foreign key constraints...")
    session = StarbaseSession()

    report = {
        'action': 'enable_foreign_keys',
        'status': 'completed' if not dry_run else 'dry_run',
        'details': []
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Check current foreign key status
        cursor.execute("PRAGMA foreign_keys")
        current_status = cursor.fetchone()[0]
        report['details'].append(f"Current foreign key status: {'ON' if current_status else 'OFF'}")

        if not dry_run:
            cursor.execute("PRAGMA foreign_keys = ON")
            cursor.close()
            session.commit()
            report['details'].append("Foreign key constraints enabled")
        else:
            report['details'].append("Would enable foreign key constraints")

        cursor.close()

    except Exception as e:
        logger.error(f"Error enabling foreign keys: {str(e)}")
        report['status'] = 'error'
        report['details'].append(f"Error: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()

    return report


def disable_foreign_keys(dry_run: bool = True) -> Dict:
    """
    Disable foreign key constraints in SQLite database.

    Args:
        dry_run (bool): If True, only report what would be done

    Returns:
        Dict: Report of the operation
    """
    logger.info("Disabling foreign key constraints...")
    session = StarbaseSession()

    report = {
        'action': 'disable_foreign_keys',
        'status': 'completed' if not dry_run else 'dry_run',
        'details': []
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Check current foreign key status
        cursor.execute("PRAGMA foreign_keys")
        current_status = cursor.fetchone()[0]
        report['details'].append(f"Current foreign key status: {'ON' if current_status else 'OFF'}")

        if not dry_run:
            cursor.execute("PRAGMA foreign_keys = OFF")
            cursor.close()
            session.commit()
            report['details'].append("Foreign key constraints disabled")
        else:
            report['details'].append("Would disable foreign key constraints")

        cursor.close()

    except Exception as e:
        logger.error(f"Error disabling foreign keys: {str(e)}")
        report['status'] = 'error'
        report['details'].append(f"Error: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()

    return report


def check_foreign_key_enforcement() -> Dict:
    """
    Check if foreign key constraints are actually being enforced in the database.

    Returns:
        Dict: Report of foreign key enforcement status
    """
    logger.info("Checking foreign key enforcement status...")
    session = StarbaseSession()

    report = {
        'foreign_keys_enabled': False,
        'tables_with_constraints': [],
        'missing_constraints': [],
        'recommendations': []
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Check if foreign keys are enabled
        cursor.execute("PRAGMA foreign_keys")
        fk_enabled = cursor.fetchone()[0]
        report['foreign_keys_enabled'] = bool(fk_enabled)

        # Get all tables
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall()]

        # Check each table for foreign key constraints
        for table in tables:
            cursor.execute(f"PRAGMA foreign_key_list({table})")
            fk_list = cursor.fetchall()

            if fk_list:
                report['tables_with_constraints'].append({
                    'table': table,
                    'constraints': len(fk_list),
                    'details': [
                        {
                            'from_col': row[3],
                            'to_table': row[2],
                            'to_col': row[4]
                        }
                        for row in fk_list
                    ]
                })
            else:
                # Check if table should have foreign keys based on schema
                should_have_fks = _table_should_have_foreign_keys(table)
                if should_have_fks:
                    report['missing_constraints'].append({
                        'table': table,
                        'expected_constraints': should_have_fks
                    })

        cursor.close()

        # Generate recommendations
        if not report['foreign_keys_enabled']:
            report['recommendations'].append("Enable foreign key constraints with PRAGMA foreign_keys = ON")

        if report['missing_constraints']:
            report['recommendations'].append(f"Add foreign key constraints to {len(report['missing_constraints'])} tables")

        if not report['tables_with_constraints']:
            report['recommendations'].append("Database has no foreign key constraints - consider adding them for data integrity")

        logger.info(f"Foreign keys enabled: {report['foreign_keys_enabled']}")
        logger.info(f"Tables with FK constraints: {len(report['tables_with_constraints'])}")
        logger.info(f"Tables missing FK constraints: {len(report['missing_constraints'])}")

    except Exception as e:
        logger.error(f"Error checking foreign key enforcement: {str(e)}")
        raise
    finally:
        session.close()

    return report


def _table_should_have_foreign_keys(table_name: str) -> List[Dict]:
    """
    Determine what foreign key constraints a table should have based on the schema.

    Args:
        table_name (str): Name of the table

    Returns:
        List[Dict]: List of expected foreign key constraints
    """
    expected_constraints = {
        'ships': [
            {'from_col': 'accession_id', 'to_table': 'accessions', 'to_col': 'id'}
        ],
        'captains': [
            {'from_col': 'ship_id', 'to_table': 'ships', 'to_col': 'id'}
        ],
        'genomes': [
            {'from_col': 'taxonomy_id', 'to_table': 'taxonomy', 'to_col': 'id'}
        ],
        'starship_features': [
            {'from_col': 'ship_id', 'to_table': 'ships', 'to_col': 'id'},
            {'from_col': 'captain_id', 'to_table': 'captains', 'to_col': 'id'}
        ],
        'navis_names': [
            {'from_col': 'ship_family_id', 'to_table': 'family_names', 'to_col': 'id'}
        ],
        'haplotype_names': [
            {'from_col': 'navis_id', 'to_table': 'navis_names', 'to_col': 'id'},
            {'from_col': 'ship_family_id', 'to_table': 'family_names', 'to_col': 'id'}
        ],
        'gff': [
            {'from_col': 'accession_id', 'to_table': 'accessions', 'to_col': 'id'},
            {'from_col': 'ship_id', 'to_table': 'ships', 'to_col': 'id'}
        ],
        'joined_ships': [
            {'from_col': 'ship_family_id', 'to_table': 'family_names', 'to_col': 'id'},
            {'from_col': 'tax_id', 'to_table': 'taxonomy', 'to_col': 'id'},
            {'from_col': 'ship_id', 'to_table': 'ships', 'to_col': 'id'},
            {'from_col': 'genome_id', 'to_table': 'genomes', 'to_col': 'id'},
            {'from_col': 'captain_id', 'to_table': 'captains', 'to_col': 'id'},
            {'from_col': 'ship_navis_id', 'to_table': 'navis_names', 'to_col': 'id'},
            {'from_col': 'ship_haplotype_id', 'to_table': 'haplotype_names', 'to_col': 'id'}
        ]
    }

    return expected_constraints.get(table_name, [])


def fix_joined_ships_ship_id_foreign_key(dry_run: bool = True) -> Dict:
    """
    CAUTION: This function fixes cases where joined_ships.ship_id contains accession IDs
    instead of ship IDs, using conservative validation.

    VALIDATION APPROACH:
    1. Check if ship_id is not a valid ship ID (doesn't exist in ships table)
    2. Verify it's a valid accession ID (exists in accessions table)
    3. Map to corresponding ship ID (accession -> ship relationship)
    4. Update joined_ships record

    VALIDATION LEVELS:
    - ✅ VALIDATED: Accession exists and has corresponding ship
    - ❌ INVALID: ship_id contains ID not found in ships or accessions tables
    - ⚠️  NO SHIP: Accession exists but no corresponding ship found

    WARNING: This function assumes that invalid ship_ids in joined_ships are actually
    accession IDs that should be mapped to ship IDs. Use --dry-run first to verify
    the mappings look correct.

    Args:
        dry_run (bool): If True, only report what would be fixed without making changes

    Returns:
        Dict: Report of fixes applied with validation status
    """
    logger.info("CAUTION: Attempting to fix joined_ships.ship_id foreign key relationship...")
    logger.warning("This function assumes ship_id contains accession IDs - verify with --dry-run first!")

    session = StarbaseSession()

    fixes = {
        'ship_ids_fixed': [],
        'no_ship_for_accession': [],
        'invalid_references_found': [],
        'summary': {}
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Step 1: Get all joined_ships entries with non-null ship_id
        cursor.execute("""
            SELECT js.id, js.starshipID, js.ship_id
            FROM joined_ships js
            WHERE js.ship_id IS NOT NULL
        """)
        joined_entries = cursor.fetchall()

        logger.info(f"Checking {len(joined_entries)} joined_ships entries for potential ship_id corrections")
        logger.info("Performing semantic validation based on starshipID and accession data...")

        fixed_count = 0

        for js_id, starshipID, current_ship_id in joined_entries:
            # First, check if current_ship_id is already a valid ship ID
            cursor.execute("SELECT id FROM ships WHERE id = ?", (current_ship_id,))
            ship_exists = cursor.fetchone()

            if ship_exists:
                # ship_id is already valid, skip this entry
                continue

            # ship_id is not a valid ship ID, check if it's an accession ID
            cursor.execute("""
                SELECT a.id, a.accession_tag, s.id as ship_id
                FROM accessions a
                LEFT JOIN ships s ON s.accession_id = a.id
                WHERE a.id = ?
            """, (current_ship_id,))
            accession_data = cursor.fetchone()

            if not accession_data:
                # current_ship_id is neither a valid ship ID nor a valid accession ID
                fixes['invalid_references_found'].append({
                    'joined_id': js_id,
                    'starshipID': starshipID,
                    'invalid_ship_id': current_ship_id,
                    'issue': 'ship_id contains invalid ID (not in ships or accessions tables)'
                })
                continue

            acc_id, acc_accession_tag, corresponding_ship_id = accession_data

            if corresponding_ship_id:
                # Accession has a corresponding ship - safe to fix
                if not dry_run:
                    cursor.execute(
                        "UPDATE joined_ships SET ship_id = ? WHERE id = ?",
                        (corresponding_ship_id, js_id)
                    )
                fixes['ship_ids_fixed'].append({
                    'joined_id': js_id,
                    'starshipID': starshipID,
                    'old_ship_id': current_ship_id,
                    'new_ship_id': corresponding_ship_id,
                    'accession_tag': acc_accession_tag,
                    'action': 'Fixed ship_id from accession to ship reference'
                })
                fixed_count += 1
            else:
                # Accession exists but no corresponding ship
                fixes['no_ship_for_accession'].append({
                    'joined_id': js_id,
                    'starshipID': starshipID,
                    'accession_id': current_ship_id,
                    'accession_tag': acc_accession_tag,
                    'issue': 'Accession exists but no corresponding ship record found'
                })

        cursor.close()

        if not dry_run:
            session.commit()
            logger.info("Applied joined_ships.ship_id foreign key fixes")
        else:
            logger.info("Dry run completed - review validation results before applying")

        # Summary
        fixes['summary'] = {
            'ship_ids_fixed': fixed_count,
            'invalid_references': len(fixes['invalid_references_found']),
            'no_ship_for_accession': len(fixes['no_ship_for_accession']),
            'recommendation': 'Review results before applying - ensure mappings look correct',
            'warning': 'This function maps accession IDs to ship IDs - verify the data relationships are correct!'
        }

        logger.info(f"Ship IDs fixed: {fixed_count}")
        logger.info(f"Invalid references found: {len(fixes['invalid_references_found'])}")
        logger.info(f"Accessions without ships: {len(fixes['no_ship_for_accession'])}")

        if fixes['invalid_references_found']:
            logger.warning("⚠️  INVALID REFERENCES FOUND - These entries have corrupted ship_id values!")
            for invalid in fixes['invalid_references_found'][:3]:  # Show first 3
                logger.warning(f"   Joined {invalid['joined_id']} ({invalid['starshipID']}): invalid ship_id {invalid['invalid_ship_id']}")

    except Exception as e:
        logger.error(f"Error fixing joined_ships.ship_id foreign key: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()

    return fixes


def validate_all_foreign_keys() -> Dict:
    """
    Comprehensive validation of all foreign key relationships in the database.

    This function checks every foreign key relationship defined in the schema
    and reports any violations or orphaned records.

    Returns:
        Dict: Comprehensive report of foreign key validation results
    """
    logger.info("Performing comprehensive foreign key validation...")
    session = StarbaseSession()

    validation_report = {
        'foreign_keys_enabled': False,
        'total_relationships_checked': 0,
        'violations_found': 0,
        'relationships': {},
        'summary': {},
        'recommendations': []
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Check if foreign keys are enabled
        cursor.execute("PRAGMA foreign_keys")
        validation_report['foreign_keys_enabled'] = bool(cursor.fetchone()[0])

        # Define all foreign key relationships to check
        fk_relationships = {
            'ships_to_accessions': {
                'from_table': 'ships',
                'from_col': 'accession_id',
                'to_table': 'accessions',
                'to_col': 'id'
            },
            'captains_to_ships': {
                'from_table': 'captains',
                'from_col': 'ship_id',
                'to_table': 'ships',
                'to_col': 'id'
            },
            'genomes_to_taxonomy': {
                'from_table': 'genomes',
                'from_col': 'taxonomy_id',
                'to_table': 'taxonomy',
                'to_col': 'id'
            },
            'starship_features_to_ships': {
                'from_table': 'starship_features',
                'from_col': 'ship_id',
                'to_table': 'ships',
                'to_col': 'id'
            },
            'starship_features_to_captains': {
                'from_table': 'starship_features',
                'from_col': 'captain_id',
                'to_table': 'captains',
                'to_col': 'id'
            },
            'navis_names_to_family_names': {
                'from_table': 'navis_names',
                'from_col': 'ship_family_id',
                'to_table': 'family_names',
                'to_col': 'id'
            },
            'haplotype_names_to_navis': {
                'from_table': 'haplotype_names',
                'from_col': 'navis_id',
                'to_table': 'navis_names',
                'to_col': 'id'
            },
            'haplotype_names_to_family_names': {
                'from_table': 'haplotype_names',
                'from_col': 'ship_family_id',
                'to_table': 'family_names',
                'to_col': 'id'
            },
            'gff_to_accessions': {
                'from_table': 'gff',
                'from_col': 'accession_id',
                'to_table': 'accessions',
                'to_col': 'id'
            },
            'gff_to_ships': {
                'from_table': 'gff',
                'from_col': 'ship_id',
                'to_table': 'ships',
                'to_col': 'id'
            },
            'joined_ships_to_family_names': {
                'from_table': 'joined_ships',
                'from_col': 'ship_family_id',
                'to_table': 'family_names',
                'to_col': 'id'
            },
            'joined_ships_to_taxonomy': {
                'from_table': 'joined_ships',
                'from_col': 'tax_id',
                'to_table': 'taxonomy',
                'to_col': 'id'
            },
            'joined_ships_to_ships': {
                'from_table': 'joined_ships',
                'from_col': 'ship_id',
                'to_table': 'ships',
                'to_col': 'id'
            },
            'joined_ships_to_genomes': {
                'from_table': 'joined_ships',
                'from_col': 'genome_id',
                'to_table': 'genomes',
                'to_col': 'id'
            },
            'joined_ships_to_captains': {
                'from_table': 'joined_ships',
                'from_col': 'captain_id',
                'to_table': 'captains',
                'to_col': 'id'
            },
            'joined_ships_to_navis_names': {
                'from_table': 'joined_ships',
                'from_col': 'ship_navis_id',
                'to_table': 'navis_names',
                'to_col': 'id'
            },
            'joined_ships_to_haplotype_names': {
                'from_table': 'joined_ships',
                'from_col': 'ship_haplotype_id',
                'to_table': 'haplotype_names',
                'to_col': 'id'
            }
        }

        total_violations = 0

        for rel_name, rel_info in fk_relationships.items():
            from_table = rel_info['from_table']
            from_col = rel_info['from_col']
            to_table = rel_info['to_table']
            to_col = rel_info['to_col']

            # Check for orphaned records (records in from_table that reference non-existent records in to_table)
            query = f"""
                SELECT COUNT(*) FROM {from_table}
                WHERE {from_col} IS NOT NULL
                AND {from_col} NOT IN (SELECT {to_col} FROM {to_table})
            """

            try:
                cursor.execute(query)
                orphaned_count = cursor.fetchone()[0]

                # Also get some examples of orphaned records
                examples_query = f"""
                    SELECT {from_table}.id, {from_table}.{from_col}
                    FROM {from_table}
                    WHERE {from_table}.{from_col} IS NOT NULL
                    AND {from_table}.{from_col} NOT IN (SELECT {to_table}.{to_col} FROM {to_table})
                    LIMIT 5
                """

                cursor.execute(examples_query)
                examples = cursor.fetchall()

                validation_report['relationships'][rel_name] = {
                    'from_table': from_table,
                    'from_col': from_col,
                    'to_table': to_table,
                    'to_col': to_col,
                    'orphaned_records': orphaned_count,
                    'examples': [
                        {'id': row[0], 'invalid_ref': row[1]}
                        for row in examples
                    ]
                }

                if orphaned_count > 0:
                    total_violations += orphaned_count
                    logger.warning(f"Found {orphaned_count} orphaned records in {from_table}.{from_col} -> {to_table}.{to_col}")

            except Exception as e:
                logger.error(f"Error checking FK relationship {rel_name}: {str(e)}")
                validation_report['relationships'][rel_name] = {
                    'error': str(e)
                }

        validation_report['total_relationships_checked'] = len(fk_relationships)
        validation_report['violations_found'] = total_violations

        # Generate summary and recommendations
        validation_report['summary'] = {
            'relationships_checked': len(fk_relationships),
            'total_violations': total_violations,
            'relationships_with_violations': sum(1 for r in validation_report['relationships'].values()
                                                if isinstance(r, dict) and r.get('orphaned_records', 0) > 0),
            'foreign_keys_enabled': validation_report['foreign_keys_enabled']
        }

        recommendations = []
        if not validation_report['foreign_keys_enabled']:
            recommendations.append("Enable foreign key constraints with PRAGMA foreign_keys = ON")

        if total_violations > 0:
            recommendations.append(f"Fix {total_violations} orphaned foreign key references")

        relationships_with_violations = [
            rel_name for rel_name, rel_info in validation_report['relationships'].items()
            if isinstance(rel_info, dict) and rel_info.get('orphaned_records', 0) > 0
        ]

        if relationships_with_violations:
            recommendations.append(f"Review violations in: {', '.join(relationships_with_violations[:3])}")

        validation_report['recommendations'] = recommendations

        logger.info(f"Foreign key validation complete: {len(fk_relationships)} relationships checked, {total_violations} violations found")

        cursor.close()

    except Exception as e:
        logger.error(f"Error during comprehensive foreign key validation: {str(e)}")
        validation_report['error'] = str(e)
        raise
    finally:
        session.close()

    return validation_report


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

    # Foreign key validation issues
    fk_validation = all_issues.get('foreign_key_validation', {})
    report.append("FOREIGN KEY VALIDATION ISSUES")
    report.append("-" * 50)
    if 'foreign_keys_enabled' in fk_validation:
        report.append(f"Foreign keys enabled: {fk_validation['foreign_keys_enabled']}")
    if 'summary' in fk_validation:
        summary = fk_validation['summary']
        report.append(f"Relationships checked: {summary.get('relationships_checked', 0)}")
        report.append(f"Total violations: {summary.get('total_violations', 0)}")
        report.append(f"Relationships with violations: {summary.get('relationships_with_violations', 0)}")
        report.append("")
    if 'relationships' in fk_validation:
        for rel_name, rel_info in fk_validation['relationships'].items():
            if isinstance(rel_info, dict) and 'orphaned_records' in rel_info:
                if rel_info['orphaned_records'] > 0:
                    report.append(f"❌ {rel_name}: {rel_info['orphaned_records']} orphaned records")
                else:
                    report.append(f"✅ {rel_name}: OK")
    if 'recommendations' in fk_validation:
        for rec in fk_validation['recommendations']:
            report.append(f"  - {rec}")
    report.append("")
    
    # Summary
    total_issues = sum(
        len(issues) for category in all_issues.values() 
        for issues in category.values() if isinstance(issues, list)
    )
    report.append(f"SUMMARY: {total_issues} total issues found across all categories")
    
    return "\n".join(report)


def create_cleanup_issues_table() -> None:
    """
    Create the cleanup_issues tracking table and indexes if they do not exist.
    """
    logger.info("Ensuring cleanup_issues table exists...")
    session = StarbaseSession()
    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS cleanup_issues (
                id INTEGER PRIMARY KEY,
                issue_type TEXT NOT NULL,
                category TEXT NOT NULL,
                table_name TEXT,
                record_id INTEGER,
                details TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'OPEN',
                source TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        # Unique index to avoid duplicate rows for the same issue
        cursor.execute(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS ux_cleanup_issues_identity
            ON cleanup_issues (category, issue_type, table_name, record_id, details)
            """
        )
        # Helpful indexes
        cursor.execute(
            """
            CREATE INDEX IF NOT EXISTS ix_cleanup_issues_status
            ON cleanup_issues (status)
            """
        )
        cursor.execute(
            """
            CREATE INDEX IF NOT EXISTS ix_cleanup_issues_table_record
            ON cleanup_issues (table_name, record_id)
            """
        )
        cursor.close()
        session.commit()
        logger.info("cleanup_issues table ready")
    except Exception as e:
        logger.error(f"Error creating cleanup_issues table: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def _infer_table_and_record_id(issue_item: Dict) -> (str, int):
    """
    Infer table_name and record_id from an issue dict by checking common *_id keys.
    Preference order: joined_id → ship_id → accession_id → genome_id → captain_id → taxonomy_id.
    """
    key_order = [
        ("joined_id", "joined_ships"),
        ("ship_id", "ships"),
        ("accession_id", "accessions"),
        ("genome_id", "genomes"),
        ("captain_id", "captains"),
        ("taxonomy_id", "taxonomy"),
    ]
    for key, table in key_order:
        if key in issue_item and issue_item[key] is not None:
            try:
                return table, int(issue_item[key])
            except Exception:
                return table, None
    # Fallbacks for alternate naming
    if "id" in issue_item and isinstance(issue_item["id"], int):
        return None, issue_item["id"]
    return None, None


def record_cleanup_issues(all_issues: Dict, source: str = "pipeline", dry_run: bool = True) -> Dict:
    """
    Persist issues discovered by the cleanup pipeline into cleanup_issues table.

    Args:
        all_issues: Dict produced by the pipeline aggregating issue lists per category.
        source: String label for where these issues originated from.
        dry_run: If True, do not write to the database, just return a summary.

    Returns:
        Dict summary with counts by category and total inserted/skipped.
    """
    logger.info("Recording cleanup issues to cleanup_issues table...")
    create_cleanup_issues_table()

    session = StarbaseSession()
    summary = {"by_category": {}, "inserted": 0, "skipped": 0, "updated": 0}

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        insert_sql = (
            """
            INSERT OR IGNORE INTO cleanup_issues
            (issue_type, category, table_name, record_id, details, status, source)
            VALUES (?, ?, ?, ?, ?, 'OPEN', ?)
            """
        )

        # Record issues found during analysis
        for category_name, category_payload in (all_issues or {}).items():
            # category_payload is expected to be a dict with lists
            inserted_this_category = 0
            if not isinstance(category_payload, dict):
                continue

            for key, value in category_payload.items():
                if not isinstance(value, list):
                    continue
                issue_type = key
                for item in value:
                    if not isinstance(item, dict):
                        # Wrap non-dict entries
                        item = {"value": item}

                    table_name, record_id = _infer_table_and_record_id(item)
                    details_text = json.dumps(item, ensure_ascii=False, sort_keys=True)

                    if dry_run:
                        summary["inserted"] += 1
                        inserted_this_category += 1
                        continue

                    try:
                        cursor.execute(
                            insert_sql,
                            (issue_type, category_name, table_name, record_id, details_text, source),
                        )
                        # Use rowcount to infer if inserted or ignored as duplicate
                        if cursor.rowcount and cursor.rowcount > 0:
                            summary["inserted"] += 1
                            inserted_this_category += 1
                        else:
                            summary["skipped"] += 1
                    except Exception as e:
                        logger.error(f"Failed to insert cleanup issue ({category_name}/{issue_type}): {str(e)}")
                        summary["skipped"] += 1

            if inserted_this_category:
                summary["by_category"][category_name] = inserted_this_category

        if not dry_run:
            cursor.close()
            session.commit()
        else:
            cursor.close()

        logger.info(
            f"Recorded issues summary: inserted={summary['inserted']} skipped={summary['skipped']}"
        )
        return summary

    except Exception as e:
        logger.error(f"Error recording cleanup issues: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def record_cleanup_fixes(all_fixes: Dict, source: str = "pipeline", dry_run: bool = True) -> Dict:
    """
    Record the results of fixes applied by the cleanup pipeline.

    This function records successful fixes and can update the status of related issues
    from 'OPEN' to 'FIXED' or 'RESOLVED'.

    Args:
        all_fixes: Dict containing results from fix operations
        source: String label for where these fixes originated from
        dry_run: If True, do not write to the database, just return a summary

    Returns:
        Dict summary with counts of recorded fixes and status updates
    """
    logger.info("Recording cleanup fixes to cleanup_issues table...")
    create_cleanup_issues_table()

    session = StarbaseSession()
    summary = {"fixes_recorded": 0, "issues_updated": 0, "skipped": 0}

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Record fix results
        fix_insert_sql = (
            """
            INSERT OR IGNORE INTO cleanup_issues
            (issue_type, category, table_name, record_id, details, status, source)
            VALUES (?, ?, ?, ?, ?, 'FIXED', ?)
            """
        )

        # Update existing issue status
        update_status_sql = (
            """
            UPDATE cleanup_issues
            SET status = 'FIXED', updated_at = CURRENT_TIMESTAMP, details = ?
            WHERE category = ? AND issue_type = ? AND table_name = ? AND record_id = ?
            """
        )

        # Process different types of fixes
        for fix_category, fix_data in (all_fixes or {}).items():
            if not isinstance(fix_data, dict) or 'summary' not in fix_data:
                continue

            # Record individual fix actions
            if 'ships_fixed' in fix_data and isinstance(fix_data['ships_fixed'], list):
                for fix in fix_data['ships_fixed']:
                    if not isinstance(fix, dict):
                        continue

                    table_name, record_id = _infer_table_and_record_id(fix)
                    details_text = json.dumps({
                        **fix,
                        'fix_category': fix_category,
                        'fix_timestamp': datetime.now().isoformat()
                    }, ensure_ascii=False, sort_keys=True)

                    if dry_run:
                        summary["fixes_recorded"] += 1
                        continue

                    try:
                        cursor.execute(
                            fix_insert_sql,
                            ('fixed_ship_reference', 'ships_accessions_joined_ships',
                             table_name, record_id, details_text, source),
                        )
                        if cursor.rowcount and cursor.rowcount > 0:
                            summary["fixes_recorded"] += 1
                    except Exception as e:
                        logger.error(f"Failed to record ship fix: {str(e)}")
                        summary["skipped"] += 1

            # Record ship_id relationship fixes
            if 'missing_ship_ids_fixed' in fix_data and isinstance(fix_data['missing_ship_ids_fixed'], list):
                for fix in fix_data['missing_ship_ids_fixed']:
                    if not isinstance(fix, dict):
                        continue

                    table_name, record_id = _infer_table_and_record_id(fix)
                    details_text = json.dumps({
                        **fix,
                        'fix_category': 'ship_id_relationships',
                        'fix_timestamp': datetime.now().isoformat()
                    }, ensure_ascii=False, sort_keys=True)

                    if dry_run:
                        summary["fixes_recorded"] += 1
                        continue

                    try:
                        cursor.execute(
                            fix_insert_sql,
                            ('fixed_ship_id', 'ship_id_relationships',
                             table_name, record_id, details_text, source),
                        )
                        if cursor.rowcount and cursor.rowcount > 0:
                            summary["fixes_recorded"] += 1
                    except Exception as e:
                        logger.error(f"Failed to record ship_id fix: {str(e)}")
                        summary["skipped"] += 1

            # Record joined_ships ship_id fixes
            if 'ship_ids_fixed' in fix_data and isinstance(fix_data['ship_ids_fixed'], list):
                for fix in fix_data['ship_ids_fixed']:
                    if not isinstance(fix, dict):
                        continue

                    table_name, record_id = _infer_table_and_record_id(fix)
                    details_text = json.dumps({
                        **fix,
                        'fix_category': 'joined_ships_ship_id',
                        'fix_timestamp': datetime.now().isoformat()
                    }, ensure_ascii=False, sort_keys=True)

                    if dry_run:
                        summary["fixes_recorded"] += 1
                        continue

                    try:
                        cursor.execute(
                            fix_insert_sql,
                            ('fixed_joined_ship_id', 'joined_ships_ship_id',
                             table_name, record_id, details_text, source),
                        )
                        if cursor.rowcount and cursor.rowcount > 0:
                            summary["fixes_recorded"] += 1
                    except Exception as e:
                        logger.error(f"Failed to record joined_ships fix: {str(e)}")
                        summary["skipped"] += 1

        if not dry_run:
            cursor.close()
            session.commit()
        else:
            cursor.close()

        logger.info(
            f"Recorded fixes summary: fixes={summary['fixes_recorded']} issues_updated={summary['issues_updated']} skipped={summary['skipped']}"
        )
        return summary

    except Exception as e:
        logger.error(f"Error recording cleanup fixes: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def update_issue_status(table_name: str, record_id: int, status: str = "FIXED", details: Dict = None, dry_run: bool = True) -> bool:
    """
    Update the status of a specific issue in the cleanup_issues table.

    Args:
        table_name: Name of the table the issue relates to
        record_id: ID of the record the issue relates to
        status: New status ('OPEN', 'FIXED', 'RESOLVED', 'CLOSED')
        details: Additional details to add to the existing details
        dry_run: If True, do not write to the database

    Returns:
        bool: True if update was successful
    """
    if dry_run:
        logger.info(f"DRY RUN: Would update issue status for {table_name}.{record_id} to {status}")
        return True

    session = StarbaseSession()
    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Get current details
        cursor.execute(
            "SELECT details FROM cleanup_issues WHERE table_name = ? AND record_id = ? ORDER BY id DESC LIMIT 1",
            (table_name, record_id)
        )
        current_details_row = cursor.fetchone()

        if current_details_row:
            current_details = json.loads(current_details_row[0])
            if details:
                current_details.update(details)
            updated_details = json.dumps(current_details, ensure_ascii=False, sort_keys=True)

            cursor.execute(
                "UPDATE cleanup_issues SET status = ?, details = ?, updated_at = CURRENT_TIMESTAMP WHERE table_name = ? AND record_id = ?",
                (status, updated_details, table_name, record_id)
            )

            session.commit()
            logger.info(f"Updated issue status for {table_name}.{record_id} to {status}")
            return True
        else:
            logger.warning(f"No issue found for {table_name}.{record_id}")
            return False

    except Exception as e:
        logger.error(f"Error updating issue status: {str(e)}")
        session.rollback()
        return False
    finally:
        session.close()


def get_cleanup_issues_summary() -> Dict:
    """
    Get a summary of issues in the cleanup_issues table.

    Returns:
        Dict: Summary statistics about cleanup issues
    """
    logger.info("Getting cleanup issues summary...")
    session = StarbaseSession()

    summary = {
        'total_issues': 0,
        'open_issues': 0,
        'fixed_issues': 0,
        'issues_by_category': {},
        'status_breakdown': {},
        'recent_issues': []
    }

    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()

        # Get total counts
        cursor.execute("SELECT COUNT(*) FROM cleanup_issues")
        summary['total_issues'] = cursor.fetchone()[0]

        # Get status breakdown
        cursor.execute("SELECT status, COUNT(*) FROM cleanup_issues GROUP BY status")
        summary['status_breakdown'] = dict(cursor.fetchall())

        summary['open_issues'] = summary['status_breakdown'].get('OPEN', 0)
        summary['fixed_issues'] = summary['status_breakdown'].get('FIXED', 0)

        # Get issues by category
        cursor.execute("SELECT category, COUNT(*) FROM cleanup_issues GROUP BY category")
        summary['issues_by_category'] = dict(cursor.fetchall())

        # Get recent issues (last 10)
        cursor.execute("""
            SELECT category, issue_type, table_name, record_id, status, created_at
            FROM cleanup_issues
            ORDER BY created_at DESC
            LIMIT 10
        """)
        recent_rows = cursor.fetchall()
        summary['recent_issues'] = [
            {
                'category': row[0],
                'issue_type': row[1],
                'table_name': row[2],
                'record_id': row[3],
                'status': row[4],
                'created_at': row[5][:19] if row[5] else 'N/A'  # Format timestamp
            }
            for row in recent_rows
        ]

        cursor.close()

        logger.info(f"Found {summary['total_issues']} total issues, {summary['open_issues']} open, {summary['fixed_issues']} fixed")

    except Exception as e:
        logger.error(f"Error getting cleanup issues summary: {str(e)}")
        raise
    finally:
        session.close()

    return summary


def fix_missing_genome_info(dry_run: bool = True) -> Dict:
    """
    Fix missing genome information, particularly missing taxonomy_id.
    
    Strategy:
    1. For genomes with missing taxonomy_id but existing ome, look for other genomes
       with the same ome that have a taxonomy_id and copy it
    2. For genomes with missing ome but existing taxonomy_id, look for other genomes
       with the same taxonomy_id that have an ome and copy it
    3. Optionally, could integrate with NCBI taxonomy database for more comprehensive fixes
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of fixes applied or would be applied
    """
    logger.info("Fixing missing genome information...")
    session = StarbaseSession()
    
    fixes = {
        'taxonomy_id_fixed': [],
        'ome_fixed': [],
        'no_match_found': [],
        'summary': {}
    }
    
    try:
        # Find genomes with missing taxonomy_id but existing ome
        genomes_missing_taxonomy = session.query(Genome).filter(
            (Genome.taxonomy_id.is_(None)) & 
            (Genome.ome.isnot(None)) & 
            (Genome.ome != '')
        ).all()
        
        logger.info(f"Found {len(genomes_missing_taxonomy)} genomes with missing taxonomy_id but existing ome")
        
        for genome in genomes_missing_taxonomy:
            # Look for other genomes with the same ome that have a taxonomy_id
            matching_genome = session.query(Genome).filter(
                (Genome.ome == genome.ome) & 
                (Genome.taxonomy_id.isnot(None)) &
                (Genome.id != genome.id)
            ).first()
            
            if matching_genome:
                if not dry_run:
                    old_taxonomy_id = genome.taxonomy_id
                    genome.taxonomy_id = matching_genome.taxonomy_id
                    session.add(genome)
                    
                    fixes['taxonomy_id_fixed'].append({
                        'genome_id': genome.id,
                        'ome': genome.ome,
                        'old_taxonomy_id': old_taxonomy_id,
                        'new_taxonomy_id': matching_genome.taxonomy_id,
                        'source_genome_id': matching_genome.id,
                        'action': f'Copied taxonomy_id {matching_genome.taxonomy_id} from genome {matching_genome.id}'
                    })
                else:
                    fixes['taxonomy_id_fixed'].append({
                        'genome_id': genome.id,
                        'ome': genome.ome,
                        'old_taxonomy_id': None,
                        'new_taxonomy_id': matching_genome.taxonomy_id,
                        'source_genome_id': matching_genome.id,
                        'action': f'Would copy taxonomy_id {matching_genome.taxonomy_id} from genome {matching_genome.id}'
                    })
            else:
                fixes['no_match_found'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'issue': 'No other genome with same ome found to copy taxonomy_id from'
                })
        
        # Find genomes with missing ome but existing taxonomy_id
        genomes_missing_ome = session.query(Genome).filter(
            (Genome.ome.is_(None) | (Genome.ome == '')) & 
            (Genome.taxonomy_id.isnot(None))
        ).all()
        
        logger.info(f"Found {len(genomes_missing_ome)} genomes with missing ome but existing taxonomy_id")
        
        for genome in genomes_missing_ome:
            # Look for other genomes with the same taxonomy_id that have an ome
            matching_genome = session.query(Genome).filter(
                (Genome.taxonomy_id == genome.taxonomy_id) & 
                (Genome.ome.isnot(None)) & 
                (Genome.ome != '') &
                (Genome.id != genome.id)
            ).first()
            
            if matching_genome:
                if not dry_run:
                    old_ome = genome.ome
                    genome.ome = matching_genome.ome
                    session.add(genome)
                    
                    fixes['ome_fixed'].append({
                        'genome_id': genome.id,
                        'taxonomy_id': genome.taxonomy_id,
                        'old_ome': old_ome,
                        'new_ome': matching_genome.ome,
                        'source_genome_id': matching_genome.id,
                        'action': f'Copied ome {matching_genome.ome} from genome {matching_genome.id}'
                    })
                else:
                    fixes['ome_fixed'].append({
                        'genome_id': genome.id,
                        'taxonomy_id': genome.taxonomy_id,
                        'old_ome': None,
                        'new_ome': matching_genome.ome,
                        'source_genome_id': matching_genome.id,
                        'action': f'Would copy ome {matching_genome.ome} from genome {matching_genome.id}'
                    })
            else:
                fixes['no_match_found'].append({
                    'genome_id': genome.id,
                    'taxonomy_id': genome.taxonomy_id,
                    'issue': 'No other genome with same taxonomy_id found to copy ome from'
                })
        
        if not dry_run:
            session.commit()
            logger.info("Applied genome information fixes")
        else:
            logger.info("Dry run completed - no changes made")
        
        # Summary
        fixes['summary'] = {
            'taxonomy_id_fixed_count': len(fixes['taxonomy_id_fixed']),
            'ome_fixed_count': len(fixes['ome_fixed']),
            'no_match_found_count': len(fixes['no_match_found']),
            'recommendation': 'Consider NCBI taxonomy lookup for remaining unfixed genomes'
        }
        
        logger.info(f"Fixed {len(fixes['taxonomy_id_fixed'])} missing taxonomy_id issues")
        logger.info(f"Fixed {len(fixes['ome_fixed'])} missing ome issues")
        logger.info(f"Could not fix {len(fixes['no_match_found'])} genomes (no matching reference found)")
        
    except Exception as e:
        logger.error(f"Error fixing missing genome information: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes


def analyze_missing_genome_info() -> Dict:
    """
    Analyze the scope of missing genome information to help prioritize fixes.
    
    Returns:
        Dict: Analysis of missing genome information
    """
    logger.info("Analyzing missing genome information...")
    session = StarbaseSession()
    
    analysis = {
        'missing_taxonomy_id': [],
        'missing_ome': [],
        'missing_both': [],
        'potential_fixes': [],
        'summary': {}
    }
    
    try:
        # Find all genomes with missing information
        all_genomes = session.query(Genome).all()
        
        for genome in all_genomes:
            missing_taxonomy = genome.taxonomy_id is None
            missing_ome = genome.ome is None or genome.ome == ''
            
            if missing_taxonomy and missing_ome:
                analysis['missing_both'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'taxonomy_id': genome.taxonomy_id,
                    'issue': 'Missing both ome and taxonomy_id'
                })
            elif missing_taxonomy:
                analysis['missing_taxonomy_id'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'taxonomy_id': genome.taxonomy_id,
                    'issue': 'Missing taxonomy_id'
                })
            elif missing_ome:
                analysis['missing_ome'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'taxonomy_id': genome.taxonomy_id,
                    'issue': 'Missing ome'
                })
        
        # Analyze potential fixes
        # Group by ome to see how many could be fixed by copying taxonomy_id
        ome_groups = {}
        for genome in analysis['missing_taxonomy_id']:
            ome = genome['ome']
            if ome not in ome_groups:
                ome_groups[ome] = []
            ome_groups[ome].append(genome)
        
        # Find omes that have some genomes with taxonomy_id
        for ome, genomes in ome_groups.items():
            # Check if any genome with this ome has a taxonomy_id
            reference_genome = session.query(Genome).filter(
                (Genome.ome == ome) & 
                (Genome.taxonomy_id.isnot(None))
            ).first()
            
            if reference_genome:
                analysis['potential_fixes'].append({
                    'fix_type': 'copy_taxonomy_id',
                    'ome': ome,
                    'genomes_to_fix': len(genomes),
                    'reference_genome_id': reference_genome.id,
                    'reference_taxonomy_id': reference_genome.taxonomy_id,
                    'description': f'Copy taxonomy_id {reference_genome.taxonomy_id} to {len(genomes)} genomes with ome {ome}'
                })
        
        # Group by taxonomy_id to see how many could be fixed by copying ome
        taxonomy_groups = {}
        for genome in analysis['missing_ome']:
            taxonomy_id = genome['taxonomy_id']
            if taxonomy_id not in taxonomy_groups:
                taxonomy_groups[taxonomy_id] = []
            taxonomy_groups[taxonomy_id].append(genome)
        
        for taxonomy_id, genomes in taxonomy_groups.items():
            # Check if any genome with this taxonomy_id has an ome
            reference_genome = session.query(Genome).filter(
                (Genome.taxonomy_id == taxonomy_id) & 
                (Genome.ome.isnot(None)) & 
                (Genome.ome != '')
            ).first()
            
            if reference_genome:
                analysis['potential_fixes'].append({
                    'fix_type': 'copy_ome',
                    'taxonomy_id': taxonomy_id,
                    'genomes_to_fix': len(genomes),
                    'reference_genome_id': reference_genome.id,
                    'reference_ome': reference_genome.ome,
                    'description': f'Copy ome {reference_genome.ome} to {len(genomes)} genomes with taxonomy_id {taxonomy_id}'
                })
        
        # Summary
        analysis['summary'] = {
            'total_genomes': len(all_genomes),
            'missing_taxonomy_id_count': len(analysis['missing_taxonomy_id']),
            'missing_ome_count': len(analysis['missing_ome']),
            'missing_both_count': len(analysis['missing_both']),
            'potential_fixes_count': len(analysis['potential_fixes']),
            'recommendation': 'Run fix_missing_genome_info_via_ome_map() to apply available fixes'
        }
        
        logger.info(f"Analysis complete:")
        logger.info(f"  - {len(analysis['missing_taxonomy_id'])} genomes missing taxonomy_id")
        logger.info(f"  - {len(analysis['missing_ome'])} genomes missing ome")
        logger.info(f"  - {len(analysis['missing_both'])} genomes missing both")
        logger.info(f"  - {len(analysis['potential_fixes'])} potential fixes identified")
        
    except Exception as e:
        logger.error(f"Error analyzing missing genome information: {str(e)}")
        raise
    finally:
        session.close()
    
    return analysis


def fix_missing_genome_taxonomy_from_joined_ships(dry_run: bool = True) -> Dict:
    """
    Fill missing Genome.taxonomy_id using the most common joined_ships.tax_id
    for entries linked to that genome.

    Args:
        dry_run: If True, do not write changes, only report.

    Returns:
        Dict: Report of updates applied or that would be applied.
    """
    logger.info("Setting missing Genome.taxonomy_id from joined_ships.tax_id...")
    session = StarbaseSession()

    report = {
        'genomes_updated': [],
        'no_tax_id_found': [],
        'summary': {}
    }

    try:
        # Get genomes missing taxonomy_id
        genomes = session.query(Genome).filter(Genome.taxonomy_id.is_(None)).all()
        logger.info(f"Found {len(genomes)} genomes missing taxonomy_id")

        total_updated = 0
        for genome in genomes:
            # Aggregate tax_id counts from joined_ships for this genome
            tax_counts = (
                session.query(JoinedShips.tax_id, func.count(JoinedShips.id))
                .filter(
                    (JoinedShips.genome_id == genome.id) &
                    (JoinedShips.tax_id.isnot(None))
                )
                .group_by(JoinedShips.tax_id)
                .order_by(func.count(JoinedShips.id).desc())
                .all()
            )

            if tax_counts:
                chosen_tax_id = tax_counts[0][0]
                if not dry_run:
                    old_taxonomy_id = genome.taxonomy_id
                    genome.taxonomy_id = chosen_tax_id
                    session.add(genome)
                report['genomes_updated'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'new_taxonomy_id': chosen_tax_id,
                    'source': 'joined_ships.mode_tax_id'
                })
                total_updated += 1
            else:
                report['no_tax_id_found'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'issue': 'No joined_ships.tax_id for this genome'
                })

        if not dry_run:
            session.commit()
            logger.info(f"Updated taxonomy_id for {total_updated} genomes")
        else:
            logger.info(f"Would update taxonomy_id for {total_updated} genomes")

        report['summary'] = {
            'total_missing': len(genomes),
            'total_updated': total_updated,
            'total_unresolved': len(report['no_tax_id_found'])
        }

    except Exception as e:
        logger.error(f"Error setting genome taxonomy from joined_ships: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()

    return report


def _ncbi_taxonomy_search(scientific_name: str, email: str = None, api_key: str = None) -> Dict:
    """
    Query NCBI Taxonomy for a scientific name and return minimal info.
    Returns dict with keys: taxid (str), scientific_name (str). Empty dict if not found.
    """
    params = {
        'db': 'taxonomy',
        'term': scientific_name,
        'retmode': 'json'
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        ids = data.get('esearchresult', {}).get('idlist', [])
        if not ids:
            return {}
        taxid = ids[0]
        # Fetch summary to get scientific name
        sum_params = {
            'db': 'taxonomy',
            'id': taxid,
            'retmode': 'json'
        }
        if email:
            sum_params['email'] = email
        if api_key:
            sum_params['api_key'] = api_key
        sum_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' + urllib.parse.urlencode(sum_params)
        with urllib.request.urlopen(sum_url, timeout=15) as resp2:
            sum_data = json.loads(resp2.read().decode('utf-8'))
        docsum = (sum_data.get('result', {}) or {}).get(taxid, {})
        sci_name = docsum.get('scientificname') or docsum.get('title') or scientific_name
        return {'taxid': str(taxid), 'scientific_name': sci_name}
    except Exception as e:
        logger.warning(f"NCBI lookup failed for '{scientific_name}': {str(e)}")
        return {}


def _find_or_create_taxonomy(name: str, taxid: str = None, session=None) -> int:
    """
    Find existing taxonomy entry or create new one without duplicating.
    
    Args:
        name (str): Scientific name to search for
        taxid (str): NCBI taxID if available
        session: Database session
        
    Returns:
        int: Taxonomy ID (existing or newly created)
    """
    if not session:
        session = StarbaseSession()
    
    # First try to find by name (exact match)
    existing = session.query(Taxonomy).filter(Taxonomy.name == name).first()
    if existing:
        return existing.id
    
    # If taxID provided, check by taxID
    if taxid:
        existing_by_taxid = session.query(Taxonomy).filter(Taxonomy.taxID == taxid).first()
        if existing_by_taxid:
            return existing_by_taxid.id
    
    # Create new entry
    new_taxonomy = Taxonomy(name=name, taxID=taxid)
    session.add(new_taxonomy)
    session.flush()  # Get the ID without committing
    return new_taxonomy.id


def fix_missing_genome_taxonomy_via_ncbi(map_csv: str = None, email: str = None, api_key: str = None, dry_run: bool = True, delay_seconds: float = 0.34) -> Dict:
    """
    Resolve missing Genome.taxonomy_id by querying NCBI Taxonomy.
    
    Workflow:
    - Optionally load an ome->scientific_name CSV map to derive names.
    - For each genome with missing taxonomy_id, get candidate name from map or skip.
    - Query NCBI Taxonomy (esearch+esummary) to get taxid and scientific name.
    - If a taxonomy row with matching taxID or name exists, reuse it.
    - Else, create a new taxonomy row (minimal fields) and link genome to it.
    - Respects dry_run; delay between requests to be polite (default ~3/sec).
    
    Args:
        map_csv: Path to CSV with columns ome,scientific_name (header required).
        email: Contact email for NCBI E-utilities etiquette.
        api_key: Optional NCBI API key.
        dry_run: If True, don't write changes.
        delay_seconds: Sleep between NCBI calls to respect rate limits.
    
    Returns:
        Dict summary and lists of actions.
    """
    logger.info("Fixing missing genome taxonomy via NCBI...")
    session = StarbaseSession()

    report = {
        'genomes_linked': [],
        'taxa_created': [],
        'skipped_no_name': [],
        'skipped_not_found': [],
        'errors': [],
        'summary': {}
    }

    ome_to_name = {}
    if map_csv:
        try:
            with open(map_csv, 'r', newline='') as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    ome = (row.get('ome') or '').strip()
                    name = (row.get('scientific_name') or '').strip()
                    if ome and name:
                        ome_to_name[ome] = name
            logger.info(f"Loaded {len(ome_to_name)} ome->name mappings from {map_csv}")
        except Exception as e:
            logger.error(f"Failed to read map CSV {map_csv}: {str(e)}")
            raise

    try:
        genomes = session.query(Genome).filter(Genome.taxonomy_id.is_(None)).all()
        logger.info(f"Processing {len(genomes)} genomes with missing taxonomy_id for NCBI lookup")

        linked = 0
        created_taxa = 0
        for genome in genomes:
            candidate_name = ome_to_name.get(genome.ome) if genome.ome else None
            if not candidate_name:
                report['skipped_no_name'].append({'genome_id': genome.id, 'ome': genome.ome, 'issue': 'No mapping to scientific name'})
                continue

            ncbi = _ncbi_taxonomy_search(candidate_name, email=email, api_key=api_key)
            if not ncbi:
                report['skipped_not_found'].append({'genome_id': genome.id, 'ome': genome.ome, 'query': candidate_name})
                time.sleep(delay_seconds)
                continue

            taxid_str = ncbi['taxid']
            sci_name = ncbi['scientific_name']

            # Try to find existing taxonomy row by taxID or name
            existing = session.query(Taxonomy).filter(
                (Taxonomy.taxID == taxid_str) | (Taxonomy.name == sci_name)
            ).first()

            taxonomy_row = existing
            if not taxonomy_row and not dry_run:
                # Use helper to prevent duplicates
                taxonomy_id = _find_or_create_taxonomy(sci_name, taxid_str, session)
                taxonomy_row = session.query(Taxonomy).filter(Taxonomy.id == taxonomy_id).first()
                
                # Only count as created if it's new
                if taxonomy_id not in [t['taxonomy_id'] for t in report['taxa_created']]:
                    report['taxa_created'].append({'taxonomy_id': taxonomy_id, 'name': sci_name, 'taxID': taxid_str})
                    created_taxa += 1

            if taxonomy_row:
                if not dry_run:
                    genome.taxonomy_id = taxonomy_row.id
                    session.add(genome)
                report['genomes_linked'].append({'genome_id': genome.id, 'ome': genome.ome, 'taxonomy_id': taxonomy_row.id, 'name': sci_name})
                linked += 1
            else:
                # dry_run path: we still need a placeholder id in report
                report['genomes_linked'].append({'genome_id': genome.id, 'ome': genome.ome, 'taxonomy_id': '(new)', 'name': sci_name})
                linked += 1

            time.sleep(delay_seconds)

        if not dry_run:
            session.commit()
            logger.info(f"Linked {linked} genomes; created {created_taxa} taxonomy rows")
        else:
            logger.info(f"Would link {linked} genomes; would create {created_taxa} taxonomy rows")

        report['summary'] = {
            'total_missing': len(genomes),
            'linked': linked,
            'taxa_created': created_taxa,
            'skipped_no_name': len(report['skipped_no_name']),
            'skipped_not_found': len(report['skipped_not_found'])
        }
        return report

    except Exception as e:
        logger.error(f"Error during NCBI taxonomy fix: {str(e)}")
        if not dry_run:
            session.rollback()
        report['errors'].append(str(e))
        raise
    finally:
        session.close()


def _parse_taxdump_names(names_dmp_path: str, include_synonyms: bool = False) -> Dict[str, str]:
    """
    Parse NCBI taxdump names.dmp and return a mapping name -> taxID (as string).
    By default only 'scientific name' entries are included; if include_synonyms is True,
    also include 'synonym' and 'equivalent name' classes.
    """
    wanted = {"scientific name"}
    if include_synonyms:
        wanted.update({"synonym", "equivalent name", "genbank common name", "common name"})

    mapping = {}
    try:
        with open(names_dmp_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                # Format: tax_id | name_txt | unique name | name class |
                parts = [p.strip() for p in line.split('|')]
                if len(parts) < 4:
                    continue
                tax_id = parts[0]
                name_txt = parts[1]
                name_class = parts[3]
                if name_class in wanted and name_txt:
                    # Prefer first seen scientific name mapping; don't overwrite
                    if name_txt not in mapping:
                        mapping[name_txt] = tax_id
        logger.info(f"Loaded {len(mapping)} names from {names_dmp_path} (include_synonyms={include_synonyms})")
        return mapping
    except Exception as e:
        logger.error(f"Failed to parse names.dmp at {names_dmp_path}: {str(e)}")
        raise


def _normalize_name(value: str) -> str:
    if value is None:
        return ''
    s = value.replace('_', ' ').strip().lower()
    # collapse multiple spaces
    return ' '.join(s.split())


def fix_missing_genome_taxonomy_via_taxdump(names_dmp_path: str, include_synonyms: bool = False, dry_run: bool = True) -> Dict:
    """
    Resolve missing Genome.taxonomy_id by matching genome.ome to names in NCBI taxdump names.dmp.
    Matching strategy:
      - names.dmp exact/case-insensitive/underscore-to-space match
      - fallback to existing Taxonomy rows by (in order):
        1) Taxonomy.name
        2) Taxonomy.genus + ' ' + Taxonomy.species
        3) Taxonomy.genus
    If a matching Taxonomy already exists (by taxID or name), reuse it; otherwise create it (minimal fields).
    """
    logger.info("Fixing missing genome taxonomy via taxdump names.dmp...")
    session = StarbaseSession()

    report = {
        'genomes_linked': [],
        'taxa_created': [],
        'skipped_no_match': [],
        'summary': {}
    }

    name_to_taxid = _parse_taxdump_names(names_dmp_path, include_synonyms=include_synonyms)
    # Build normalized lookup maps for robustness
    lower_map = {k.lower(): v for k, v in name_to_taxid.items()}
    underspaced_map = {k.replace('_', ' ').lower(): v for k, v in name_to_taxid.items()}

    try:
        genomes = session.query(Genome).filter(Genome.taxonomy_id.is_(None)).all()
        logger.info(f"Processing {len(genomes)} genomes with missing taxonomy_id using names.dmp + taxonomy fallbacks")

        linked = 0
        created = 0
        for genome in genomes:
            ome = genome.ome or ''
            if ome == '':
                report['skipped_no_match'].append({'genome_id': genome.id, 'ome': genome.ome})
                continue

            norm_ome = _normalize_name(ome)
            taxid = None
            sci_name = None
            taxonomy_row = None

            # 1) Try names.dmp mappings
            if ome in name_to_taxid:
                taxid = name_to_taxid[ome]
                sci_name = ome
            else:
                low = ome.lower()
                if low in lower_map:
                    taxid = lower_map[low]
                    sci_name = next((k for k, v in name_to_taxid.items() if k.lower() == low and v == taxid), ome)
                else:
                    norm = ome.replace('_', ' ').lower()
                    if norm in underspaced_map:
                        taxid = underspaced_map[norm]
                        sci_name = next((k for k, v in name_to_taxid.items() if k.replace('_',' ').lower() == norm and v == taxid), ome.replace('_', ' '))

            if taxid:
                taxonomy_row = session.query(Taxonomy).filter(
                    (Taxonomy.taxID == str(taxid)) | (Taxonomy.name == sci_name)
                ).first()
                if not taxonomy_row and not dry_run:
                    # Use helper to prevent duplicates
                    taxonomy_id = _find_or_create_taxonomy(sci_name, str(taxid), session)
                    taxonomy_row = session.query(Taxonomy).filter(Taxonomy.id == taxonomy_id).first()
                    
                    # Only count as created if it's new
                    if taxonomy_id not in [t['taxonomy_id'] for t in report['taxa_created']]:
                        report['taxa_created'].append({'taxonomy_id': taxonomy_id, 'name': sci_name, 'taxID': str(taxid)})
                        created += 1

            # 2) Fallback: Existing Taxonomy by normalized name
            if taxonomy_row is None:
                taxonomy_row = (
                    session.query(Taxonomy)
                    .filter(Taxonomy.name.isnot(None))
                    .all()
                )
                taxonomy_row = next((t for t in taxonomy_row if _normalize_name(t.name) == norm_ome), None)

            # 3) Fallback: genus + species
            if taxonomy_row is None:
                candidates = session.query(Taxonomy).filter(
                    (Taxonomy.genus.isnot(None)) & (Taxonomy.genus != '') &
                    (Taxonomy.species.isnot(None)) & (Taxonomy.species != '')
                ).all()
                taxonomy_row = next((t for t in candidates if _normalize_name(f"{t.genus} {t.species}") == norm_ome), None)

            # 4) Fallback: genus only
            if taxonomy_row is None:
                candidates = session.query(Taxonomy).filter(
                    (Taxonomy.genus.isnot(None)) & (Taxonomy.genus != '')
                ).all()
                taxonomy_row = next((t for t in candidates if _normalize_name(t.genus) == norm_ome), None)

            if taxonomy_row:
                if not dry_run:
                    genome.taxonomy_id = taxonomy_row.id
                    session.add(genome)
                report['genomes_linked'].append({'genome_id': genome.id, 'ome': genome.ome, 'taxonomy_id': taxonomy_row.id, 'name': taxonomy_row.name})
                linked += 1
            else:
                report['skipped_no_match'].append({'genome_id': genome.id, 'ome': genome.ome})

        if not dry_run:
            session.commit()
            logger.info(f"Linked {linked} genomes; created {created} taxonomy rows (names.dmp + fallbacks)")
        else:
            logger.info(f"Would link {linked} genomes; would create {created} taxonomy rows (names.dmp + fallbacks)")

        report['summary'] = {
            'total_missing': len(genomes),
            'linked': linked,
            'taxa_created': created,
            'skipped_no_match': len(report['skipped_no_match'])
        }
        return report

    except Exception as e:
        logger.error(f"Error during taxdump-based taxonomy fix: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()


def _parse_ome_map(map_path: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a tab-delimited mapping file with columns like:
    ome\tgenus\tspecies_epithet\tstrain\tjson_tax\t...
    Returns map: ome -> {genus, species, strain, name}
    Where name is 'Genus species' (and optionally strain ignored for taxonomy name).
    """
    mapping: Dict[str, Dict[str, str]] = {}
    try:
        with open(map_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                cols = line.split('\t')
                if len(cols) < 3:
                    continue
                ome = cols[0].strip()
                genus = (cols[1] or '').strip()
                species_part = (cols[2] or '').strip()
                strain = (cols[3] or '').strip() if len(cols) > 3 else ''
                json_blob = (cols[4] or '').strip() if len(cols) > 4 else ''
                assembly_accession = (cols[10] or '').strip() if len(cols) > 10 else ''

                # Prefer JSON values when available
                try:
                    if json_blob:
                        jd = json.loads(json_blob)
                        genus_json = (jd.get('genus') or '').strip()
                        species_json = (jd.get('species') or '').strip()
                        strain_json = (jd.get('strain') or '').strip()
                        if genus_json:
                            genus = genus_json
                        if species_json:
                            # JSON 'species' might include full 'Genus species'; prefer epithet if genus present
                            if genus and species_json.lower().startswith(genus.lower() + ' '):
                                species_part = species_json[len(genus) + 1:].strip()
                            else:
                                species_part = species_json
                        if strain_json:
                            strain = strain_json
                except Exception:
                    pass

                if not ome or not genus:
                    continue
                species_epithet = species_part
                # Build scientific name
                sci_name = genus if not species_epithet else f"{genus} {species_epithet}".strip()
                mapping[ome] = {
                    'genus': genus,
                    'species_epithet': species_epithet,
                    'strain': strain,
                    'name': sci_name,
                    'assembly_accession': assembly_accession
                }
        logger.info(f"Loaded {len(mapping)} OME taxonomy mappings from {map_path}")
        return mapping
    except Exception as e:
        logger.error(f"Failed to parse OME taxonomy map at {map_path}: {str(e)}")
        raise

def reconcile_genome_taxonomy_via_map(map_path: str, dry_run: bool = True) -> Dict:
    """
    Reconcile existing genome entries against OME mapping to ensure consistency.
    - If mapped assembly accession doesn't match genome.assembly_accession, update to the mapped assembly accession
    - If mapped name (Genus species) doesn't match taxonomy.name or (genus,species), update to the mapped taxonomy
        - (reuse existing row by name/genus+species or create a new one).
    Only processes genomes that already have taxonomy_id assigned.
    """
    logger.info("Reconciling genome info according to OME mapping file...")
    session = StarbaseSession()

    report = {
        'genomes_linked': [],
        'genomes_updated': [],
        'taxa_created': [],
        'skipped_no_map': [],
        'skipped_ambiguous_genus': [],
        'assembly_accessions_updated': [],
        'summary': {}
    }

    ome_map = _parse_ome_map(map_path)

    try:
        genomes = session.query(Genome).filter(Genome.taxonomy_id.isnot(None)).all()
        logger.info(f"Checking {len(genomes)} genomes with existing taxonomy_id against OME map")

        linked = 0
        created = 0
        updated = 0
        assembly_accessions_updated = 0

        for genome in genomes:
            ome = genome.ome or ''
            assembly_accession = genome.assembly_accession or ''
            if ome == '' or ome not in ome_map:
                report['skipped_no_map'].append({'genome_id': genome.id, 'ome': genome.ome})
                continue

            mapped_genus = ome_map[ome]['genus']
            mapped_species_epithet = ome_map[ome]['species_epithet']
            mapped_name = ome_map[ome]['name']
            mapped_assembly_accession = ome_map[ome]['assembly_accession']

            current_tax = session.query(Taxonomy).get(genome.taxonomy_id)
            current_name = current_tax.name if current_tax else None
            current_genus = current_tax.genus if current_tax else None
            current_species = current_tax.species if current_tax else None

            matches = (
                (current_name == mapped_name) or
                (current_genus == mapped_genus and (current_species or '') == (mapped_species_epithet or ''))
            )
            if matches:
                continue

            # Find or create the mapped taxonomy row
            taxonomy_row = session.query(Taxonomy).filter(Taxonomy.name == mapped_name).first()
            if taxonomy_row is None and mapped_genus and mapped_species_epithet:
                taxonomy_row = session.query(Taxonomy).filter(
                    (Taxonomy.genus == mapped_genus) & (Taxonomy.species == mapped_species_epithet)
                ).first()

            # If mapped assembly accession doesn't match genome.assembly_accession, update to the mapped assembly accession
            if assembly_accession != mapped_assembly_accession and not dry_run:
                genome.assembly_accession = mapped_assembly_accession
                session.add(genome)
                report['assembly_accessions_updated'].append({'genome_id': genome.id, 'ome': genome.ome, 'old_assembly_accession': assembly_accession, 'new_assembly_accession': mapped_assembly_accession})
                assembly_accessions_updated += 1

            if taxonomy_row is None and not dry_run:
                # Use helper to prevent duplicates
                taxonomy_id = _find_or_create_taxonomy(mapped_name, None, session)
                taxonomy_row = session.query(Taxonomy).filter(Taxonomy.id == taxonomy_id).first()
                
                # Update genus and species if available
                if mapped_genus and not taxonomy_row.genus:
                    taxonomy_row.genus = mapped_genus
                if mapped_species_epithet and not taxonomy_row.species:
                    taxonomy_row.species = mapped_species_epithet
                session.add(taxonomy_row)
                
                # Only count as created if it's new
                if taxonomy_id not in [t['taxonomy_id'] for t in report['taxa_created']]:
                    report['taxa_created'].append({'taxonomy_id': taxonomy_id, 'name': mapped_name})
                    created += 1

            if taxonomy_row:
                if not dry_run:
                    genome.taxonomy_id = taxonomy_row.id
                    session.add(genome)
                report['genomes_updated'].append({
                    'genome_id': genome.id,
                    'ome': genome.ome,
                    'old_taxonomy_id': current_tax.id if current_tax else None,
                    'new_taxonomy_id': taxonomy_row.id,
                    'name': taxonomy_row.name
                })
                updated += 1

        if not dry_run:
            session.commit()
            logger.info(f"Updated taxonomy for {updated} genomes; created {created} taxonomy rows (reconcile)")
        else:
            logger.info(f"Would update taxonomy for {updated} genomes; would create {created} taxonomy rows (reconcile)")

        report['summary'] = {
            'checked': len(genomes),
            'updated': updated,
            'taxa_created': created,
            'skipped_no_map': len(report['skipped_no_map']),
            'assembly_accessions_updated': assembly_accessions_updated
        }
        return report

    except Exception as e:
        logger.error(f"Error during OME map taxonomy reconciliation: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()


def _normalize_contig_id(raw_contig_id: str) -> str:
    """
    Normalize contigID for lookup:
    - Remove leading genome OME code prefix before the first underscore (e.g., altals1_)
    - Also remove a literal 'ome' prefix if present
    - Trim whitespace
    """
    if raw_contig_id is None:
        return ''
    s = raw_contig_id.strip()
    # Remove literal 'ome' prefix if the field was mistakenly prefixed that way
    if s.startswith('ome'):
        s = s[3:]
    s = s.strip()
    # Remove leading OME token up to first underscore (e.g., altals1_)
    if '_' in s:
        # Only strip if it looks like a compact token (letters/digits/dots) followed by underscore
        # Allow stripping even if there's a leading underscore from 'ome' removal
        m = re.match(r'^_?[A-Za-z0-9\.]+_(.+)$', s)
        if m:
            s = m.group(1)
    return s.strip()


def _ncbi_esearch_assembly(term: str, email: str = None, api_key: str = None) -> str:
    """
    Query NCBI E-utilities esearch for an assembly UID given a term.
    Returns assembly UID as string or '' if none.
    """
    params = {
        'db': 'assembly',
        'term': term,
        'retmode': 'json',
        'retmax': '1'
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?' + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        ids = data.get('esearchresult', {}).get('idlist', [])
        return ids[0] if ids else ''
    except Exception as e:
        logger.warning(f"esearch failure for term '{term}': {str(e)}")
        return ''


def _ncbi_assembly_esummary(uid: str, email: str = None, api_key: str = None) -> Dict:
    """
    Query NCBI E-utilities esummary for an assembly UID.
    Returns dict with keys: assemblyaccession, organism, bioproject, biosample.
    """
    if not uid:
        return {}
    params = {
        'db': 'assembly',
        'id': uid,
        'retmode': 'json'
    }
    if email:
        params['email'] = email
    if api_key:
        params['api_key'] = api_key
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?' + urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        result = (data.get('result') or {})
        doc = result.get(uid) or {}
        return {
            'assemblyaccession': doc.get('assemblyaccession') or '',
            'organism': doc.get('organism') or '',
            'bioproject': doc.get('bioproject') or '',
            'biosample': doc.get('biosample') or ''
        }
    except Exception as e:
        logger.warning(f"esummary failure for UID '{uid}': {str(e)}")
        return {}


def identify_taxonomy_duplicates() -> Dict:
    """
    Identify duplicate taxonomy entries based on taxID and consistent information.
    
    Duplicates are defined as:
    - Multiple entries with the same taxID (not NULL)
    - OR entries with consistent information across all non-NULL columns
      (excluding potentially variable fields: species, section, species_group, subgenus, strain)
    
    Returns:
        Dict: Report of duplicate groups found
    """
    logger.info("Identifying taxonomy duplicates...")
    session = StarbaseSession()
    
    duplicates = {
        'taxid_duplicates': [],
        'consistent_duplicates': [],
        'summary': {}
    }
    
    try:
        # Find duplicates by taxID (excluding NULL taxIDs)
        taxid_groups = session.query(
            Taxonomy.taxID,
            func.count(Taxonomy.id).label('count')
        ).filter(
            Taxonomy.taxID.isnot(None),
            Taxonomy.taxID != ''
        ).group_by(Taxonomy.taxID).having(
            func.count(Taxonomy.id) > 1
        ).all()
        
        for taxid, count in taxid_groups:
            # Get all entries for this taxID
            entries = session.query(Taxonomy).filter(Taxonomy.taxID == taxid).order_by(Taxonomy.id).all()
            
            # Check if entries are consistent across non-variable fields
            primary_entry = entries[0]  # Lowest ID (oldest)
            consistent = True
            inconsistent_fields = []
            
            # Fields to check for consistency (excluding variable fields)
            check_fields = [
                'name', 'superkingdom', 'clade', 'kingdom', 'subkingdom', 
                'phylum', 'subphylum', 'class_', 'subclass', 'order', 
                'suborder', 'family', 'genus'
            ]
            
            for entry in entries[1:]:
                for field in check_fields:
                    primary_val = getattr(primary_entry, field)
                    entry_val = getattr(entry, field)
                    
                    # Only compare if both values are not None/empty
                    if (primary_val and primary_val.strip() and 
                        entry_val and entry_val.strip() and 
                        primary_val.strip() != entry_val.strip()):
                        consistent = False
                        inconsistent_fields.append({
                            'field': field,
                            'primary_value': primary_val,
                            'conflicting_value': entry_val,
                            'conflicting_id': entry.id
                        })
            
            duplicates['taxid_duplicates'].append({
                'taxID': taxid,
                'count': count,
                'primary_id': primary_entry.id,
                'duplicate_ids': [e.id for e in entries[1:]],
                'consistent': consistent,
                'inconsistent_fields': inconsistent_fields,
                'entries': [
                    {
                        'id': e.id,
                        'name': e.name,
                        'taxID': e.taxID,
                        'genus': e.genus,
                        'species': e.species,
                        'section': e.section
                    } for e in entries
                ]
            })
        
        # Find potential semantic duplicates (same core taxonomy but different/missing taxID)
        # Group by (name, genus, superkingdom, kingdom, phylum, family) combination
        all_taxonomies = session.query(Taxonomy).all()
        semantic_groups = {}
        
        for tax in all_taxonomies:
            # Create a key from core taxonomic fields
            key_parts = []
            for field in ['name', 'genus', 'superkingdom', 'kingdom', 'phylum', 'family']:
                val = getattr(tax, field)
                key_parts.append((val or '').strip().lower())
            key = tuple(key_parts)
            
            if key not in semantic_groups:
                semantic_groups[key] = []
            semantic_groups[key].append(tax)
        
        # Find groups with multiple entries
        for key, entries in semantic_groups.items():
            if len(entries) > 1:
                # Skip if all entries have the same taxID (already covered above)
                taxids = set((e.taxID or '').strip() for e in entries)
                if len(taxids) <= 1:
                    continue
                
                # Check if entries are consistent
                primary_entry = min(entries, key=lambda x: x.id)  # Lowest ID
                consistent = True
                inconsistent_fields = []
                
                check_fields = [
                    'name', 'superkingdom', 'clade', 'kingdom', 'subkingdom', 
                    'phylum', 'subphylum', 'class_', 'subclass', 'order', 
                    'suborder', 'family', 'genus'
                ]
                
                for entry in entries:
                    if entry.id == primary_entry.id:
                        continue
                    for field in check_fields:
                        primary_val = getattr(primary_entry, field)
                        entry_val = getattr(entry, field)
                        
                        # Only compare if both values are not None/empty
                        if (primary_val and primary_val.strip() and 
                            entry_val and entry_val.strip() and 
                            primary_val.strip() != entry_val.strip()):
                            consistent = False
                            inconsistent_fields.append({
                                'field': field,
                                'primary_value': primary_val,
                                'conflicting_value': entry_val,
                                'conflicting_id': entry.id
                            })
                
                if consistent:  # Only report if consistent
                    duplicates['consistent_duplicates'].append({
                        'key': str(key),
                        'count': len(entries),
                        'primary_id': primary_entry.id,
                        'duplicate_ids': [e.id for e in entries if e.id != primary_entry.id],
                        'taxids': list(taxids),
                        'entries': [
                            {
                                'id': e.id,
                                'name': e.name,
                                'taxID': e.taxID,
                                'genus': e.genus,
                                'species': e.species,
                                'section': e.section
                            } for e in entries
                        ]
                    })
        
        # Summary
        total_taxid_duplicates = sum(d['count'] - 1 for d in duplicates['taxid_duplicates'])
        total_consistent_duplicates = sum(d['count'] - 1 for d in duplicates['consistent_duplicates'])
        
        duplicates['summary'] = {
            'taxid_duplicate_groups': len(duplicates['taxid_duplicates']),
            'total_taxid_duplicates': total_taxid_duplicates,
            'consistent_duplicate_groups': len(duplicates['consistent_duplicates']),
            'total_consistent_duplicates': total_consistent_duplicates,
            'total_entries_to_remove': total_taxid_duplicates + total_consistent_duplicates
        }
        
        logger.info(f"Found {len(duplicates['taxid_duplicates'])} taxID duplicate groups ({total_taxid_duplicates} entries to remove)")
        logger.info(f"Found {len(duplicates['consistent_duplicates'])} consistent duplicate groups ({total_consistent_duplicates} entries to remove)")
        
    except Exception as e:
        logger.error(f"Error identifying taxonomy duplicates: {str(e)}")
        raise
    finally:
        session.close()
    
    return duplicates


def consolidate_taxonomy_duplicates(dry_run: bool = True) -> Dict:
    """
    Consolidate duplicate taxonomy entries by:
    1. Keeping the entry with the lowest ID (first created)
    2. Updating all references to point to the retained entry
    3. Deleting the duplicate entries
    
    Args:
        dry_run (bool): If True, only analyze and report what would be done
        
    Returns:
        Dict: Report of consolidations performed
    """
    logger.info("Consolidating taxonomy duplicates...")
    session = StarbaseSession()
    
    consolidation = {
        'entries_consolidated': [],
        'references_updated': [],
        'entries_deleted': [],
        'summary': {}
    }
    
    try:
        # First identify duplicates
        duplicates_report = identify_taxonomy_duplicates()
        
        total_updated_refs = 0
        total_deleted = 0
        
        # Process taxID duplicates
        for duplicate_group in duplicates_report['taxid_duplicates']:
            if not duplicate_group['consistent']:
                logger.warning(f"Skipping inconsistent taxID group {duplicate_group['taxID']} - manual review needed")
                continue
            
            primary_id = duplicate_group['primary_id']
            duplicate_ids = duplicate_group['duplicate_ids']
            
            # Update all references to point to primary_id
            for dup_id in duplicate_ids:
                # Update genomes.taxonomy_id
                genomes_updated = session.query(Genome).filter(Genome.taxonomy_id == dup_id).all()
                for genome in genomes_updated:
                    if not dry_run:
                        genome.taxonomy_id = primary_id
                        session.add(genome)
                    
                    consolidation['references_updated'].append({
                        'table': 'genomes',
                        'record_id': genome.id,
                        'old_taxonomy_id': dup_id,
                        'new_taxonomy_id': primary_id,
                        'action': f'Updated genome {genome.id} taxonomy_id from {dup_id} to {primary_id}'
                    })
                    total_updated_refs += 1
                
                # Update joined_ships.tax_id
                joined_ships_updated = session.query(JoinedShips).filter(JoinedShips.tax_id == dup_id).all()
                for joined in joined_ships_updated:
                    if not dry_run:
                        joined.tax_id = primary_id
                        session.add(joined)
                    
                    consolidation['references_updated'].append({
                        'table': 'joined_ships',
                        'record_id': joined.id,
                        'old_taxonomy_id': dup_id,
                        'new_taxonomy_id': primary_id,
                        'action': f'Updated joined_ships {joined.id} tax_id from {dup_id} to {primary_id}'
                    })
                    total_updated_refs += 1
                
                # Delete the duplicate taxonomy entry
                if not dry_run:
                    duplicate_entry = session.query(Taxonomy).filter(Taxonomy.id == dup_id).first()
                    if duplicate_entry:
                        session.delete(duplicate_entry)
                        consolidation['entries_deleted'].append({
                            'taxonomy_id': dup_id,
                            'name': duplicate_entry.name,
                            'taxID': duplicate_entry.taxID,
                            'action': f'Deleted duplicate taxonomy entry {dup_id}'
                        })
                        total_deleted += 1
                else:
                    duplicate_entry = session.query(Taxonomy).filter(Taxonomy.id == dup_id).first()
                    if duplicate_entry:
                        consolidation['entries_deleted'].append({
                            'taxonomy_id': dup_id,
                            'name': duplicate_entry.name,
                            'taxID': duplicate_entry.taxID,
                            'action': f'Would delete duplicate taxonomy entry {dup_id}'
                        })
                        total_deleted += 1
            
            consolidation['entries_consolidated'].append({
                'taxID': duplicate_group['taxID'],
                'primary_id': primary_id,
                'duplicate_ids': duplicate_ids,
                'duplicates_removed': len(duplicate_ids),
                'references_updated': len([r for r in consolidation['references_updated'] 
                                         if r['new_taxonomy_id'] == primary_id]),
                'action': f'Consolidated {len(duplicate_ids)} duplicates into taxonomy {primary_id}'
            })
        
        # Process consistent semantic duplicates
        for duplicate_group in duplicates_report['consistent_duplicates']:
            primary_id = duplicate_group['primary_id']
            duplicate_ids = duplicate_group['duplicate_ids']
            
            # Update all references to point to primary_id
            for dup_id in duplicate_ids:
                # Update genomes.taxonomy_id
                genomes_updated = session.query(Genome).filter(Genome.taxonomy_id == dup_id).all()
                for genome in genomes_updated:
                    if not dry_run:
                        genome.taxonomy_id = primary_id
                        session.add(genome)
                    
                    consolidation['references_updated'].append({
                        'table': 'genomes',
                        'record_id': genome.id,
                        'old_taxonomy_id': dup_id,
                        'new_taxonomy_id': primary_id,
                        'action': f'Updated genome {genome.id} taxonomy_id from {dup_id} to {primary_id}'
                    })
                    total_updated_refs += 1
                
                # Update joined_ships.tax_id
                joined_ships_updated = session.query(JoinedShips).filter(JoinedShips.tax_id == dup_id).all()
                for joined in joined_ships_updated:
                    if not dry_run:
                        joined.tax_id = primary_id
                        session.add(joined)
                    
                    consolidation['references_updated'].append({
                        'table': 'joined_ships',
                        'record_id': joined.id,
                        'old_taxonomy_id': dup_id,
                        'new_taxonomy_id': primary_id,
                        'action': f'Updated joined_ships {joined.id} tax_id from {dup_id} to {primary_id}'
                    })
                    total_updated_refs += 1
                
                # Delete the duplicate taxonomy entry
                if not dry_run:
                    duplicate_entry = session.query(Taxonomy).filter(Taxonomy.id == dup_id).first()
                    if duplicate_entry:
                        session.delete(duplicate_entry)
                        consolidation['entries_deleted'].append({
                            'taxonomy_id': dup_id,
                            'name': duplicate_entry.name,
                            'taxID': duplicate_entry.taxID,
                            'action': f'Deleted duplicate taxonomy entry {dup_id}'
                        })
                        total_deleted += 1
                else:
                    duplicate_entry = session.query(Taxonomy).filter(Taxonomy.id == dup_id).first()
                    if duplicate_entry:
                        consolidation['entries_deleted'].append({
                            'taxonomy_id': dup_id,
                            'name': duplicate_entry.name,
                            'taxID': duplicate_entry.taxID,
                            'action': f'Would delete duplicate taxonomy entry {dup_id}'
                        })
                        total_deleted += 1
            
            consolidation['entries_consolidated'].append({
                'key': duplicate_group['key'],
                'primary_id': primary_id,
                'duplicate_ids': duplicate_ids,
                'duplicates_removed': len(duplicate_ids),
                'references_updated': len([r for r in consolidation['references_updated'] 
                                         if r['new_taxonomy_id'] == primary_id]),
                'action': f'Consolidated {len(duplicate_ids)} semantic duplicates into taxonomy {primary_id}'
            })
        
        if not dry_run:
            session.commit()
            logger.info(f"Consolidated {len(consolidation['entries_consolidated'])} duplicate groups")
            logger.info(f"Updated {total_updated_refs} references")
            logger.info(f"Deleted {total_deleted} duplicate entries")
        else:
            logger.info(f"Would consolidate {len(consolidation['entries_consolidated'])} duplicate groups")
            logger.info(f"Would update {total_updated_refs} references")
            logger.info(f"Would delete {total_deleted} duplicate entries")
        
        consolidation['summary'] = {
            'groups_consolidated': len(consolidation['entries_consolidated']),
            'total_references_updated': total_updated_refs,
            'total_entries_deleted': total_deleted,
            'recommendation': 'Review inconsistent groups manually before processing'
        }
        
    except Exception as e:
        logger.error(f"Error consolidating taxonomy duplicates: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return consolidation


def _find_or_create_genome_from_ome(ome_code: str, ome_map: Dict, session) -> int:
    """
    Find existing genome by ome code or create new one using OME map data.
    
    Args:
        ome_code (str): The OME code to find/create
        ome_map (Dict): OME mapping data from _parse_ome_map()
        session: Database session
        
    Returns:
        int: Genome ID (existing or newly created)
    """
    # First try to find existing genome
    existing_genome = session.query(Genome).filter(Genome.ome == ome_code).first()
    if existing_genome:
        return existing_genome.id
    
    # Check if we have OME map data for this code
    if ome_code not in ome_map:
        raise ValueError(f"No OME map data found for {ome_code}")
    
    ome_data = ome_map[ome_code]
    
    # Create or find taxonomy entry
    taxonomy_id = _find_or_create_taxonomy(
        name=ome_data['name'],
        taxid=None,  # OME map doesn't typically have taxIDs
        session=session
    )
    
    # Update taxonomy with additional fields from OME map
    taxonomy = session.query(Taxonomy).filter(Taxonomy.id == taxonomy_id).first()
    if taxonomy:
        if ome_data['genus'] and not taxonomy.genus:
            taxonomy.genus = ome_data['genus']
        if ome_data['species_epithet'] and not taxonomy.species:
            taxonomy.species = ome_data['species_epithet']
        session.add(taxonomy)
    
    # Create new genome entry
    new_genome = Genome(
        ome=ome_code,
        taxonomy_id=taxonomy_id,
        assembly_accession=ome_data.get('assembly_accession') or None
    )
    session.add(new_genome)
    session.flush()  # Get the ID
    
    return new_genome.id


def fix_missing_tax_id_via_ome_consistency(dry_run: bool = True, ome_map_path: str = None) -> Dict:
    """
    Fix missing tax_id in joined_ships table by using ome code consistency.
    
    For entries with starshipID that have "ome" codes as prefix (e.g., 'altbur1_sequence1'),
    extract the ome code and ensure all entries with the same ome code have consistent tax_id.
    
    If no existing tax_id is available and ome_map_path is provided, will create
    missing genome and taxonomy entries from the OME map data.
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        ome_map_path (str): Optional path to OME mapping file for creating missing entries
        
    Returns:
        Dict: Report of fixes applied
    """
    logger.info("Fixing missing tax_id in joined_ships using ome code consistency...")
    session = StarbaseSession()
    
    fixes = {
        'tax_ids_filled': [],
        'ome_groups_analyzed': [],
        'inconsistent_ome_groups': [],
        'genomes_created': [],
        'taxonomies_created': [],
        'warnings': [],
        'summary': {}
    }
    
    # Load OME map if provided
    ome_map = {}
    if ome_map_path:
        try:
            ome_map = _parse_ome_map(ome_map_path)
            logger.info(f"Loaded {len(ome_map)} OME mappings from {ome_map_path}")
        except Exception as e:
            logger.error(f"Failed to load OME map from {ome_map_path}: {str(e)}")
            raise
    
    try:
        # Get all joined_ships entries with starshipIDs that look like ome codes
        all_joined_ships = session.query(JoinedShips).all()
        
        # Group by ome code extracted from starshipID
        ome_groups = {}
        
        for js in all_joined_ships:
            if not js.starshipID:
                continue
                
            # Extract ome code - look for pattern like 'abcdef1_' at start
            import re
            ome_match = re.match(r'^([a-z]{6}\d+)_', js.starshipID.lower())
            if not ome_match:
                continue
                
            ome_code = ome_match.group(1)
            
            if ome_code not in ome_groups:
                ome_groups[ome_code] = []
            ome_groups[ome_code].append(js)
        
        # Analyze each ome group for tax_id consistency
        for ome_code, entries in ome_groups.items():
            if len(entries) < 2:
                continue  # Skip single entries
            
            # Collect tax_ids for this ome group
            tax_ids = []
            entries_with_tax_id = []
            entries_without_tax_id = []
            
            for entry in entries:
                if entry.tax_id:
                    tax_ids.append(entry.tax_id)
                    entries_with_tax_id.append(entry)
                else:
                    entries_without_tax_id.append(entry)
            
            # Check consistency
            unique_tax_ids = list(set(tax_ids))
            
            if len(unique_tax_ids) == 0:
                # No tax_id assignments for this ome group
                # Try to create genome and taxonomy from OME map if available
                if ome_map and ome_code in ome_map and not dry_run:
                    try:
                        # Find or create genome entry
                        genome_id = _find_or_create_genome_from_ome(ome_code, ome_map, session)
                        
                        # Get the taxonomy_id from the genome
                        genome = session.query(Genome).filter(Genome.id == genome_id).first()
                        if genome and genome.taxonomy_id:
                            consensus_tax_id = genome.taxonomy_id
                            
                            # Fill tax_id for all entries in this ome group
                            for entry in entries_without_tax_id:
                                # Update joined_ships to link to the genome
                                entry.tax_id = consensus_tax_id
                                entry.genome_id = genome_id
                                session.add(entry)
                                
                                fixes['tax_ids_filled'].append({
                                    'joined_ships_id': entry.id,
                                    'starshipID': entry.starshipID,
                                    'ome_code': ome_code,
                                    'assigned_tax_id': consensus_tax_id,
                                    'genome_id': genome_id,
                                    'action': f'Created genome and taxonomy for ome {ome_code}, assigned tax_id {consensus_tax_id}'
                                })
                            
                            # Track creation
                            if genome_id not in [g['genome_id'] for g in fixes['genomes_created']]:
                                fixes['genomes_created'].append({
                                    'genome_id': genome_id,
                                    'ome_code': ome_code,
                                    'taxonomy_id': consensus_tax_id,
                                    'name': ome_map[ome_code]['name']
                                })
                            
                            fixes['ome_groups_analyzed'].append({
                                'ome_code': ome_code,
                                'total_entries': len(entries),
                                'with_tax_id': 0,
                                'without_tax_id': len(entries_without_tax_id),
                                'consensus_tax_id': consensus_tax_id,
                                'genome_id': genome_id,
                                'filled': len(entries_without_tax_id),
                                'status': 'created_from_ome_map'
                            })
                            continue
                        else:
                            logger.warning(f"Created genome {genome_id} for ome {ome_code} but no taxonomy_id assigned")
                    except Exception as e:
                        logger.error(f"Failed to create genome/taxonomy for ome {ome_code}: {str(e)}")
                        fixes['warnings'].append(f"Failed to create entries for ome {ome_code}: {str(e)}")
                
                # No tax_id available and either no OME map or creation failed
                fixes['ome_groups_analyzed'].append({
                    'ome_code': ome_code,
                    'total_entries': len(entries),
                    'with_tax_id': 0,
                    'without_tax_id': len(entries_without_tax_id),
                    'status': 'no_tax_id_available' + (' (no_ome_map)' if not ome_map else ' (creation_failed)' if ome_code in ome_map else ' (not_in_ome_map)')
                })
                continue
                
            elif len(unique_tax_ids) == 1:
                # Consistent tax_id - fill in missing ones
                consensus_tax_id = unique_tax_ids[0]
                
                for entry in entries_without_tax_id:
                    if not dry_run:
                        entry.tax_id = consensus_tax_id
                        session.add(entry)
                    
                    fixes['tax_ids_filled'].append({
                        'joined_ships_id': entry.id,
                        'starshipID': entry.starshipID,
                        'ome_code': ome_code,
                        'assigned_tax_id': consensus_tax_id,
                        'action': f'Assigned tax_id {consensus_tax_id} to {entry.starshipID} based on ome code {ome_code}'
                    })
                
                fixes['ome_groups_analyzed'].append({
                    'ome_code': ome_code,
                    'total_entries': len(entries),
                    'with_tax_id': len(entries_with_tax_id),
                    'without_tax_id': len(entries_without_tax_id),
                    'consensus_tax_id': consensus_tax_id,
                    'filled': len(entries_without_tax_id),
                    'status': 'consistent'
                })
                
            else:
                # Inconsistent tax_ids - needs manual review
                fixes['inconsistent_ome_groups'].append({
                    'ome_code': ome_code,
                    'total_entries': len(entries),
                    'tax_ids': unique_tax_ids,
                    'entries_by_tax_id': {
                        tax_id: [e.starshipID for e in entries_with_tax_id if e.tax_id == tax_id]
                        for tax_id in unique_tax_ids
                    },
                    'entries_without_tax_id': [e.starshipID for e in entries_without_tax_id],
                    'status': 'inconsistent'
                })
                
                fixes['warnings'].append(f"OME code {ome_code} has inconsistent tax_ids: {unique_tax_ids}")
        
        if not dry_run:
            session.commit()
            logger.info(f"Applied tax_id fixes for {len(fixes['tax_ids_filled'])} entries")
        else:
            logger.info(f"Would apply tax_id fixes for {len(fixes['tax_ids_filled'])} entries")
        
        fixes['summary'] = {
            'ome_groups_found': len(ome_groups),
            'ome_groups_analyzed': len(fixes['ome_groups_analyzed']),
            'tax_ids_filled': len(fixes['tax_ids_filled']),
            'inconsistent_groups': len(fixes['inconsistent_ome_groups']),
            'genomes_created': len(fixes['genomes_created']),
            'warnings': len(fixes['warnings']),
            'ome_map_used': bool(ome_map),
            'recommendation': 'Review inconsistent ome groups manually' + ('; Use --ome-map to create missing genomes/taxonomies' if not ome_map else '')
        }
        
        logger.info(f"Analyzed {len(ome_groups)} ome groups")
        logger.info(f"Filled {len(fixes['tax_ids_filled'])} missing tax_ids")
        logger.info(f"Found {len(fixes['inconsistent_ome_groups'])} inconsistent ome groups")
        
    except Exception as e:
        logger.error(f"Error fixing tax_id via ome consistency: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return fixes


def lookup_assembly_accessions_from_contigs(email: str = None, api_key: str = None, limit: int = 0, dry_run: bool = True) -> Dict:
    """
    For each unique contigID in starship_features, remove any leading 'ome' prefix, attempt to resolve
    an assembly accession via NCBI E-utilities (esearch+esummary), and record results to cleanup_issues.

    Args:
        email: Optional contact email for NCBI etiquette
        api_key: Optional API key
        limit: Optional cap on number of distinct contigIDs to query (0 = no limit)
        dry_run: If True, do not write to cleanup_issues; just return a summary

    Returns:
        Dict summary with counts and first few mappings
    """
    logger.info("Looking up assembly accessions from contigIDs via NCBI esearch...")
    session = StarbaseSession()

    summary = {
        'total_contigs': 0,
        'queried': 0,
        'resolved': 0,
        'skipped_empty': 0,
        'examples': []
    }

    try:
        # Collect unique contigIDs
        contigs = session.query(StarshipFeatures.contigID).filter(
            StarshipFeatures.contigID.isnot(None),
            StarshipFeatures.contigID != ''
        ).distinct().all()
        unique_contigs = [row[0] for row in contigs]
        summary['total_contigs'] = len(unique_contigs)

        # Prepare cleanup_issues table
        create_cleanup_issues_table()

        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()
        insert_sql = (
            """
            INSERT OR IGNORE INTO cleanup_issues
            (issue_type, category, table_name, record_id, details, status, source)
            VALUES (?, ?, ?, ?, ?, 'OPEN', ?)
            """
        )

        queried = 0
        resolved = 0
        for contig in unique_contigs:
            if limit and queried >= limit:
                break
            term_raw = _normalize_contig_id(contig)
            if not term_raw:
                summary['skipped_empty'] += 1
                continue
            queried += 1

            uid = _ncbi_esearch_assembly(term_raw, email=email, api_key=api_key)
            meta = _ncbi_assembly_esummary(uid, email=email, api_key=api_key) if uid else {}

            details_obj = {
                'contigID': contig,
                'normalized_term': term_raw,
                'assembly_uid': uid,
                'assemblyaccession': meta.get('assemblyaccession') or '',
                'organism': meta.get('organism') or '',
                'bioproject': meta.get('bioproject') or '',
                'biosample': meta.get('biosample') or ''
            }

            if details_obj['assemblyaccession']:
                resolved += 1

            if not dry_run:
                try:
                    cursor.execute(
                        insert_sql,
                        (
                            'assembly_mapping',
                            'assembly_lookup',
                            'starship_features',
                            None,
                            json.dumps(details_obj, ensure_ascii=False, sort_keys=True),
                            'pipeline'
                        )
                    )
                except Exception as e:
                    logger.warning(f"Failed to record assembly mapping for contig '{contig}': {str(e)}")
            else:
                # Capture a few examples in dry-run mode
                if len(summary['examples']) < 5:
                    summary['examples'].append(details_obj)

        if not dry_run:
            cursor.close()
            session.commit()
        else:
            cursor.close()

        summary['queried'] = queried
        summary['resolved'] = resolved
        logger.info(f"Assembly lookup complete: queried={queried}, resolved={resolved}")
        return summary

    except Exception as e:
        logger.error(f"Error during assembly lookup: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def remove_suffixed_joined_ships_duplicates(dry_run: bool = True) -> Dict:
    r"""
    Remove duplicate joined_ships entries that have _ACC_XXXX or _SHIP_XXXX suffixes in their starshipID.
    
    Strategy:
    1. Identify duplicates by finding starshipIDs with _ACC_\d+ or _SHIP_\d+ suffixes
    2. Group duplicates with their base starshipID (without suffix)
    3. For each group:
       - Keep the entry with the lowest ID (the original)
       - Coalesce NULL values in original with non-NULL values from duplicates for:
         accession_id, ship_id, ship_family_id, tax_id, genome_id, captain_id, 
         ship_navis_id, ship_haplotype_id
       - Delete the duplicate entries
    
    Args:
        dry_run (bool): If True, only analyze and report what would be done
        
    Returns:
        Dict: Report of duplicates removed and values coalesced
    """
    logger.info("Removing suffixed duplicate joined_ships entries...")
    session = StarbaseSession()
    
    report = {
        'duplicates_found': [],
        'entries_updated': [],
        'entries_deleted': [],
        'summary': {}
    }
    
    try:
        # Get all joined_ships entries
        all_entries = session.query(JoinedShips).all()
        
        # Group entries by base starshipID
        groups = {}
        suffix_pattern = re.compile(r'^(.+?)_(ACC|SHIP)_(\d+)$')
        
        for entry in all_entries:
            if not entry.starshipID:
                continue
            
            # Check if this starshipID has a suffix
            match = suffix_pattern.match(entry.starshipID)
            
            if match:
                # This is a duplicate with suffix
                base_starship_id = match.group(1)
                suffix_type = match.group(2)
                suffix_num = match.group(3)
                
                if base_starship_id not in groups:
                    groups[base_starship_id] = {
                        'base': None,
                        'duplicates': []
                    }
                
                groups[base_starship_id]['duplicates'].append({
                    'entry': entry,
                    'suffix_type': suffix_type,
                    'suffix_num': suffix_num
                })
            else:
                # This might be an original entry
                if entry.starshipID not in groups:
                    groups[entry.starshipID] = {
                        'base': None,
                        'duplicates': []
                    }
                
                # Store as base if we don't have one yet, or if this has a lower ID
                if (groups[entry.starshipID]['base'] is None or 
                    entry.id < groups[entry.starshipID]['base'].id):
                    groups[entry.starshipID]['base'] = entry
        
        # Process each group
        total_updated = 0
        total_deleted = 0
        
        for base_starship_id, group_data in groups.items():
            if not group_data['duplicates']:
                continue  # No duplicates for this starshipID
            
            base_entry = group_data['base']
            duplicates = group_data['duplicates']
            
            # If no base entry found, use the duplicate with lowest ID as base
            if base_entry is None:
                all_in_group = [d['entry'] for d in duplicates]
                all_in_group.sort(key=lambda x: x.id)
                base_entry = all_in_group[0]
                duplicates = [{'entry': e, 'suffix_type': 'UNKNOWN', 'suffix_num': ''} 
                             for e in all_in_group[1:]]
            
            # Fields to coalesce
            fields_to_coalesce = [
                'accession_id', 'ship_id', 'ship_family_id', 'tax_id', 
                'genome_id', 'captain_id', 'ship_navis_id', 'ship_haplotype_id'
            ]
            
            # Track what we're updating
            fields_updated = {}
            
            # Coalesce NULL values from base with non-NULL values from duplicates
            for field in fields_to_coalesce:
                base_value = getattr(base_entry, field)
                
                if base_value is None:
                    # Look for a non-NULL value in duplicates
                    for dup_data in duplicates:
                        dup_entry = dup_data['entry']
                        dup_value = getattr(dup_entry, field)
                        
                        if dup_value is not None:
                            if not dry_run:
                                setattr(base_entry, field, dup_value)
                                session.add(base_entry)
                            
                            fields_updated[field] = {
                                'old_value': base_value,
                                'new_value': dup_value,
                                'source_id': dup_entry.id,
                                'source_starship_id': dup_entry.starshipID
                            }
                            break  # Use the first non-NULL value found
            
            # Record the update if any fields were coalesced
            if fields_updated:
                report['entries_updated'].append({
                    'joined_ships_id': base_entry.id,
                    'starshipID': base_entry.starshipID,
                    'fields_updated': fields_updated,
                    'action': f'Coalesced {len(fields_updated)} NULL fields from duplicates'
                })
                total_updated += 1
            
            # Delete the duplicate entries
            for dup_data in duplicates:
                dup_entry = dup_data['entry']
                
                if not dry_run:
                    session.delete(dup_entry)
                
                report['entries_deleted'].append({
                    'joined_ships_id': dup_entry.id,
                    'starshipID': dup_entry.starshipID,
                    'suffix_type': dup_data['suffix_type'],
                    'suffix_num': dup_data['suffix_num'],
                    'base_id': base_entry.id,
                    'base_starshipID': base_entry.starshipID,
                    'action': f'Deleted duplicate entry (merged into {base_entry.id})'
                })
                total_deleted += 1
            
            # Record the duplicate group
            report['duplicates_found'].append({
                'base_starship_id': base_starship_id,
                'base_id': base_entry.id,
                'duplicate_count': len(duplicates),
                'duplicate_ids': [d['entry'].id for d in duplicates],
                'fields_coalesced': list(fields_updated.keys()) if fields_updated else []
            })
        
        if not dry_run:
            session.commit()
            logger.info(f"Removed {total_deleted} duplicate entries, updated {total_updated} base entries")
        else:
            logger.info(f"Would remove {total_deleted} duplicate entries, would update {total_updated} base entries")
        
        report['summary'] = {
            'duplicate_groups_found': len(report['duplicates_found']),
            'entries_updated': total_updated,
            'entries_deleted': total_deleted,
            'total_fields_coalesced': sum(len(e['fields_updated']) for e in report['entries_updated'])
        }
        
        logger.info(f"Found {len(report['duplicates_found'])} duplicate groups")
        logger.info(f"Updated {total_updated} base entries with coalesced values")
        logger.info(f"Deleted {total_deleted} duplicate entries")
        
    except Exception as e:
        logger.error(f"Error removing suffixed duplicates: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return report


def fix_ships_primary_key_issues(dry_run: bool = True) -> Dict:
    """
    Fix ships table primary key issues including NULL ids and missing entries.
    
    This function:
    1. Assigns proper sequential IDs to ships with NULL id
    2. Ensures all joined_ships entries have corresponding ships entries
    3. Creates missing ships entries if needed (with empty sequence - can be filled later)
    
    Args:
        dry_run (bool): If True, only analyze and report what would be fixed
        
    Returns:
        Dict: Report of fixes applied
    """
    logger.info("Fixing ships table primary key issues...")
    session = StarbaseSession()
    
    report = {
        'null_ids_fixed': [],
        'ships_created': [],
        'joined_ships_linked': [],
        'warnings': [],
        'summary': {}
    }
    
    try:
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()
        
        # Step 1: Find the current max ID
        cursor.execute("SELECT MAX(id) FROM ships WHERE id IS NOT NULL")
        max_id_result = cursor.fetchone()
        next_id = (max_id_result[0] or 0) + 1 if max_id_result else 1
        
        logger.info(f"Current max ship id: {max_id_result[0] if max_id_result else 0}, starting new IDs from {next_id}")
        
        # Step 2: Fix ships with NULL id
        cursor.execute("SELECT rowid, accession_id, sequence, md5 FROM ships WHERE id IS NULL")
        null_id_ships = cursor.fetchall()
        
        logger.info(f"Found {len(null_id_ships)} ships with NULL id")
        
        for rowid, accession_id, sequence, md5 in null_id_ships:
            if not dry_run:
                # Assign a new ID to this ship
                cursor.execute("UPDATE ships SET id = ? WHERE rowid = ?", (next_id, rowid))
                
                report['null_ids_fixed'].append({
                    'rowid': rowid,
                    'new_id': next_id,
                    'accession_id': accession_id,
                    'has_sequence': sequence is not None and sequence != '',
                    'has_md5': md5 is not None and md5 != '',
                    'action': f'Assigned id {next_id} to ship with rowid {rowid}'
                })
            else:
                report['null_ids_fixed'].append({
                    'rowid': rowid,
                    'new_id': next_id,
                    'accession_id': accession_id,
                    'has_sequence': sequence is not None and sequence != '',
                    'has_md5': md5 is not None and md5 != '',
                    'action': f'Would assign id {next_id} to ship with rowid {rowid}'
                })
            
            next_id += 1
        
        # Step 3: Check for joined_ships without corresponding ships entries
        # Get all unique ship_id values from joined_ships that should exist in ships
        cursor.execute("""
            SELECT DISTINCT js.ship_id 
            FROM joined_ships js
            WHERE js.ship_id IS NOT NULL
            AND NOT EXISTS (SELECT 1 FROM ships s WHERE s.id = js.ship_id)
        """)
        missing_ship_ids = [row[0] for row in cursor.fetchall()]
        
        logger.info(f"Found {len(missing_ship_ids)} joined_ships entries pointing to non-existent ships")
        
        for ship_id in missing_ship_ids:
            # Get info from joined_ships to help create the ship
            cursor.execute("""
                SELECT accession_id, starshipID
                FROM joined_ships
                WHERE ship_id = ?
                LIMIT 1
            """, (ship_id,))
            js_info = cursor.fetchone()
            
            if js_info:
                accession_id, starshipID = js_info
                
                if not dry_run:
                    # Create a new ship entry
                    cursor.execute("""
                        INSERT INTO ships (id, accession_id, sequence, md5)
                        VALUES (?, ?, '', NULL)
                    """, (ship_id, accession_id))
                    
                    report['ships_created'].append({
                        'ship_id': ship_id,
                        'accession_id': accession_id,
                        'starshipID': starshipID,
                        'action': f'Created ship {ship_id} (empty sequence, can be filled later)'
                    })
                else:
                    report['ships_created'].append({
                        'ship_id': ship_id,
                        'accession_id': accession_id,
                        'starshipID': starshipID,
                        'action': f'Would create ship {ship_id} (empty sequence)'
                    })
        
        # Step 4: Find joined_ships without ship_id that need a ship created
        # NOTE: ships.accession_id is deprecated - joined_ships.accession_id is the primary link
        # We only use ships.accession_id temporarily to find existing ships, but new ships 
        # won't have accession_id set (it's redundant)
        cursor.execute("""
            SELECT js.id, js.starshipID, js.accession_id
            FROM joined_ships js
            WHERE js.ship_id IS NULL
            AND js.accession_id IS NOT NULL
        """)
        unlinked_joined_ships = cursor.fetchall()
        
        logger.info(f"Found {len(unlinked_joined_ships)} joined_ships entries without ship_id but with accession_id")
        
        for js_id, starshipID, accession_id in unlinked_joined_ships:
            # Try to find a ship with this accession_id (for backward compatibility with old data)
            # Note: This is temporary - ships.accession_id should be deprecated
            cursor.execute("SELECT id FROM ships WHERE accession_id = ? LIMIT 1", (accession_id,))
            ship_result = cursor.fetchone()
            
            if ship_result:
                ship_id = ship_result[0]
                
                if not dry_run:
                    cursor.execute("UPDATE joined_ships SET ship_id = ? WHERE id = ?", (ship_id, js_id))
                    
                    report['joined_ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': ship_id,
                        'accession_id': accession_id,
                        'action': f'Linked joined_ships {js_id} to existing ship {ship_id} (via ships.accession_id)'
                    })
                else:
                    report['joined_ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': ship_id,
                        'accession_id': accession_id,
                        'action': f'Would link joined_ships {js_id} to existing ship {ship_id} (via ships.accession_id)'
                    })
            else:
                # No ship exists, create one WITHOUT accession_id (following new architecture)
                # The accession link is in joined_ships.accession_id, not ships.accession_id
                if not dry_run:
                    # Use the next available ID
                    cursor.execute("""
                        INSERT INTO ships (id, sequence, md5)
                        VALUES (?, '', NULL)
                    """, (next_id,))
                    
                    cursor.execute("UPDATE joined_ships SET ship_id = ? WHERE id = ?", (next_id, js_id))
                    
                    report['ships_created'].append({
                        'ship_id': next_id,
                        'accession_id': None,  # Not storing in ships anymore
                        'starshipID': starshipID,
                        'action': f'Created ship {next_id} for joined_ships {js_id} (accession in joined_ships only)'
                    })
                    
                    report['joined_ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': next_id,
                        'accession_id': accession_id,
                        'action': f'Linked joined_ships {js_id} to newly created ship {next_id}'
                    })
                    
                    next_id += 1
                else:
                    report['ships_created'].append({
                        'ship_id': next_id,
                        'accession_id': None,  # Not storing in ships anymore
                        'starshipID': starshipID,
                        'action': f'Would create ship {next_id} for joined_ships {js_id} (accession in joined_ships only)'
                    })
                    
                    report['joined_ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': next_id,
                        'accession_id': accession_id,
                        'action': f'Would link joined_ships {js_id} to newly created ship {next_id}'
                    })
                    
                    next_id += 1
        
        cursor.close()
        
        if not dry_run:
            session.commit()
            logger.info("Applied primary key fixes to database")
        else:
            logger.info("Dry run - no changes applied")
        
        report['summary'] = {
            'null_ids_fixed': len(report['null_ids_fixed']),
            'ships_created': len(report['ships_created']),
            'joined_ships_linked': len(report['joined_ships_linked']),
            'next_available_id': next_id,
            'recommendation': 'After fixing, update ship sequences from FASTA file'
        }
        
        logger.info(f"Fixed {len(report['null_ids_fixed'])} NULL ids")
        logger.info(f"Created {len(report['ships_created'])} new ship entries")
        logger.info(f"Linked {len(report['joined_ships_linked'])} joined_ships entries")
        
    except Exception as e:
        logger.error(f"Error fixing ships primary key issues: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return report


def link_joined_ships_via_fasta(fasta_path: str, dry_run: bool = True) -> Dict:
    """
    Link joined_ships to ships by matching starshipID to FASTA headers and comparing sequences.
    
    This is the proper way to link joined_ships.ship_id:
    1. For each joined_ships without ship_id, search for starshipID in FASTA headers
    2. Extract sequence from matching FASTA entry
    3. Calculate MD5 hash and find ship with matching sequence
    4. Link joined_ships.ship_id to the matching ship.id
    5. If no ship found with that sequence, create a new one
    
    Args:
        fasta_path (str): Path to FASTA file with starship sequences
        dry_run (bool): If True, only analyze and report what would be done
        
    Returns:
        Dict: Report of links created
    """
    logger.info(f"Linking joined_ships to ships via FASTA file: {fasta_path}")
    session = StarbaseSession()
    
    report = {
        'ships_linked': [],
        'ships_created': [],
        'fasta_entries_parsed': 0,
        'fallback_fna_matches': [],
        'no_fasta_match': [],
        'warnings': [],
        'summary': {}
    }
    
    try:
        # Import sequence utilities
        try:
            from ...utils.seq_utils import clean_sequence, revcomp
            from ...utils.classification_utils import generate_md5_hash
        except ImportError:
            import sys
            import os
            sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
            from src.utils.seq_utils import clean_sequence, revcomp
            from src.utils.classification_utils import generate_md5_hash
        
        # Parse FASTA file into a dict: header -> sequence
        logger.info("Parsing FASTA file...")
        fasta_sequences = {}
        
        with open(fasta_path, 'r') as f:
            current_header = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous entry
                    if current_header and current_seq:
                        fasta_sequences[current_header] = ''.join(current_seq)
                    
                    # Start new entry
                    current_header = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last entry
            if current_header and current_seq:
                fasta_sequences[current_header] = ''.join(current_seq)
        
        report['fasta_entries_parsed'] = len(fasta_sequences)
        logger.info(f"Parsed {len(fasta_sequences)} sequences from FASTA file")
        
        # Get joined_ships entries without ship_id
        raw_conn = session.connection().connection
        cursor = raw_conn.cursor()
        
        cursor.execute("""
            SELECT id, starshipID, accession_id
            FROM joined_ships
            WHERE ship_id IS NULL
        """)
        unlinked_joined_ships = cursor.fetchall()
        
        logger.info(f"Found {len(unlinked_joined_ships)} joined_ships entries without ship_id")
        
        # Get current max ship ID for creating new entries
        cursor.execute("SELECT MAX(id) FROM ships WHERE id IS NOT NULL")
        max_id_result = cursor.fetchone()
        next_id = (max_id_result[0] or 0) + 1 if max_id_result else 1
        
        for js_id, starshipID, accession_id in unlinked_joined_ships:
            # Find matching FASTA header (partial match)
            matching_header = None
            matching_seq = None
            
            for header, seq in fasta_sequences.items():
                # Check if starshipID appears in the header
                if starshipID in header:
                    matching_header = header
                    matching_seq = seq
                    break
            
            if not matching_seq:
                # Fallback: Check for individual FASTA file in ships/fna directory
                # Note: Check both fna/ and fna/fna/ subdirectories
                base_fna_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'db', 'ships', 'fna')
                fna_dirs = [
                    base_fna_dir,
                    os.path.join(base_fna_dir, 'fna')  # Additional nested fna directory
                ]
                
                for fna_dir in fna_dirs:
                    if matching_seq:
                        break  # Already found a match
                    
                    if os.path.exists(fna_dir):
                        # Try to find a file matching the starshipID
                        for filename in os.listdir(fna_dir):
                            if starshipID in filename and (filename.endswith('.fna') or filename.endswith('.fasta') or filename.endswith('.fa')):
                                fna_path = os.path.join(fna_dir, filename)
                                try:
                                    # Read the sequence from the individual file
                                    with open(fna_path, 'r') as fna_file:
                                        fna_header = None
                                        fna_seq = []
                                        for line in fna_file:
                                            line = line.strip()
                                            if not line:
                                                continue
                                            if line.startswith('>'):
                                                if fna_header:  # Already have a header, stop
                                                    break
                                                fna_header = line[1:]
                                            else:
                                                fna_seq.append(line)
                                        
                                        if fna_seq:
                                            matching_header = f"{filename} (from {os.path.relpath(fna_dir, os.path.dirname(__file__))})"
                                            matching_seq = ''.join(fna_seq)
                                            logger.info(f"Found fallback sequence for {starshipID} in {fna_path}")
                                            report['fallback_fna_matches'].append({
                                                'joined_ships_id': js_id,
                                                'starshipID': starshipID,
                                                'filename': filename,
                                                'directory': os.path.relpath(fna_dir, os.path.dirname(__file__))
                                            })
                                            break
                                except Exception as e:
                                    logger.warning(f"Error reading {filename}: {str(e)}")
                
                if not matching_seq:
                    report['no_fasta_match'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'issue': 'No matching FASTA header found (checked main FASTA and fna directory)'
                    })
                    continue
            
            # Clean sequence and calculate MD5
            clean_seq = clean_sequence(matching_seq)
            if not clean_seq:
                report['warnings'].append({
                    'joined_ships_id': js_id,
                    'starshipID': starshipID,
                    'issue': 'Could not clean sequence from FASTA'
                })
                continue
            
            md5_hash = generate_md5_hash(clean_seq)
            md5_hash_revcomp = generate_md5_hash(revcomp(clean_seq))
            
            if not md5_hash:
                report['warnings'].append({
                    'joined_ships_id': js_id,
                    'starshipID': starshipID,
                    'issue': 'Could not generate MD5 hash'
                })
                continue
            
            # Find ship with matching MD5 (or reverse complement)
            cursor.execute("""
                SELECT id FROM ships 
                WHERE md5 = ? OR rev_comp_md5 = ? OR md5 = ? OR rev_comp_md5 = ?
                LIMIT 1
            """, (md5_hash, md5_hash, md5_hash_revcomp, md5_hash_revcomp))
            ship_result = cursor.fetchone()
            
            if ship_result:
                # Found existing ship with matching sequence
                ship_id = ship_result[0]
                
                if not dry_run:
                    cursor.execute("""
                        UPDATE joined_ships SET ship_id = ? WHERE id = ?
                    """, (ship_id, js_id))
                    
                    report['ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': ship_id,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'action': f'Linked to existing ship {ship_id} (matched by MD5)'
                    })
                else:
                    report['ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': ship_id,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'action': f'Would link to existing ship {ship_id} (matched by MD5)'
                    })
            else:
                # No ship with this sequence exists - create new one
                if not dry_run:
                    cursor.execute("""
                        INSERT INTO ships (id, sequence, md5, rev_comp_md5)
                        VALUES (?, ?, ?, ?)
                    """, (next_id, clean_seq, md5_hash, md5_hash_revcomp))
                    
                    cursor.execute("""
                        UPDATE joined_ships SET ship_id = ? WHERE id = ?
                    """, (next_id, js_id))
                    
                    report['ships_created'].append({
                        'ship_id': next_id,
                        'starshipID': starshipID,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'sequence_length': len(clean_seq),
                        'action': f'Created new ship {next_id} with sequence from FASTA'
                    })
                    
                    report['ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': next_id,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'action': f'Linked to newly created ship {next_id}'
                    })
                    
                    next_id += 1
                else:
                    report['ships_created'].append({
                        'ship_id': next_id,
                        'starshipID': starshipID,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'sequence_length': len(clean_seq),
                        'action': f'Would create new ship {next_id} with sequence from FASTA'
                    })
                    
                    report['ships_linked'].append({
                        'joined_ships_id': js_id,
                        'starshipID': starshipID,
                        'ship_id': next_id,
                        'md5': md5_hash,
                        'fasta_header': matching_header,
                        'action': f'Would link to newly created ship {next_id}'
                    })
                    
                    next_id += 1
        
        cursor.close()
        
        if not dry_run:
            session.commit()
            logger.info("Applied FASTA-based links to database")
        else:
            logger.info("Dry run - no changes applied")
        
        report['summary'] = {
            'fasta_entries': report['fasta_entries_parsed'],
            'joined_ships_processed': len(unlinked_joined_ships),
            'ships_linked': len(report['ships_linked']),
            'ships_created': len(report['ships_created']),
            'fallback_fna_matches': len(report['fallback_fna_matches']),
            'no_fasta_match': len(report['no_fasta_match']),
            'warnings': len(report['warnings']),
            'recommendation': 'Verify sequence matches and MD5 hashes are correct'
        }
        
        logger.info(f"Linked {len(report['ships_linked'])} joined_ships entries")
        logger.info(f"Created {len(report['ships_created'])} new ship entries")
        logger.info(f"No FASTA match for {len(report['no_fasta_match'])} entries")
        
    except Exception as e:
        logger.error(f"Error linking joined_ships via FASTA: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return report


def fill_missing_family_ids(dry_run: bool = True) -> Dict:
    """
    Fill in missing ship_family_id values in joined_ships table.
    
    This function uses a multi-step approach to assign family IDs:
    1. Inherit from existing classifications (if starshipID already has family elsewhere)
    2. Fill based on shared navis haplotypes (if navis/haplotype has consensus family)
    3. Classify using captain sequences (using HMMER/BLAST on captain sequences)
    
    Args:
        dry_run (bool): If True, only report what would be changed without making changes
        
    Returns:
        Dict: Report of all fixes applied or that would be applied
    """
    logger.info("Filling missing ship_family_id in joined_ships table...")
    session = StarbaseSession()
    
    report = {
        'inherited_from_existing': [],
        'filled_from_navis': [],
        'filled_from_haplotype': [],
        'filled_from_captain': [],
        'no_classification_found': [],
        'errors': [],
        'summary': {}
    }
    
    try:
        # Get all joined_ships entries without family_id
        missing_family = session.query(JoinedShips).filter(
            JoinedShips.ship_family_id.is_(None)
        ).all()
        
        logger.info(f"Found {len(missing_family)} joined_ships entries without family_id")
        
        for entry in missing_family:
            try:
                family_id = None
                source = None
                details = {}
                
                # Strategy 1: Inherit from existing classifications
                # Look for other entries with the same starshipID that have a family_id
                if entry.starshipID:
                    existing_with_family = session.query(JoinedShips).filter(
                        JoinedShips.starshipID == entry.starshipID,
                        JoinedShips.ship_family_id.isnot(None),
                        JoinedShips.id != entry.id
                    ).first()
                    
                    if existing_with_family:
                        family_id = existing_with_family.ship_family_id
                        source = 'inherited_from_existing'
                        details = {
                            'joined_ships_id': entry.id,
                            'starshipID': entry.starshipID,
                            'assigned_family_id': family_id,
                            'source_entry_id': existing_with_family.id,
                            'method': 'Inherited from another entry with same starshipID'
                        }
                        report['inherited_from_existing'].append(details)
                
                # Strategy 2: Fill based on shared navis (if they have consensus family)
                if not family_id and entry.ship_navis_id:
                    # Get all entries with the same navis
                    navis_entries = session.query(JoinedShips).filter(
                        JoinedShips.ship_navis_id == entry.ship_navis_id,
                        JoinedShips.ship_family_id.isnot(None)
                    ).all()
                    
                    if navis_entries:
                        # Count family occurrences
                        family_counts = {}
                        for nav_entry in navis_entries:
                            fam_id = nav_entry.ship_family_id
                            family_counts[fam_id] = family_counts.get(fam_id, 0) + 1
                        
                        # Get consensus family (most common)
                        consensus_family = max(family_counts, key=family_counts.get)
                        total_count = sum(family_counts.values())
                        consensus_ratio = family_counts[consensus_family] / total_count
                        
                        # Use consensus if >70% agreement
                        if consensus_ratio > 0.7:
                            family_id = consensus_family
                            source = 'filled_from_navis'
                            
                            # Get navis name for reporting
                            navis = session.query(Navis).filter(Navis.id == entry.ship_navis_id).first()
                            navis_name = navis.navis_name if navis else f"ID:{entry.ship_navis_id}"
                            
                            details = {
                                'joined_ships_id': entry.id,
                                'starshipID': entry.starshipID,
                                'assigned_family_id': family_id,
                                'navis_id': entry.ship_navis_id,
                                'navis_name': navis_name,
                                'consensus_ratio': consensus_ratio,
                                'method': f'Filled from navis consensus ({consensus_ratio:.1%} agreement)'
                            }
                            report['filled_from_navis'].append(details)
                
                # Strategy 3: Fill based on shared haplotype (if they have consensus family)
                if not family_id and entry.ship_haplotype_id:
                    # Get all entries with the same haplotype
                    haplotype_entries = session.query(JoinedShips).filter(
                        JoinedShips.ship_haplotype_id == entry.ship_haplotype_id,
                        JoinedShips.ship_family_id.isnot(None)
                    ).all()
                    
                    if haplotype_entries:
                        # Count family occurrences
                        family_counts = {}
                        for hap_entry in haplotype_entries:
                            fam_id = hap_entry.ship_family_id
                            family_counts[fam_id] = family_counts.get(fam_id, 0) + 1
                        
                        # Get consensus family (most common)
                        consensus_family = max(family_counts, key=family_counts.get)
                        total_count = sum(family_counts.values())
                        consensus_ratio = family_counts[consensus_family] / total_count
                        
                        # Use consensus if >80% agreement (higher threshold for haplotype)
                        if consensus_ratio > 0.8:
                            family_id = consensus_family
                            source = 'filled_from_haplotype'
                            
                            # Get haplotype name for reporting
                            haplotype = session.query(Haplotype).filter(
                                Haplotype.id == entry.ship_haplotype_id
                            ).first()
                            haplotype_name = haplotype.haplotype_name if haplotype else f"ID:{entry.ship_haplotype_id}"
                            
                            details = {
                                'joined_ships_id': entry.id,
                                'starshipID': entry.starshipID,
                                'assigned_family_id': family_id,
                                'haplotype_id': entry.ship_haplotype_id,
                                'haplotype_name': haplotype_name,
                                'consensus_ratio': consensus_ratio,
                                'method': f'Filled from haplotype consensus ({consensus_ratio:.1%} agreement)'
                            }
                            report['filled_from_haplotype'].append(details)
                
                # Strategy 4: Classify using captain sequences (if captain_id exists)
                if not family_id and entry.captain_id:
                    # Get the captain sequence
                    captain = session.query(Captains).filter(
                        Captains.id == entry.captain_id
                    ).first()
                    
                    if captain and captain.sequence:
                        try:
                            # Use HMMER-based family classification
                            from src.utils.classification_utils import classify_family
                            import tempfile
                            
                            # Create temporary FASTA file with captain sequence
                            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                                f.write(f">{entry.starshipID}\n{captain.sequence}\n")
                                temp_fasta = f.name
                            
                            try:
                                # Run family classification
                                family_result, _ = classify_family(
                                    fasta=temp_fasta,
                                    seq_type='prot',  # Captain sequences are protein
                                    pident_thresh=90,
                                    input_eval=0.001,
                                    threads=1
                                )
                                
                                if family_result and 'family' in family_result:
                                    family_name = family_result['family']
                                    
                                    # Look up family_id from family name
                                    family_obj = session.query(FamilyNames).filter(
                                        FamilyNames.familyName == family_name
                                    ).first()
                                    
                                    if family_obj:
                                        family_id = family_obj.id
                                        source = 'filled_from_captain'
                                        details = {
                                            'joined_ships_id': entry.id,
                                            'starshipID': entry.starshipID,
                                            'assigned_family_id': family_id,
                                            'family_name': family_name,
                                            'captain_id': entry.captain_id,
                                            'evalue': family_result.get('evalue'),
                                            'method': 'Classified using captain sequence with HMMER'
                                        }
                                        report['filled_from_captain'].append(details)
                                
                            finally:
                                # Clean up temp file
                                import os
                                if os.path.exists(temp_fasta):
                                    os.remove(temp_fasta)
                                    
                        except Exception as classify_error:
                            logger.warning(
                                f"Failed to classify captain sequence for {entry.starshipID}: {str(classify_error)}"
                            )
                            report['errors'].append({
                                'joined_ships_id': entry.id,
                                'starshipID': entry.starshipID,
                                'error': f'Captain classification failed: {str(classify_error)}'
                            })
                
                # Apply the family_id if found
                if family_id and not dry_run:
                    entry.ship_family_id = family_id
                    session.add(entry)
                    logger.debug(f"Set family_id={family_id} for joined_ships.id={entry.id} via {source}")
                elif not family_id:
                    # No classification method worked
                    report['no_classification_found'].append({
                        'joined_ships_id': entry.id,
                        'starshipID': entry.starshipID,
                        'has_navis': entry.ship_navis_id is not None,
                        'has_haplotype': entry.ship_haplotype_id is not None,
                        'has_captain': entry.captain_id is not None,
                        'reason': 'No classification method succeeded'
                    })
                    
            except Exception as e:
                logger.error(f"Error processing joined_ships.id={entry.id}: {str(e)}")
                report['errors'].append({
                    'joined_ships_id': entry.id,
                    'starshipID': entry.starshipID if entry.starshipID else 'Unknown',
                    'error': str(e)
                })
        
        # Commit changes if not dry run
        if not dry_run:
            session.commit()
            logger.info("Applied family_id assignments to database")
        else:
            logger.info("Dry run - no changes applied")
        
        # Generate summary
        report['summary'] = {
            'total_missing': len(missing_family),
            'inherited_from_existing': len(report['inherited_from_existing']),
            'filled_from_navis': len(report['filled_from_navis']),
            'filled_from_haplotype': len(report['filled_from_haplotype']),
            'filled_from_captain': len(report['filled_from_captain']),
            'no_classification_found': len(report['no_classification_found']),
            'errors': len(report['errors']),
            'total_filled': (
                len(report['inherited_from_existing']) + 
                len(report['filled_from_navis']) + 
                len(report['filled_from_haplotype']) + 
                len(report['filled_from_captain'])
            )
        }
        
        logger.info(f"Summary: {report['summary']}")
        
    except Exception as e:
        logger.error(f"Error filling missing family IDs: {str(e)}")
        if not dry_run:
            session.rollback()
        raise
    finally:
        session.close()
    
    return report
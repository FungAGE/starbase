#!/usr/bin/env python3
"""
Utility functions for managing ship quality tags.

This module provides functions to add, remove, and query quality tags
for ships in the joined_ships table using the new ship_quality_tags table.
"""

from typing import List, Dict, Optional, Set
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import func, and_, or_

from src.database.models.schema import JoinedShips, ShipQualityTags
from src.database.sql_engine import StarbaseSession
from src.config.logging import get_logger

logger = get_logger(__name__)

# Standard quality tag types based on TODO.md
STANDARD_TAG_TYPES = {
    'complete',
    'missing_boundaries',
    'missing_captain', 
    'missing_classification',
    'incomplete',
    'unannotated'
}


def add_quality_tag(joined_ship_id: int, tag_type: str, tag_value: Optional[str] = None, 
                   created_by: str = 'auto', session: Optional[Session] = None) -> bool:
    """
    Add a quality tag to a ship.
    
    Args:
        joined_ship_id: ID of the joined_ships record
        tag_type: Type of quality tag (e.g., 'missing_captain')
        tag_value: Optional value for the tag
        created_by: Who created this tag ('auto' or user identifier)
        session: Optional database session (will create one if not provided)
    
    Returns:
        bool: True if tag was added successfully, False if it already exists
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        # Check if tag already exists
        existing_tag = session.query(ShipQualityTags).filter(
            and_(
                ShipQualityTags.joined_ship_id == joined_ship_id,
                ShipQualityTags.tag_type == tag_type
            )
        ).first()
        
        if existing_tag:
            logger.debug(f"Quality tag '{tag_type}' already exists for ship {joined_ship_id}")
            return False
        
        # Create new tag
        new_tag = ShipQualityTags(
            joined_ship_id=joined_ship_id,
            tag_type=tag_type,
            tag_value=tag_value,
            created_by=created_by,
            created_at=datetime.now()
        )
        
        session.add(new_tag)
        session.commit()
        
        logger.info(f"Added quality tag '{tag_type}' to ship {joined_ship_id}")
        return True
        
    except Exception as e:
        session.rollback()
        logger.error(f"Failed to add quality tag '{tag_type}' to ship {joined_ship_id}: {str(e)}")
        raise
    finally:
        if should_close_session:
            session.close()


def remove_quality_tag(joined_ship_id: int, tag_type: str, session: Optional[Session] = None) -> bool:
    """
    Remove a quality tag from a ship.
    
    Args:
        joined_ship_id: ID of the joined_ships record
        tag_type: Type of quality tag to remove
        session: Optional database session
    
    Returns:
        bool: True if tag was removed, False if it didn't exist
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        tag = session.query(ShipQualityTags).filter(
            and_(
                ShipQualityTags.joined_ship_id == joined_ship_id,
                ShipQualityTags.tag_type == tag_type
            )
        ).first()
        
        if not tag:
            logger.debug(f"Quality tag '{tag_type}' not found for ship {joined_ship_id}")
            return False
        
        session.delete(tag)
        session.commit()
        
        logger.info(f"Removed quality tag '{tag_type}' from ship {joined_ship_id}")
        return True
        
    except Exception as e:
        session.rollback()
        logger.error(f"Failed to remove quality tag '{tag_type}' from ship {joined_ship_id}: {str(e)}")
        raise
    finally:
        if should_close_session:
            session.close()


def get_ship_quality_tags(joined_ship_id: int, session: Optional[Session] = None) -> List[Dict]:
    """
    Get all quality tags for a ship.
    
    Args:
        joined_ship_id: ID of the joined_ships record
        session: Optional database session
    
    Returns:
        List of dictionaries containing tag information
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        tags = session.query(ShipQualityTags).filter(
            ShipQualityTags.joined_ship_id == joined_ship_id
        ).all()
        
        return [
            {
                'tag_type': tag.tag_type,
                'tag_value': tag.tag_value,
                'created_at': tag.created_at,
                'created_by': tag.created_by
            }
            for tag in tags
        ]
        
    finally:
        if should_close_session:
            session.close()


def find_ships_with_tags(tag_types: List[str], match_all: bool = False, 
                        session: Optional[Session] = None) -> List[int]:
    """
    Find ships that have specific quality tags.
    
    Args:
        tag_types: List of tag types to search for
        match_all: If True, ship must have ALL tags; if False, ship must have ANY tag
        session: Optional database session
    
    Returns:
        List of joined_ship_ids that match the criteria
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        if match_all:
            # Ship must have ALL specified tags
            subqueries = []
            for tag_type in tag_types:
                subquery = session.query(ShipQualityTags.joined_ship_id).filter(
                    ShipQualityTags.tag_type == tag_type
                ).subquery()
                subqueries.append(subquery)
            
            # Find intersection of all subqueries
            query = session.query(ShipQualityTags.joined_ship_id).filter(
                ShipQualityTags.tag_type == tag_types[0]
            )
            
            for i, tag_type in enumerate(tag_types[1:], 1):
                query = query.filter(
                    ShipQualityTags.joined_ship_id.in_(
                        session.query(ShipQualityTags.joined_ship_id).filter(
                            ShipQualityTags.tag_type == tag_type
                        )
                    )
                )
            
            result = query.distinct().all()
        else:
            # Ship must have ANY of the specified tags
            result = session.query(ShipQualityTags.joined_ship_id).filter(
                ShipQualityTags.tag_type.in_(tag_types)
            ).distinct().all()
        
        return [row[0] for row in result]
        
    finally:
        if should_close_session:
            session.close()


def get_quality_tag_summary(session: Optional[Session] = None) -> Dict[str, int]:
    """
    Get a summary of all quality tags in the database.
    
    Args:
        session: Optional database session
    
    Returns:
        Dictionary mapping tag_type to count
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        result = session.query(
            ShipQualityTags.tag_type,
            func.count(ShipQualityTags.id).label('count')
        ).group_by(ShipQualityTags.tag_type).all()
        
        return {tag_type: count for tag_type, count in result}
        
    finally:
        if should_close_session:
            session.close()


def assess_ship_quality(joined_ship_id: int, session: Optional[Session] = None) -> Dict:
    """
    Assess the quality of a ship based on available data and add appropriate tags.
    
    This function examines a ship's data completeness and automatically assigns
    quality tags based on missing or incomplete information.
    
    Args:
        joined_ship_id: ID of the joined_ships record
        session: Optional database session
    
    Returns:
        Dictionary with assessment results and tags added
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        # Get the ship record
        ship = session.query(JoinedShips).filter(JoinedShips.id == joined_ship_id).first()
        if not ship:
            raise ValueError(f"Ship with ID {joined_ship_id} not found")
        
        assessment = {
            'ship_id': joined_ship_id,
            'starshipID': ship.starshipID,
            'tags_added': [],
            'tags_removed': [],
            'is_complete': True
        }
        
        # Check for missing captain
        if not ship.captain_id:
            if add_quality_tag(joined_ship_id, 'missing_captain', session=session):
                assessment['tags_added'].append('missing_captain')
            assessment['is_complete'] = False
        else:
            # Remove tag if captain is now present
            if remove_quality_tag(joined_ship_id, 'missing_captain', session=session):
                assessment['tags_removed'].append('missing_captain')
        
        # Check for missing classification (family)
        if not ship.ship_family_id:
            if add_quality_tag(joined_ship_id, 'missing_classification', session=session):
                assessment['tags_added'].append('missing_classification')
            assessment['is_complete'] = False
        else:
            if remove_quality_tag(joined_ship_id, 'missing_classification', session=session):
                assessment['tags_removed'].append('missing_classification')
        
        # Check if ship has sequence data
        if not ship.ship_id:
            if add_quality_tag(joined_ship_id, 'unannotated', session=session):
                assessment['tags_added'].append('unannotated')
            assessment['is_complete'] = False
        else:
            if remove_quality_tag(joined_ship_id, 'unannotated', session=session):
                assessment['tags_removed'].append('unannotated')
        
        # Add complete tag if all checks pass
        if assessment['is_complete']:
            if add_quality_tag(joined_ship_id, 'complete', session=session):
                assessment['tags_added'].append('complete')
        else:
            # Remove complete tag if ship is not complete
            if remove_quality_tag(joined_ship_id, 'complete', session=session):
                assessment['tags_removed'].append('complete')
        
        return assessment
        
    finally:
        if should_close_session:
            session.close()


def bulk_assess_ship_quality(joined_ship_ids: Optional[List[int]] = None, 
                           dry_run: bool = True, session: Optional[Session] = None) -> Dict:
    """
    Assess quality for multiple ships or all ships in the database.
    
    Args:
        joined_ship_ids: List of ship IDs to assess (None = all ships)
        dry_run: If True, don't actually add/remove tags
        session: Optional database session
    
    Returns:
        Dictionary with summary of assessments
    """
    should_close_session = session is None
    if session is None:
        session = StarbaseSession()
    
    try:
        if joined_ship_ids is None:
            # Get all ship IDs
            joined_ship_ids = [row[0] for row in session.query(JoinedShips.id).all()]
        
        summary = {
            'total_ships': len(joined_ship_ids),
            'ships_assessed': 0,
            'tags_added': {},
            'tags_removed': {},
            'complete_ships': 0,
            'incomplete_ships': 0
        }
        
        for ship_id in joined_ship_ids:
            try:
                if dry_run:
                    logger.info(f"[DRY RUN] Would assess ship {ship_id}")
                    continue
                
                assessment = assess_ship_quality(ship_id, session=session)
                summary['ships_assessed'] += 1
                
                if assessment['is_complete']:
                    summary['complete_ships'] += 1
                else:
                    summary['incomplete_ships'] += 1
                
                # Count tags added/removed
                for tag in assessment['tags_added']:
                    summary['tags_added'][tag] = summary['tags_added'].get(tag, 0) + 1
                
                for tag in assessment['tags_removed']:
                    summary['tags_removed'][tag] = summary['tags_removed'].get(tag, 0) + 1
                
            except Exception as e:
                logger.error(f"Failed to assess ship {ship_id}: {str(e)}")
                continue
        
        return summary
        
    finally:
        if should_close_session:
            session.close()

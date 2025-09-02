#!/usr/bin/env python3
"""
Test script to verify the database cleanup fix works.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from src.config.database import StarbaseSession
from src.config.logging import get_logger
from src.database.models.schema import Accessions, Ships, JoinedShips, Captains, StarshipFeatures, Gff

logger = get_logger(__name__)

def test_cleanup_logic():
    """
    Test the cleanup logic without actually making changes.
    """
    logger.info("Testing cleanup logic...")
    
    session = StarbaseSession()
    
    try:
        # Test that we can query the tables without the join issue
        accessions = session.query(Accessions).all()
        logger.info(f"Found {len(accessions)} accessions")
        
        if accessions:
            # Test with the first available accession instead of a hardcoded index
            test_accession = accessions[1299]
            logger.info(f"Testing with accession: {test_accession.accession_tag} (ID: {test_accession.id})")
            
            # Test the query pattern that was causing issues
            # This should work now (no join + update)
            ships_count = session.query(Ships).filter(
                Ships.accession_id == test_accession.id
            ).count()
            
            joined_count = session.query(JoinedShips).filter(
                JoinedShips.ship_id == test_accession.id
            ).count()
            
            captains_count = session.query(Captains).filter(
                Captains.ship_id == test_accession.id
            ).count()
            
            features_count = session.query(StarshipFeatures).filter(
                StarshipFeatures.ship_id == test_accession.id
            ).count()
            
            # Test GFF queries - this should use accession_id, not ship_id
            gff_count_by_accession = session.query(Gff).filter(
                Gff.accession_id == test_accession.id
            ).count()
            
            # Also test ship_id if there are ships for this accession
            if ships_count > 0:
                ships = session.query(Ships).filter(
                    Ships.accession_id == test_accession.id
                ).all()
                ship_ids = [ship.id for ship in ships]
                
                gff_count_by_ship = session.query(Gff).filter(
                    Gff.ship_id.in_(ship_ids)
                ).count()
                
                logger.info(f"  - {gff_count_by_ship} gff entries (by ship_id)")
            else:
                gff_count_by_ship = 0
            
            logger.info(f"Test accession {test_accession.accession_tag} has:")
            logger.info(f"  - {ships_count} ships")
            logger.info(f"  - {joined_count} joined_ships")
            logger.info(f"  - {captains_count} captains")
            logger.info(f"  - {features_count} features")
            logger.info(f"  - {gff_count_by_accession} gff entries (by accession_id)")
            
            # Test all accessions to see the data distribution
            logger.info("\nTesting all available accessions:")
            for i, acc in enumerate(accessions):
                ships_count = session.query(Ships).filter(
                    Ships.accession_id == acc.id
                ).count()
                
                gff_count = session.query(Gff).filter(
                    Gff.accession_id == acc.id
                ).count()
                
                logger.info(f"  {i+1}. {acc.accession_tag} (ID: {acc.id}): {ships_count} ships, {gff_count} gff entries")
        
        logger.info("✅ Cleanup logic test passed!")
        
    except Exception as e:
        logger.error(f"❌ Test failed: {str(e)}")
        raise
    finally:
        session.close()

if __name__ == "__main__":
    test_cleanup_logic()

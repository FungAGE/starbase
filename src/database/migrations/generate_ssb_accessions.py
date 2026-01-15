"""
Script to generate SSB (Ship) accessions for all existing ships in the database.

This script should be run after the ship_accessions table has been created via migration.
It will assign a unique SSB accession to each ship that doesn't already have one.

Usage:
    python -m src.database.migrations.generate_ssb_accessions
"""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from sqlalchemy import text
from src.config.database import StarbaseSession
from src.config.logging import get_logger

logger = get_logger(__name__)


def get_next_ssb_number(session):
    """Get the next available SSB accession number."""
    result = session.execute(
        text("""
        SELECT ship_accession_tag FROM ship_accessions
        WHERE ship_accession_tag LIKE 'SSB%'
        ORDER BY ship_accession_tag DESC
        LIMIT 1
        """)
    ).fetchone()
    
    if result:
        # Extract number from SSB0000001 format
        last_accession = result[0]
        last_number = int(last_accession.replace('SSB', ''))
        return last_number + 1
    else:
        return 1


def generate_ssb_accessions(dry_run=False):
    """
    Generate SSB accessions for all ships that don't have one.
    
    Args:
        dry_run (bool): If True, don't actually insert into database, just show what would happen
    """
    session = StarbaseSession()
    
    try:
        # Get all ships that don't have a ship_accession yet
        ships_without_ssb = session.execute(
            text("""
            SELECT s.id, s.md5, s.sequence_length, a.accession_tag as ssa_accession
            FROM ships s
            LEFT JOIN ship_accessions sa ON sa.ship_id = s.id
            LEFT JOIN accessions a ON a.id = s.accession_id
            WHERE sa.id IS NULL
            ORDER BY s.id
            """)
        ).fetchall()
        
        total_ships = len(ships_without_ssb)
        logger.info(f"Found {total_ships} ships without SSB accessions")
        
        if total_ships == 0:
            logger.info("All ships already have SSB accessions!")
            return
        
        # Get the starting SSB number
        next_ssb_number = get_next_ssb_number(session)
        logger.info(f"Starting from SSB{next_ssb_number:07d}")
        
        # Generate SSB accessions
        inserted_count = 0
        for ship_id, md5, seq_length, ssa_accession in ships_without_ssb:
            ssb_tag = f"SSB{next_ssb_number:07d}"
            
            if dry_run:
                logger.info(
                    f"[DRY RUN] Would create: {ssb_tag} for ship_id={ship_id}, "
                    f"SSA={ssa_accession}, md5={md5[:8]}..., length={seq_length}"
                )
            else:
                # Insert the new ship_accession
                session.execute(
                    text("""
                    INSERT INTO ship_accessions (ship_accession_tag, version_tag, ship_id)
                    VALUES (:tag, :version, :ship_id)
                    """),
                    {"tag": ssb_tag, "version": 1, "ship_id": ship_id}
                )
                
                if inserted_count % 100 == 0:
                    logger.info(f"Created {inserted_count}/{total_ships} SSB accessions...")
                    session.commit()  # Commit in batches
                
                inserted_count += 1
            
            next_ssb_number += 1
        
        if not dry_run:
            session.commit()
            logger.info(f"Successfully created {inserted_count} SSB accessions")
        else:
            logger.info(f"[DRY RUN] Would have created {total_ships} SSB accessions")
        
        # Display summary statistics
        if not dry_run:
            total_ssb = session.execute(
                text("SELECT COUNT(*) FROM ship_accessions")
            ).fetchone()[0]
            logger.info(f"Total SSB accessions in database: {total_ssb}")
        
    except Exception as e:
        session.rollback()
        logger.error(f"Error generating SSB accessions: {str(e)}")
        raise
    finally:
        session.close()


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate SSB accessions for ships in the database"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually modifying the database"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Skip confirmation prompt"
    )
    
    args = parser.parse_args()
    
    if not args.force and not args.dry_run:
        response = input(
            "This will generate SSB accessions for all ships without one. "
            "Continue? (yes/no): "
        )
        if response.lower() not in ['yes', 'y']:
            logger.info("Aborted by user")
            return
    
    logger.info("Starting SSB accession generation...")
    generate_ssb_accessions(dry_run=args.dry_run)
    logger.info("Done!")


if __name__ == "__main__":
    main()

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
from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session
from src.utils.classification_utils import get_next_available_accession

logger = get_logger(__name__)


def generate_ssb_accessions(dry_run=False):
    """
    Generate SSB accessions for all ships that don't have one.

    Use for backfilling existing data. New ships get an SSB at insert time via
    ensure_ship_has_ssb() in the submission flow.

    Args:
        session: SQLAlchemy session.
        dry_run (bool): If True, don't actually insert into database, just show what would happen
    """
    with get_starbase_session() as session:
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

            # get_next_available_accession returns full tag string (e.g. "SSB0000002")
            next_ssb_tag = get_next_available_accession(session, "SSB")
            next_ssb_num = int(next_ssb_tag.replace("SSB", ""))
            logger.info(f"Starting from {next_ssb_tag}")

            # Generate SSB accessions
            inserted_count = 0
            for ship_id, md5, seq_length, ssa_accession in ships_without_ssb:
                ssb_tag = f"SSB{next_ssb_num:07d}"

                if dry_run:
                    logger.info(
                        f"[DRY RUN] Would create: {ssb_tag} for ship_id={ship_id}, "
                        f"SSA={ssa_accession}, md5={md5[:8] if md5 else 'N/A'}..., length={seq_length}"
                    )
                else:
                    # Insert the new ship_accession
                    session.execute(
                        text("""
                        INSERT INTO ship_accessions (ship_accession_tag, version_tag, ship_id)
                        VALUES (:tag, :version, :ship_id)
                        """),
                        {"tag": ssb_tag, "version": 1, "ship_id": ship_id},
                    )

                    if inserted_count % 100 == 0:
                        logger.info(
                            f"Created {inserted_count}/{total_ships} SSB accessions..."
                        )
                        session.commit()  # Commit in batches

                    inserted_count += 1

                next_ssb_num += 1

            if not dry_run:
                session.commit()
                logger.info(f"Successfully created {inserted_count} SSB accessions")
            else:
                logger.info(
                    f"[DRY RUN] Would have created {total_ships} SSB accessions"
                )

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


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate SSB accessions for ships in the database"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without actually modifying the database",
    )
    parser.add_argument("--force", action="store_true", help="Skip confirmation prompt")

    args = parser.parse_args()

    if not args.force and not args.dry_run:
        response = input(
            "This will generate SSB accessions for all ships without one. "
            "Continue? (yes/no): "
        )
        if response.lower() not in ["yes", "y"]:
            logger.info("Aborted by user")
            return

    logger.info("Starting SSB accession generation...")
    with get_starbase_session() as session:
        generate_ssb_accessions(session, dry_run=args.dry_run)
    logger.info("Done!")


if __name__ == "__main__":
    main()

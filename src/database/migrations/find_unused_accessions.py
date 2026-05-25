"""
Find accessions that are not referenced by any other table.

An accession is considered unused if its id does not appear in:
- ships.accession_id
- gff.accession_id
- joined_ships.accession_id
- starship_features.accession_id

Usage:
    python -m src.database.migrations.find_unused_accessions
    python -m src.database.migrations.find_unused_accessions --csv  # output as CSV
"""

import argparse
import csv
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from sqlalchemy import text
from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

UNUSED_ACCESSIONS_SQL = """
SELECT a.id, a.accession_tag, a.version_tag
FROM accessions a
LEFT JOIN ships s ON s.accession_id = a.id
LEFT JOIN gff g ON g.accession_id = a.id
LEFT JOIN joined_ships j ON j.accession_id = a.id
LEFT JOIN starship_features sf ON sf.accession_id = a.id
WHERE s.accession_id IS NULL
  AND g.accession_id IS NULL
  AND j.accession_id IS NULL
  AND sf.accession_id IS NULL
ORDER BY a.id
"""


def find_unused_accessions(csv_output=False):
    """Query the database for accessions not referenced by any other table."""
    with get_starbase_session() as session:
        rows = session.execute(text(UNUSED_ACCESSIONS_SQL)).fetchall()

    if csv_output:
        writer = csv.writer(sys.stdout)
        writer.writerow(("id", "accession_tag", "version_tag"))
        for row in rows:
            writer.writerow(row)
        return

    if not rows:
        logger.info("No unused accessions found.")
        return

    logger.info("Found %d unused accession(s):", len(rows))
    for row in rows:
        print(f"  id={row[0]}  accession_tag={row[1]}  version_tag={row[2]}")


def main():
    parser = argparse.ArgumentParser(
        description="Find accessions not referenced by ships, gff, joined_ships, or starship_features."
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Output as CSV (id, accession_tag, version_tag)",
    )
    args = parser.parse_args()
    find_unused_accessions(csv_output=args.csv)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Populate joined_ships.ship_id using ships via accession_id mapping.

Strategy:
- For joined_ships rows where accession_id is not NULL and ship_id is NULL,
  find a ships.id with ships.accession_id = joined_ships.accession_id.
- If exactly one ship matches, set ship_id to that ships.id.
- If multiple ships match for the accession, skip (ambiguous) and count.
- Report counts and a few examples.
"""

import sqlite3
import sys
from datetime import datetime

def populate(db_path: str, dry_run: bool = True):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    report = {
        'candidates': 0,
        'updated': 0,
        'ambiguous': 0,
        'no_match': 0,
        'examples_updated': [],
        'examples_ambiguous': [],
        'examples_no_match': []
    }

    try:
        # Find candidate joined_ships rows
        cur.execute(
            """
            SELECT js.id, js.accession_id
            FROM joined_ships js
            WHERE js.ship_id IS NULL AND js.accession_id IS NOT NULL
            """
        )
        rows = cur.fetchall()
        report['candidates'] = len(rows)

        for js_id, accession_id in rows:
            # Count ships for this accession
            cur.execute(
                "SELECT id FROM ships WHERE accession_id = ?",
                (accession_id,)
            )
            ship_rows = cur.fetchall()
            if not ship_rows:
                report['no_match'] += 1
                if len(report['examples_no_match']) < 5:
                    report['examples_no_match'].append({'joined_id': js_id, 'accession_id': accession_id})
                continue

            if len(ship_rows) > 1:
                report['ambiguous'] += 1
                if len(report['examples_ambiguous']) < 5:
                    report['examples_ambiguous'].append({
                        'joined_id': js_id,
                        'accession_id': accession_id,
                        'ship_ids': [r[0] for r in ship_rows]
                    })
                continue

            ship_id = ship_rows[0][0]
            if not dry_run:
                cur.execute(
                    "UPDATE joined_ships SET ship_id = ?, updated_at = ? WHERE id = ?",
                    (ship_id, datetime.now(), js_id)
                )
            report['updated'] += 1
            if len(report['examples_updated']) < 5:
                report['examples_updated'].append({'joined_id': js_id, 'accession_id': accession_id, 'ship_id': ship_id})

        if not dry_run:
            conn.commit()
        return report
    finally:
        conn.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python populate_ship_id.py <db_path> [--execute]")
        sys.exit(1)
    db_path = sys.argv[1]
    dry_run = '--execute' not in sys.argv[2:]
    rep = populate(db_path, dry_run=dry_run)
    mode = 'DRY RUN' if dry_run else 'APPLIED'
    print(f"\nPopulate ship_id: {mode}")
    print(f"Candidates: {rep['candidates']}")
    print(f"Updated: {rep['updated']}")
    print(f"Ambiguous: {rep['ambiguous']}")
    print(f"No match: {rep['no_match']}")
    if rep['examples_updated']:
        print("\nExamples updated:")
        for ex in rep['examples_updated']:
            print(f"  joined {ex['joined_id']}: accession {ex['accession_id']} -> ship {ex['ship_id']}")
    if rep['examples_ambiguous']:
        print("\nAmbiguous examples:")
        for ex in rep['examples_ambiguous']:
            print(f"  joined {ex['joined_id']}: accession {ex['accession_id']} ships={ex['ship_ids']}")
    if rep['examples_no_match']:
        print("\nNo match examples:")
        for ex in rep['examples_no_match']:
            print(f"  joined {ex['joined_id']}: accession {ex['accession_id']}")



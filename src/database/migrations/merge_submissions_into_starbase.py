"""
One-time data migration: copy rows from the standalone submissions.sqlite
into starbase.sqlite's submissions table (created by Alembic revision
m6n7o8p9q0r1). Not an Alembic migration -- it reads from a second SQLite
file, which Alembic isn't set up to do.

Run this AFTER applying Alembic migrations (so the target submissions
table exists), and only once per environment. Safe to re-run: skips rows
whose id already exists in the target table.

Usage:
    python src/database/migrations/merge_submissions_into_starbase.py \\
        --old-db /path/to/old/submissions.sqlite \\
        --new-db /path/to/starbase.sqlite
"""

import argparse
import sqlite3
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--old-db", required=True, help="Path to standalone submissions.sqlite"
    )
    parser.add_argument(
        "--new-db", required=True, help="Path to unified starbase.sqlite"
    )
    args = parser.parse_args()

    old_con = sqlite3.connect(args.old_db)
    new_con = sqlite3.connect(args.new_db)

    old_cur = old_con.execute("PRAGMA table_info(submissions)")
    old_cols = [row[1] for row in old_cur.fetchall()]

    new_cur = new_con.execute("PRAGMA table_info(submissions)")
    new_cols = {row[1] for row in new_cur.fetchall()}
    missing = set(old_cols) - new_cols
    if missing:
        raise SystemExit(
            f"Target submissions table is missing columns present in source: {missing}. "
            "Run Alembic migrations on --new-db first."
        )

    existing_ids = {row[0] for row in new_con.execute("SELECT id FROM submissions")}

    col_list = ", ".join(old_cols)
    placeholders = ", ".join("?" for _ in old_cols)
    insert_sql = f"INSERT INTO submissions ({col_list}) VALUES ({placeholders})"

    rows = old_con.execute(f"SELECT {col_list} FROM submissions").fetchall()
    to_insert = [row for row in rows if row[old_cols.index("id")] not in existing_ids]

    if not to_insert:
        print(
            f"Nothing to migrate: {len(rows)} source rows, all already present in target."
        )
        return

    new_con.executemany(insert_sql, to_insert)
    new_con.commit()
    print(
        f"Migrated {len(to_insert)} of {len(rows)} submission rows "
        f"({len(rows) - len(to_insert)} already present, skipped)."
    )

    old_con.close()
    new_con.close()


if __name__ == "__main__":
    main()

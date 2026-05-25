"""
Add classification columns to submissions table.

Run this migration to add columns for storing BLAST classification data:
- classification_source
- classification_family
- classification_navis
- classification_haplotype
- closest_match
- classification_confidence

Usage:
    python src/database/migrations/add_submission_classification_columns.py
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


def main():
    from sqlalchemy import text
    from src.database.sql_engine import get_submissions_session

    columns = [
        ("classification_source", "VARCHAR(50)"),
        ("classification_family", "VARCHAR(100)"),
        ("classification_navis", "VARCHAR(100)"),
        ("classification_haplotype", "VARCHAR(100)"),
        ("closest_match", "VARCHAR(50)"),
        ("classification_confidence", "VARCHAR(20)"),
    ]

    with get_submissions_session() as session:
        for col_name, col_type in columns:
            try:
                session.execute(
                    text(f"ALTER TABLE submissions ADD COLUMN {col_name} {col_type}")
                )
                print(f"Added column: {col_name}")
            except Exception as e:
                if "duplicate column name" in str(e).lower():
                    print(f"Column {col_name} already exists, skipping")
                else:
                    raise


if __name__ == "__main__":
    main()

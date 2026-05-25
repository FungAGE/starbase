"""
Apply Submission schema to the submissions database.
Creates the submissions table if it doesn't exist, or adds any missing columns
to match the current schema in src/database/models/schema.py.
Usage:
    python src/database/migrations/apply_submission_schema.py
"""

import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


def main():
    from sqlalchemy import text, inspect
    from src.database.sql_engine import get_submissions_session
    from src.config.database import submissions_engine

    # Full schema: (column_name, sqlite_type)
    columns_spec = [
        ("id", "INTEGER PRIMARY KEY AUTOINCREMENT"),
        ("seq_contents", "TEXT NOT NULL"),
        ("seq_filename", "VARCHAR(255) NOT NULL"),
        ("seq_date", "VARCHAR(50)"),
        ("anno_contents", "TEXT"),
        ("anno_filename", "VARCHAR(255)"),
        ("anno_date", "VARCHAR(50)"),
        ("uploader", "VARCHAR(255)"),
        ("evidence", "VARCHAR(100)"),
        ("genus", "VARCHAR(255)"),
        ("species", "VARCHAR(255)"),
        ("hostchr", "VARCHAR(255)"),
        ("shipstart", "INTEGER"),
        ("shipend", "INTEGER"),
        ("shipstrand", "VARCHAR(10)"),
        ("comment", "TEXT"),
        ("ship_accession_tag", "VARCHAR(50)"),
        ("accession_tag", "VARCHAR(50)"),
        ("needs_review", "BOOLEAN"),
        ("classification_source", "VARCHAR(50)"),
        ("classification_family", "VARCHAR(100)"),
        ("classification_navis", "VARCHAR(100)"),
        ("classification_haplotype", "VARCHAR(100)"),
        ("closest_match", "VARCHAR(50)"),
        ("classification_confidence", "VARCHAR(20)"),
    ]
    inspector = inspect(submissions_engine)
    with get_submissions_session() as session:
        if "submissions" not in inspector.get_table_names():
            # Create table from scratch
            col_defs = ", ".join(f"{name} {typ}" for name, typ in columns_spec)
            session.execute(text(f"CREATE TABLE submissions ({col_defs})"))
            session.commit()
            print("Created submissions table with all columns.")
        else:
            # Table exists: add any missing columns
            existing = {c["name"] for c in inspector.get_columns("submissions")}
            add_columns = [
                (name, typ)
                for name, typ in columns_spec
                if name != "id" and name not in existing
            ]
            for col_name, col_type in add_columns:
                try:
                    session.execute(
                        text(
                            f"ALTER TABLE submissions ADD COLUMN {col_name} {col_type}"
                        )
                    )
                    session.commit()
                    print(f"Added column: {col_name}")
                except Exception as e:
                    session.rollback()
                    if "duplicate column name" in str(e).lower():
                        print(f"Column {col_name} already exists, skipping")
                    else:
                        raise
            if not add_columns:
                print("Submissions table already has all columns.")


if __name__ == "__main__":
    main()

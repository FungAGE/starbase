"""
Update the type_ship flag on ships based on selection criteria.

A type_ship is "1" when it is the canonical representative for its accession_tag.
Run this script after data changes to refresh type_ship values.

Ensures the type_ship column exists on ships before updating (adds it if missing).

Usage:
    python src/database/migrations/update_type_ship.py
    python src/database/migrations/update_type_ship.py --seed 42  # reproducible tie-breaking
"""

import argparse
import sys
from pathlib import Path

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))


def main():
    parser = argparse.ArgumentParser(
        description="Update ships.type_ship based on accession_tag grouping and quality criteria."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for tie-breaking (optional, for reproducibility)",
    )
    args = parser.parse_args()

    from src.database.migrations import ensure_type_ship_column
    from src.database.type_ship_utils import update_type_ship_in_db

    ensure_type_ship_column()
    n = update_type_ship_in_db(seed=args.seed)
    print(f"Updated {n} ships as type_ship=1")


if __name__ == "__main__":
    main()

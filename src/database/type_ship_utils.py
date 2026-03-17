# coding: utf-8
"""
Utilities for computing and updating the type_ship flag.

A type_ship is "1" when it is the canonical representative for its accession_tag.
- If only one ship has that accession_tag, it gets type_ship="1".
- If multiple ships share the same accession_tag, we select one using these criteria
  (in order, until one remains; ties broken randomly):
  1. Most non-NULL columns among: accession_id, ship_id, ship_family_id, tax_id,
     genome_id, captain_id, ship_navis_id, ship_haplotype_id
  2. curated_status = "curated" preferred
  3. Fewest ship_quality_tags
  4. Earliest created_at
  5. Random tie-breaker
"""

import random
from typing import Set

import pandas as pd

from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

# Columns used for "completeness" score (more non-NULL = higher quality)
COMPLETENESS_COLUMNS = [
    "accession_id",
    "ship_id",
    "ship_family_id",
    "tax_id",
    "genome_id",
    "captain_id",
    "ship_navis_id",
    "ship_haplotype_id",
]


def _count_non_null(row: pd.Series, columns: list) -> int:
    """Count how many of the given columns are non-NULL in the row."""
    return sum(
        1
        for c in columns
        if c in row.index and pd.notna(row.get(c)) and row[c] is not None
    )


def _select_type_ship_for_group(group_df: pd.DataFrame, seed: int = 42):
    """
    Select the single best ship_id from a group of JoinedShips with the same accession_tag.

    Applies criteria in order until one remains.
    """
    if group_df.empty:
        return None
    if len(group_df) == 1:
        return int(group_df["ship_id"].iloc[0])

    df = group_df.copy()

    # 1. Most non-NULL completeness columns (higher = better)
    df["_completeness"] = df.apply(
        lambda r: _count_non_null(r, COMPLETENESS_COLUMNS), axis=1
    )
    max_completeness = df["_completeness"].max()
    df = df[df["_completeness"] == max_completeness]
    if len(df) == 1:
        return int(df["ship_id"].iloc[0])

    # 2. curated_status = "curated" preferred
    curated = df[df["curated_status"] == "curated"]
    if not curated.empty:
        df = curated
    if len(df) == 1:
        return int(df["ship_id"].iloc[0])

    # 3. Fewest ship_quality_tags (lower = better)
    df["_n_quality_tags"] = df["quality_tag_count"].fillna(0).astype(int)
    min_tags = df["_n_quality_tags"].min()
    df = df[df["_n_quality_tags"] == min_tags]
    if len(df) == 1:
        return int(df["ship_id"].iloc[0])

    # 4. Earliest created_at (skip if all NaN)
    if df["created_at"].notna().any():
        df = df.sort_values("created_at", na_position="last")
        earliest = df["created_at"].dropna().iloc[0]
        df = df[df["created_at"] == earliest]
        if len(df) == 1:
            return int(df["ship_id"].iloc[0])

    # 5. Random tie-breaker
    if df.empty:
        return None
    rng = random.Random(seed)
    idx = rng.randint(0, len(df) - 1)
    return int(df["ship_id"].iloc[idx])


def compute_type_ship_ship_ids(seed: int = 42) -> Set[int]:
    """
    Compute the set of ship_ids that should have type_ship="1".

    Returns:
        Set of ship IDs that are the type ship for their accession_tag.
    """
    query = """
    SELECT
        j.id AS joined_ship_id,
        j.ship_id,
        j.accession_id,
        j.ship_family_id,
        j.tax_id,
        j.genome_id,
        j.captain_id,
        j.ship_navis_id,
        j.ship_haplotype_id,
        j.curated_status,
        j.created_at,
        a.accession_tag,
        COALESCE(tag_counts.n_tags, 0) AS quality_tag_count
    FROM joined_ships j
    LEFT JOIN accessions a ON j.accession_id = a.id
    LEFT JOIN (
        SELECT joined_ship_id, COUNT(*) AS n_tags
        FROM ship_quality_tags
        GROUP BY joined_ship_id
    ) tag_counts ON tag_counts.joined_ship_id = j.id
    WHERE j.accession_id IS NOT NULL AND a.accession_tag IS NOT NULL
    """
    # Use one row per (accession_tag, ship_id) - if same ship has multiple joined_ships
    # rows for same accession, pick the best joined_ships row for that ship
    with get_starbase_session() as session:
        df = pd.read_sql_query(query, session.bind)

    if df.empty:
        return set()

    # For each (accession_tag, ship_id), keep the row with best metadata for that ship
    # (we'll pick one ship per accession_tag, so we need one row per ship per accession)
    def best_row_per_ship(g):
        # Prefer curated, fewer tags, earlier created_at
        # curated_status: "curated" < "uncurated" alphabetically, so ascending puts curated first
        g = g.sort_values(
            by=["curated_status", "quality_tag_count", "created_at"],
            ascending=[True, True, True],
            na_position="last",
        )
        return g.iloc[[0]]

    df = (
        df.groupby(["accession_tag", "ship_id"], group_keys=False)
        .apply(best_row_per_ship)
        .reset_index(drop=True)
    )

    type_ship_ids = set()
    for accession_tag, group in df.groupby("accession_tag"):
        ship_id = _select_type_ship_for_group(group, seed=seed)
        if ship_id is not None:
            type_ship_ids.add(ship_id)

    return type_ship_ids


def update_type_ship_in_db(seed: int = 42):
    """
    Update ships.type_ship in the database based on the selection criteria.

    Sets type_ship="1" for type ships, type_ship="0" for others in the same accession groups.
    Ships not in JoinedShips with an accession are left unchanged.

    Args:
        seed: Random seed for tie-breaking (optional, for reproducibility)
    Returns:
        Number of ships updated to type_ship="1".
    """
    from sqlalchemy import text

    from src.database.models.schema import Ships

    type_ship_ids = compute_type_ship_ship_ids(seed=seed)

    with get_starbase_session() as session:
        try:
            # Get all ship_ids that appear in joined_ships with an accession (candidates)
            result = session.execute(
                text(
                    "SELECT DISTINCT ship_id FROM joined_ships j "
                    "JOIN accessions a ON j.accession_id = a.id "
                    "WHERE j.accession_id IS NOT NULL AND a.accession_tag IS NOT NULL"
                )
            )
            candidate_ship_ids = [r[0] for r in result if r[0] is not None]

            # Reset candidates to "0"
            if candidate_ship_ids:
                session.query(Ships).filter(Ships.id.in_(candidate_ship_ids)).update(
                    {Ships.type_ship: "0"}, synchronize_session=False
                )
            # Set type ships to "1"
            if type_ship_ids:
                session.query(Ships).filter(Ships.id.in_(type_ship_ids)).update(
                    {Ships.type_ship: "1"}, synchronize_session=False
                )
            session.commit()
            logger.info(f"Updated type_ship: {len(type_ship_ids)} ships marked as type")
            return len(type_ship_ids)
        except Exception as e:
            session.rollback()
            logger.error(f"Error updating type_ship: {e}")
            raise


if __name__ == "__main__":
    update_type_ship_in_db()

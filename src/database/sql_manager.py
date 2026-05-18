"""
sql_manager — HTTP facade over the starbase compute backend.

All public functions match the original SQLAlchemy signatures exactly
(same args, same return types — DataFrames where applicable) so no
call sites need to change.

On the backend itself, routers import from sql_manager_impl directly.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from src.config import backend_client
from src.config.logging import get_logger

logger = get_logger(__name__)


# ── helpers ──────────────────────────────────────────────────────────────────

def _to_df(records: list[dict] | None) -> pd.DataFrame:
    if not records:
        return pd.DataFrame()
    return pd.DataFrame.from_records(records)


# ── data functions ───────────────────────────────────────────────────────────

def fetch_meta_data(curated: bool = False, accession_tags=None) -> pd.DataFrame:
    records = backend_client.fetch_meta_data(curated=curated, accession_tags=accession_tags)
    return _to_df(records)


def fetch_paper_data() -> pd.DataFrame:
    return _to_df(backend_client.fetch_paper_data())


def dereplicate_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """Pure-Python deduplicate by md5/rev_comp_md5; no network call needed."""
    seen_sequences: set[str] = set()
    indices_to_keep = []

    for idx, row in df.iterrows():
        md5_val = row.get("md5", "")
        rev_comp_md5_val = row.get("rev_comp_md5", "")

        if not md5_val or not rev_comp_md5_val:
            indices_to_keep.append(idx)
            continue
        if md5_val not in seen_sequences and rev_comp_md5_val not in seen_sequences:
            indices_to_keep.append(idx)
            seen_sequences.add(md5_val)
            seen_sequences.add(rev_comp_md5_val)

    return df.loc[indices_to_keep]


def fetch_ships(
    accession_tags=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> pd.DataFrame:
    records = backend_client.fetch_ships(
        accession_tags=accession_tags,
        curated=curated,
        dereplicate=dereplicate,
        with_sequence=with_sequence,
    )
    return _to_df(records)


def fetch_ship_table(
    curated: bool = True,
    with_sequence: bool = False,
    with_gff_entries: bool = False,
) -> pd.DataFrame:
    records = backend_client.fetch_ship_table(
        curated=curated,
        with_sequence=with_sequence,
        with_gff_entries=with_gff_entries,
    )
    return _to_df(records)


def fetch_accession_ship(accession_tag: str) -> dict[str, pd.DataFrame | None]:
    result = backend_client.fetch_accession_ship(accession_tag)
    seq = result.get("sequence")
    gff = result.get("gff")
    return {
        "sequence": _to_df(seq) if seq is not None else None,
        "gff": _to_df(gff) if gff is not None else None,
    }


def fetch_captains(
    accession_tags=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> pd.DataFrame:
    records = backend_client.fetch_captains(
        accession_tags=accession_tags,
        curated=curated,
        dereplicate=dereplicate,
        with_sequence=with_sequence,
    )
    return _to_df(records)


def fetch_captain_tree() -> str:
    return backend_client.fetch_captain_tree()


def fetch_sf_data() -> pd.DataFrame:
    return _to_df(backend_client.fetch_sf_data())


def get_database_version() -> str:
    return backend_client.get_database_version()


def set_database_version(
    semantic_version: str, description: str = "", created_by: str = "manual"
) -> bool:
    return backend_client.set_database_version(
        semantic_version, description=description, created_by=created_by
    )


def get_alembic_schema_version() -> str:
    return backend_client.get_alembic_schema_version()


def get_database_stats() -> dict[str, Any]:
    return backend_client.get_database_stats()


def add_quality_tag(
    joined_ship_id: int,
    tag_type: str,
    tag_value: str | None = None,
    created_by: str = "auto",
) -> int:
    return backend_client.add_quality_tag(
        joined_ship_id, tag_type, tag_value=tag_value, created_by=created_by
    )


def remove_quality_tag(joined_ship_id: int, tag_type: str) -> bool:
    return backend_client.remove_quality_tag(joined_ship_id, tag_type)


def get_quality_tags(joined_ship_id: int) -> list[dict[str, Any]]:
    return backend_client.get_quality_tags(joined_ship_id)

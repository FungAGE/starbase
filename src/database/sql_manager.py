"""
sql_manager — dual-mode data access facade.

When BACKEND_API_URL is set: HTTP to compute backend (with optional response cache).
When unset: direct SQLAlchemy via sql_manager_impl (local DB debug / monolith).

Backend routers import sql_manager_impl directly.
"""

from __future__ import annotations

from typing import Any, Callable, TypeVar

import pandas as pd

from src.config import backend_client
from src.config.cache import smart_cache
from src.config.logging import get_logger

logger = get_logger(__name__)

T = TypeVar("T")

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.sql_manager_impl as impl

        _impl = impl
    return _impl


def _to_df(records: list[dict] | None) -> pd.DataFrame:
    if not records:
        return pd.DataFrame()
    return pd.DataFrame.from_records(records)


def _resolve_accession_tags(accession_tags=None, accessions=None):
    """Accept both ``accession_tags`` (API/split) and ``accessions`` (legacy call sites)."""
    return accession_tags if accession_tags is not None else accessions


def _via_backend(
    http_fn: Callable[..., T], local_fn: Callable[..., T], *args, **kwargs
) -> T:
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def _filter_meta_df(
    df: pd.DataFrame,
    curated: bool,
    accessions,
) -> pd.DataFrame:
    impl = _get_impl()
    accession_mode = impl._get_accession_mode(accessions)
    filtered_df = df.copy()

    if curated:
        filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]
    if accessions:
        formatted_values = [str(tag).strip("'\"") for tag in accessions]
        if accession_mode == "SSA":
            filtered_df = filtered_df[
                filtered_df["accession_tag"].isin(formatted_values)
            ]
        elif accession_mode == "SSB":
            filtered_df = filtered_df[
                filtered_df["ship_accession_tag"].isin(formatted_values)
            ]
        elif accession_mode is None:
            mask = filtered_df["accession_tag"].isin(formatted_values) | filtered_df[
                "ship_accession_tag"
            ].isin(formatted_values)
            filtered_df = filtered_df[mask]
        else:
            raise ValueError(f"Invalid accession mode: {accession_mode}")
    return filtered_df


@smart_cache(timeout=None, key_prefix="http_fetch_meta_data:full:v2")
def _http_fetch_meta_data_full() -> pd.DataFrame:
    records = backend_client.fetch_meta_data(curated=False, accession_tags=None)
    return _to_df(records)


@smart_cache(timeout=None, key_prefix="http_fetch_paper_data")
def _http_fetch_paper_data() -> pd.DataFrame:
    return _to_df(backend_client.fetch_paper_data())


@smart_cache(timeout=None, key_prefix="http_fetch_ship_table:full")
def _http_fetch_ship_table_full() -> pd.DataFrame:
    records = backend_client.fetch_ship_table(
        curated=False, with_sequence=False, with_gff_entries=False
    )
    return _to_df(records)


@smart_cache(timeout=None, key_prefix="http_get_database_stats")
def _http_get_database_stats() -> dict[str, Any]:
    return backend_client.get_database_stats()


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


def fetch_meta_data(
    curated: bool = False, accession_tags=None, accessions=None
) -> pd.DataFrame:
    tags = _resolve_accession_tags(accession_tags, accessions)
    normalized = _get_impl()._normalize_accession_tags(tags, None)

    def _local():
        return _get_impl().fetch_meta_data(
            curated=curated, accession_tags=tags, accessions=accessions
        )

    def _http():
        full_df = _http_fetch_meta_data_full()
        return _filter_meta_df(full_df, curated, normalized)

    return _via_backend(_http, _local)


def fetch_paper_data() -> pd.DataFrame:
    return _via_backend(_http_fetch_paper_data, _get_impl().fetch_paper_data)


def fetch_ships(
    accession_tags=None,
    accessions=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> pd.DataFrame:
    tags = _resolve_accession_tags(accession_tags, accessions)

    def _local():
        return _get_impl().fetch_ships(
            accession_tags=tags,
            accessions=accessions,
            curated=curated,
            dereplicate=dereplicate,
            with_sequence=with_sequence,
        )

    def _http():
        records = backend_client.fetch_ships(
            accession_tags=tags,
            curated=curated,
            dereplicate=dereplicate,
            with_sequence=with_sequence,
        )
        return _to_df(records)

    return _via_backend(_http, _local)


def fetch_ship_table(
    curated: bool = True,
    with_sequence: bool = False,
    with_gff_entries: bool = False,
) -> pd.DataFrame:
    def _local():
        return _get_impl().fetch_ship_table(
            curated=curated,
            with_sequence=with_sequence,
            with_gff_entries=with_gff_entries,
        )

    def _http():
        if with_gff_entries:
            records = backend_client.fetch_ship_table(
                curated=curated,
                with_sequence=with_sequence,
                with_gff_entries=True,
            )
            return _to_df(records)

        filtered_df = _http_fetch_ship_table_full().copy()
        if with_sequence:
            filtered_df = filtered_df[filtered_df["ship_id"].notna()]
        if curated:
            filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]
        if "familyName" in filtered_df.columns:
            filtered_df = filtered_df.sort_values(by="familyName")
        return filtered_df

    return _via_backend(_http, _local)


def fetch_accession_ship(accession_tag: str) -> dict[str, pd.DataFrame | None]:
    def _local():
        return _get_impl().fetch_accession_ship(accession_tag)

    def _http():
        result = backend_client.fetch_accession_ship(accession_tag)
        seq = result.get("sequence")
        gff = result.get("gff")
        return {
            "sequence": _to_df(seq) if seq is not None else None,
            "gff": _to_df(gff) if gff is not None else None,
        }

    return _via_backend(_http, _local)


def fetch_captains(
    accession_tags=None,
    accessions=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> pd.DataFrame:
    tags = _resolve_accession_tags(accession_tags, accessions)

    def _local():
        return _get_impl().fetch_captains(
            accession_tags=tags,
            accessions=accessions,
            curated=curated,
            dereplicate=dereplicate,
            with_sequence=with_sequence,
        )

    def _http():
        records = backend_client.fetch_captains(
            accession_tags=tags,
            curated=curated,
            dereplicate=dereplicate,
            with_sequence=with_sequence,
        )
        return _to_df(records)

    return _via_backend(_http, _local)


def fetch_captain_tree() -> str:
    return _via_backend(
        backend_client.fetch_captain_tree, _get_impl().fetch_captain_tree
    )


def fetch_sf_data() -> pd.DataFrame:
    def _http():
        return _to_df(backend_client.fetch_sf_data())

    return _via_backend(_http, _get_impl().fetch_sf_data)


def get_database_version() -> str:
    return _via_backend(
        backend_client.get_database_version, _get_impl().get_database_version
    )


def set_database_version(
    semantic_version: str, description: str = "", created_by: str = "manual"
) -> bool:
    def _http():
        return backend_client.set_database_version(
            semantic_version, description=description, created_by=created_by
        )

    def _local():
        return _get_impl().set_database_version(
            semantic_version, description=description, created_by=created_by
        )

    return _via_backend(_http, _local)


def get_alembic_schema_version() -> str:
    return _via_backend(
        backend_client.get_alembic_schema_version,
        _get_impl().get_alembic_schema_version,
    )


def get_database_stats() -> dict[str, Any]:
    return _via_backend(_http_get_database_stats, _get_impl().get_database_stats)


def add_quality_tag(
    joined_ship_id: int,
    tag_type: str,
    tag_value: str | None = None,
    created_by: str = "auto",
) -> int:
    def _http():
        return backend_client.add_quality_tag(
            joined_ship_id, tag_type, tag_value=tag_value, created_by=created_by
        )

    def _local():
        return _get_impl().add_quality_tag(
            joined_ship_id, tag_type, tag_value=tag_value, created_by=created_by
        )

    return _via_backend(_http, _local)


def remove_quality_tag(joined_ship_id: int, tag_type: str) -> bool:
    def _http():
        return backend_client.remove_quality_tag(joined_ship_id, tag_type)

    def _local():
        return _get_impl().remove_quality_tag(joined_ship_id, tag_type)

    return _via_backend(_http, _local)


def get_quality_tags(joined_ship_id: int) -> list[dict[str, Any]]:
    def _http():
        return backend_client.get_quality_tags(joined_ship_id)

    def _local():
        return _get_impl().get_quality_tags(joined_ship_id)

    return _via_backend(_http, _local)

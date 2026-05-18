"""REST routers that delegate to src.database.sql_manager (runs on backend only)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from fastapi import APIRouter
from pydantic import BaseModel

from backend.dependencies import CacheContext, RequireApiKey
from src.database import sql_manager

router = APIRouter(prefix="/api/data", tags=["data"], dependencies=[RequireApiKey, CacheContext])


def _df_to_records(df: pd.DataFrame | None) -> list[dict[str, Any]] | None:
    if df is None:
        return None
    if df.empty:
        return []
    out = df.replace({np.nan: None})
    return out.to_dict(orient="records")


class FetchMetaBody(BaseModel):
    curated: bool = False
    accession_tags: str | list[str] | None = None


@router.post("/meta")
def fetch_meta(body: FetchMetaBody) -> list[dict[str, Any]]:
    df = sql_manager.fetch_meta_data(
        curated=body.curated, accession_tags=body.accession_tags
    )
    return _df_to_records(df) or []


@router.get("/papers")
def fetch_papers() -> list[dict[str, Any]]:
    df = sql_manager.fetch_paper_data()
    return _df_to_records(df) or []


class FetchShipsBody(BaseModel):
    accession_tags: list[str] | None = None
    curated: bool = False
    dereplicate: bool = True
    with_sequence: bool = False


@router.post("/ships")
def fetch_ships(body: FetchShipsBody) -> list[dict[str, Any]]:
    df = sql_manager.fetch_ships(
        accession_tags=body.accession_tags,
        curated=body.curated,
        dereplicate=body.dereplicate,
        with_sequence=body.with_sequence,
    )
    return _df_to_records(df) or []


class FetchShipTableBody(BaseModel):
    curated: bool = True
    with_sequence: bool = False
    with_gff_entries: bool = False


@router.post("/ship-table")
def fetch_ship_table(body: FetchShipTableBody) -> list[dict[str, Any]]:
    df = sql_manager.fetch_ship_table(
        curated=body.curated,
        with_sequence=body.with_sequence,
        with_gff_entries=body.with_gff_entries,
    )
    return _df_to_records(df) or []


@router.get("/accession/{accession_tag}/ship")
def fetch_accession_ship(accession_tag: str) -> dict[str, Any]:
    result = sql_manager.fetch_accession_ship(accession_tag)
    seq = result.get("sequence")
    gff = result.get("gff")
    return {
        "sequence": _df_to_records(seq) if seq is not None else None,
        "gff": _df_to_records(gff) if gff is not None else None,
    }


class FetchCaptainsBody(BaseModel):
    accession_tags: list[str] | None = None
    curated: bool = False
    dereplicate: bool = True
    with_sequence: bool = False


@router.post("/captains")
def fetch_captains(body: FetchCaptainsBody) -> list[dict[str, Any]]:
    df = sql_manager.fetch_captains(
        accession_tags=body.accession_tags,
        curated=body.curated,
        dereplicate=body.dereplicate,
        with_sequence=body.with_sequence,
    )
    return _df_to_records(df) or []


@router.get("/captain-tree")
def fetch_captain_tree() -> dict[str, str]:
    return {"newick": sql_manager.fetch_captain_tree()}


@router.get("/sf-data")
def fetch_sf_data() -> list[dict[str, Any]]:
    df = sql_manager.fetch_sf_data()
    return _df_to_records(df) or []


@router.get("/database-version")
def database_version() -> dict[str, str]:
    return {"semantic_version": sql_manager.get_database_version()}


class SetDatabaseVersionBody(BaseModel):
    semantic_version: str
    description: str = ""
    created_by: str = "manual"


@router.post("/database-version")
def set_database_version(body: SetDatabaseVersionBody) -> dict[str, bool]:
    sql_manager.set_database_version(
        body.semantic_version,
        description=body.description,
        created_by=body.created_by,
    )
    return {"ok": True}


@router.get("/alembic-version")
def alembic_version() -> dict[str, str]:
    return {"revision": sql_manager.get_alembic_schema_version()}


@router.get("/stats")
def database_stats() -> dict[str, Any]:
    return sql_manager.get_database_stats()


class AddQualityTagBody(BaseModel):
    joined_ship_id: int
    tag_type: str
    tag_value: str | None = None
    created_by: str = "api"


@router.post("/quality-tags")
def add_quality_tag(body: AddQualityTagBody) -> dict[str, int]:
    tag_id = sql_manager.add_quality_tag(
        body.joined_ship_id,
        body.tag_type,
        tag_value=body.tag_value,
        created_by=body.created_by,
    )
    return {"id": tag_id}


class RemoveQualityTagBody(BaseModel):
    joined_ship_id: int
    tag_type: str


@router.post("/quality-tags/remove")
def remove_quality_tag(body: RemoveQualityTagBody) -> dict[str, bool]:
    ok = sql_manager.remove_quality_tag(body.joined_ship_id, body.tag_type)
    return {"removed": ok}


@router.get("/joined-ships/{joined_ship_id}/quality-tags")
def get_quality_tags(joined_ship_id: int) -> list[dict[str, Any]]:
    tags = sql_manager.get_quality_tags(joined_ship_id)
    out: list[dict[str, Any]] = []
    for t in tags:
        row = dict(t)
        ca = row.get("created_at")
        if ca is not None and hasattr(ca, "isoformat"):
            row["created_at"] = ca.isoformat()
        out.append(row)
    return out

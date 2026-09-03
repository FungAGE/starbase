"""REST router for admin grid read/write operations (runs on backend only)."""

from __future__ import annotations

import glob
import os
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from backend.dependencies import RequireApiKey
from src.database import admin_jobs_manager_impl as admin_jobs_manager
from src.database import admin_manager_impl as admin_manager
from src.database import admin_submissions_manager_impl as admin_submissions_manager

router = APIRouter(prefix="/api/admin", tags=["admin"], dependencies=[RequireApiKey])


def _df_to_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    if df.empty:
        return []
    out = df.replace({np.nan: None})
    return out.to_dict(orient="records")


@router.get("/table/{table_key}")
def fetch_admin_table(table_key: str) -> list[dict[str, Any]]:
    try:
        df = admin_manager.fetch_admin_table(table_key)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    return _df_to_records(df)


class AdminInsertBody(BaseModel):
    table_key: str
    col_values: Dict[str, Any]


@router.post("/insert")
def admin_insert(body: AdminInsertBody) -> dict[str, Any]:
    success, error, new_id = admin_manager.do_insert(body.table_key, body.col_values)
    return {"success": success, "error": error, "new_id": new_id}


class AdminUpdateBody(BaseModel):
    table_key: str
    row_id: Union[int, str]
    col_id: str
    new_value: Optional[Any] = None


@router.post("/update")
def admin_update(body: AdminUpdateBody) -> dict[str, Any]:
    success, error = admin_manager.do_update(
        body.table_key, body.row_id, body.col_id, body.new_value
    )
    return {"success": success, "error": error}


@router.post("/jobs/{job_key}/run")
def run_admin_job(job_key: str) -> dict[str, Any]:
    try:
        return admin_jobs_manager.run_job(job_key)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


class ProcessSubmissionsBody(BaseModel):
    sub_ids: list[int]


@router.post("/submissions/process")
def process_submissions(body: ProcessSubmissionsBody) -> list[dict[str, Any]]:
    return admin_submissions_manager.process_submissions(body.sub_ids)


class PromoteSubmissionBody(BaseModel):
    sub_id: int


@router.post("/submissions/promote")
def promote_submission(body: PromoteSubmissionBody) -> dict[str, Any]:
    success, accession, error = admin_submissions_manager.promote_submission(
        body.sub_id
    )
    return {"success": success, "accession": accession, "error": error}


def _blastdb_status() -> dict[str, Any]:
    from src.config.settings import BLAST_DB_PATHS
    from src.database.blastdb import blast_db_exists, sourmash_sig_exists

    paths = {
        "ships_all": BLAST_DB_PATHS["ship"]["all"]["nucl"],
        "ships_curated": BLAST_DB_PATHS["ship"]["curated"]["nucl"],
        "captains": BLAST_DB_PATHS["gene"]["tyr"]["prot"],
    }
    status = {}
    for name, path in paths.items():
        prefix, _ = os.path.splitext(path)
        index_files = sorted(glob.glob(prefix + ".*"))
        newest = max((os.path.getmtime(p) for p in index_files), default=None)
        status[name] = {
            "path": path,
            "built": blast_db_exists(path),
            "signatures": sourmash_sig_exists(path),
            "updated_at": newest,
            "file_count": len(index_files),
        }
    return status


@router.get("/blastdb/status")
def blastdb_status() -> dict[str, Any]:
    """Build state of the managed BLAST databases (ships all/curated, captains)."""
    return _blastdb_status()


@router.post("/blastdb/rebuild")
def blastdb_rebuild() -> dict[str, Any]:
    """Dispatch a full BLAST DB rebuild to the worker (runs in the background)."""
    from backend.tasks.blastdb import dispatch_blastdb_rebuild

    return dispatch_blastdb_rebuild()

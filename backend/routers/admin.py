"""REST router for admin grid read/write operations (runs on backend only)."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from backend.dependencies import RequireApiKey
from src.database import admin_manager_impl as admin_manager

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

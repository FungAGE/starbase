"""REST router for the curation/annotation-review workflow (runs on backend only)."""

from __future__ import annotations

from typing import Any, Dict, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from backend.dependencies import RequireApiKey
from src.database import curation_manager_impl as curation_manager

router = APIRouter(
    prefix="/api/curation", tags=["curation"], dependencies=[RequireApiKey]
)


class QueueBody(BaseModel):
    flag: Optional[int] = None
    assigned_to: Optional[str] = None
    limit: int = 50


@router.post("/queue")
def fetch_annotation_queue(body: QueueBody) -> list[dict[str, Any]]:
    return curation_manager.fetch_annotation_queue(
        flag=body.flag, assigned_to=body.assigned_to, limit=body.limit
    )


@router.get("/annotations/{annotation_id}")
def fetch_annotation(annotation_id: int) -> dict[str, Any]:
    try:
        return curation_manager.fetch_annotation(annotation_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


class UpdateAnnotationBody(BaseModel):
    changes: Dict[str, Any]
    changed_by: str = "curator"


@router.post("/annotations/{annotation_id}/update")
def update_annotation(annotation_id: int, body: UpdateAnnotationBody) -> dict[str, Any]:
    try:
        return curation_manager.update_annotation(
            annotation_id, body.changes, body.changed_by
        )
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))

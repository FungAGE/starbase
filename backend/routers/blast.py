"""BLAST / HMMER endpoints — compute runs on backend via src.tasks."""

from __future__ import annotations

from typing import Any, Literal, Optional

from fastapi import APIRouter
from pydantic import BaseModel, Field

from backend.dependencies import RequireApiKey
from src.config.logging import get_logger
from src.tasks import run_blast_search_task, run_hmmer_search_task

logger = get_logger(__name__)


def _run_blast_search(
    query_header: str,
    query_seq: str,
    query_type: str,
    *,
    eval_threshold: float,
    curated: Optional[bool],
):
    fn = run_blast_search_task
    if hasattr(fn, "apply"):
        return fn.apply(
            args=(query_header, query_seq, query_type),
            kwargs={"eval_threshold": eval_threshold, "curated": curated},
        ).get(timeout=360)
    return fn(
        query_header,
        query_seq,
        query_type,
        eval_threshold=eval_threshold,
        curated=curated,
    )


def _run_hmmer_search(
    query_header: str,
    query_seq: str,
    query_type: str,
    *,
    eval_threshold: float,
):
    fn = run_hmmer_search_task
    if hasattr(fn, "apply"):
        return fn.apply(
            args=(query_header, query_seq, query_type),
            kwargs={"eval_threshold": eval_threshold},
        ).get(timeout=360)
    return fn(query_header, query_seq, query_type, eval_threshold=eval_threshold)


router = APIRouter(
    prefix="/api/blast",
    tags=["blast"],
    dependencies=[RequireApiKey],
)


@router.post("/blast-submit")
def blast_submit() -> dict[str, bool]:
    """Compatibility shim with Flask route; rate limiting can move here later."""
    return {"allowed": True}


class BlastSearchBody(BaseModel):
    query_header: str = Field(..., description="FASTA header line (without >)")
    query_seq: str = Field(..., description="Raw sequence string")
    query_type: Literal["nucl", "prot"] = "nucl"
    eval_threshold: float = 0.01
    curated: Optional[bool] = None


@router.post("/search")
def blast_search(body: BlastSearchBody) -> dict[str, Any]:
    result = _run_blast_search(
        body.query_header,
        body.query_seq,
        body.query_type,
        eval_threshold=body.eval_threshold,
        curated=body.curated,
    )
    if result is None:
        logger.warning("BLAST search returned no results")
        return {"ok": False, "error": "BLAST failed or no output", "result": None}
    return {"ok": True, "result": result}


class HmmerSearchBody(BaseModel):
    query_header: str
    query_seq: str
    query_type: Literal["nucl", "prot"] = "nucl"
    eval_threshold: float = 0.01


@router.post("/hmmer")
def hmmer_search(body: HmmerSearchBody) -> dict[str, Any]:
    result = _run_hmmer_search(
        body.query_header,
        body.query_seq,
        body.query_type,
        eval_threshold=body.eval_threshold,
    )
    if result is None:
        logger.warning("HMMER search returned no results")
        return {"ok": False, "error": "HMMER failed or no output", "result": None}
    return {"ok": True, "result": result}

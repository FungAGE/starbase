"""Classification workflow endpoint — heavy compute runs on backend.

The workflow (exact → contained → similar → family → navis → haplotype) shells
out to minimap2/sourmash/HMMER/mmseqs, which only exist in the backend env.
Frontend pages call this synchronously through backend_client, mirroring the
blast/hmmer pattern.
"""

from __future__ import annotations

from typing import Any, Optional

from fastapi import APIRouter
from pydantic import BaseModel

from backend.dependencies import RequireApiKey
from backend.task_runner import run_task_sync
from src.config.logging import get_logger
from src.tasks import run_classification_workflow_task

logger = get_logger(__name__)

# Sequential subprocess stages against the full ship DB can take a while on
# large genomes. The frontend client mirrors this timeout.
WORKFLOW_TIMEOUT = 1800


class ClassificationWorkflowBody(BaseModel):
    workflow_state: dict[str, Any]
    blast_data: Optional[dict[str, Any]] = None
    classification_data: Optional[dict[str, Any]] = None
    meta_dict: Optional[list[dict[str, Any]]] = None


router = APIRouter(
    prefix="/api/classification",
    tags=["classification"],
    dependencies=[RequireApiKey],
)


@router.post("/workflow")
def classification_workflow(body: ClassificationWorkflowBody) -> dict[str, Any]:
    try:
        result = run_task_sync(
            run_classification_workflow_task,
            args=(body.workflow_state,),
            kwargs={
                "blast_data": body.blast_data,
                "classification_data": body.classification_data,
                "meta_dict": body.meta_dict,
            },
            timeout=WORKFLOW_TIMEOUT,
        )
    except Exception as exc:
        logger.error(f"Classification workflow failed: {exc}")
        logger.exception("Full traceback:")
        return {"ok": False, "error": str(exc), "result": None}

    if result is None:
        logger.warning("Classification workflow returned no result")
        return {"ok": False, "error": "Workflow returned no result", "result": None}
    return {"ok": True, "result": result}

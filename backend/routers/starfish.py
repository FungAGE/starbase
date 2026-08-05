"""REST router for Starfish pipeline run management (runs on backend only).

Slice 4a: create/list/get only. Start/cancel/rerun/resume land in slice 4b
alongside the actual nextflow subprocess launch.
"""

from __future__ import annotations

from typing import Any, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from backend.dependencies import RequireApiKey
from src.database import starfish_manager_impl as starfish_manager

router = APIRouter(
    prefix="/api/starfish", tags=["starfish"], dependencies=[RequireApiKey]
)


class ListRunsBody(BaseModel):
    status: Optional[str] = None
    limit: int = 50


@router.post("/runs/list")
def list_runs(body: ListRunsBody) -> list[dict[str, Any]]:
    return starfish_manager.list_runs(status=body.status, limit=body.limit)


@router.get("/runs/{run_id}")
def get_run(run_id: int) -> dict[str, Any]:
    try:
        return starfish_manager.get_run(run_id)
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


class GenomeInput(BaseModel):
    genome_id: str
    fna_path: str
    gff3_path: str
    tax_id: Optional[str] = None
    emapper_path: Optional[str] = None
    cds_path: Optional[str] = None
    faa_path: Optional[str] = None


class CreateRunBody(BaseModel):
    run_name: str
    genomes: List[GenomeInput]
    description: str = ""
    created_by: str = "curator"
    model: str = "tyr"
    threads: int = 20
    missing: int = 1
    maxcopy: int = 5
    pid_threshold: int = 90
    hsp: int = 1000
    flank: int = 6
    neighbourhood: int = 10000


@router.post("/runs")
def create_run(body: CreateRunBody) -> dict[str, Any]:
    genomes = [g.model_dump() for g in body.genomes]
    try:
        return starfish_manager.create_run(
            body.run_name,
            genomes,
            description=body.description,
            created_by=body.created_by,
            model=body.model,
            threads=body.threads,
            missing=body.missing,
            maxcopy=body.maxcopy,
            pid_threshold=body.pid_threshold,
            hsp=body.hsp,
            flank=body.flank,
            neighbourhood=body.neighbourhood,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@router.post("/runs/{run_id}/start")
def start_run(run_id: int) -> dict[str, Any]:
    try:
        return starfish_manager.start_run(run_id)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@router.post("/runs/{run_id}/cancel")
def cancel_run(run_id: int) -> dict[str, Any]:
    try:
        return starfish_manager.cancel_run(run_id)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@router.post("/runs/{run_id}/rerun")
def rerun_run(run_id: int) -> dict[str, Any]:
    try:
        return starfish_manager.rerun_run(run_id)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@router.post("/runs/{run_id}/resume")
def resume_run(run_id: int) -> dict[str, Any]:
    try:
        return starfish_manager.resume_run(run_id)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


class ImportElementBody(BaseModel):
    uploader: str = "starfish-pipeline"


@router.post("/elements/{element_id}/import")
def import_element(element_id: int, body: ImportElementBody) -> dict[str, Any]:
    try:
        return starfish_manager.import_element_to_submission(
            element_id, uploader=body.uploader
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))


@router.get("/runs/{run_id}/log")
def get_run_log(run_id: int) -> dict[str, str]:
    try:
        return {"log": starfish_manager.get_run_log(run_id)}
    except ValueError as exc:
        raise HTTPException(status_code=404, detail=str(exc))

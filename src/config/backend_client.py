"""
HTTP client for the starbase compute backend (FastAPI on institute machine).

All calls are synchronous (httpx sync client) so they drop in where
sql_manager calls currently live — no async plumbing needed in Dash callbacks.

Configuration (env vars):
    BACKEND_API_URL   - base URL, e.g. http://100.x.y.z:8001 or http://backend:8001
                        When unset, BackendClient.is_configured() returns False.
    BACKEND_API_KEY   - shared secret sent as X-API-Key header (optional but recommended)
    BACKEND_TIMEOUT   - per-request timeout in seconds (default 60)
    BACKEND_BLAST_TIMEOUT - timeout for BLAST/HMMER calls in seconds (default 360)
"""

from __future__ import annotations

import os
from typing import Any

import httpx

from src.config.logging import get_logger

logger = get_logger(__name__)

_BACKEND_API_URL = os.getenv("BACKEND_API_URL", "").rstrip("/")
_BACKEND_API_KEY = os.getenv("BACKEND_API_KEY", "")
_DEFAULT_TIMEOUT = float(os.getenv("BACKEND_TIMEOUT", "60"))
_BLAST_TIMEOUT = float(os.getenv("BACKEND_BLAST_TIMEOUT", "360"))
_CLASSIFICATION_TIMEOUT = float(os.getenv("BACKEND_CLASSIFICATION_TIMEOUT", "1800"))

# Set by start-script.sh when TS_AUTHKEY is present
_TAILSCALE_PROXY = os.getenv("TAILSCALE_PROXY", "") or None

_MAX_RETRIES = 3
_RETRY_STATUS_CODES = {429, 500, 502, 503, 504}


def _build_headers() -> dict[str, str]:
    headers: dict[str, str] = {"Content-Type": "application/json"}
    if _BACKEND_API_KEY:
        headers["X-API-Key"] = _BACKEND_API_KEY
    return headers


def _client(timeout: float = _DEFAULT_TIMEOUT) -> httpx.Client:
    """Return a configured httpx.Client. Caller is responsible for .close() / context manager."""
    return httpx.Client(
        base_url=_BACKEND_API_URL,
        headers=_build_headers(),
        timeout=timeout,
        follow_redirects=True,
        proxy=_TAILSCALE_PROXY,
    )


def is_configured() -> bool:
    """Return True when BACKEND_API_URL is set — i.e. backend split is active."""
    return bool(_BACKEND_API_URL)


def health_check() -> dict[str, str]:
    """Verify reachability and API-key authentication against the backend."""
    return _request("GET", "/api/health")


def _request(
    method: str,
    path: str,
    *,
    json: Any = None,
    timeout: float = _DEFAULT_TIMEOUT,
) -> Any:
    """
    Execute a request with simple linear retry on transient errors.

    Returns parsed JSON payload on success.
    Raises httpx.HTTPStatusError on non-retriable 4xx, RuntimeError on exhausted retries.
    """
    if not is_configured():
        raise RuntimeError(
            "BACKEND_API_URL is not set; use sql_manager local impl or set BACKEND_API_URL"
        )

    last_exc: Exception | None = None
    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            with _client(timeout) as client:
                response = client.request(method, path, json=json)
                if (
                    response.status_code in _RETRY_STATUS_CODES
                    and attempt < _MAX_RETRIES
                ):
                    logger.warning(
                        "Backend %s %s → %s, retry %d/%d",
                        method,
                        path,
                        response.status_code,
                        attempt,
                        _MAX_RETRIES,
                    )
                    continue
                response.raise_for_status()
                return response.json()
        except httpx.TimeoutException as exc:
            last_exc = exc
            logger.warning(
                "Backend %s %s timed out (attempt %d/%d)",
                method,
                path,
                attempt,
                _MAX_RETRIES,
            )
        except httpx.HTTPStatusError:
            raise
        except Exception as exc:
            last_exc = exc
            logger.warning(
                "Backend %s %s error (attempt %d/%d): %s",
                method,
                path,
                attempt,
                _MAX_RETRIES,
                exc,
            )

    raise RuntimeError(
        f"Backend request {method} {path} failed after {_MAX_RETRIES} attempts"
    ) from last_exc


# ── Data endpoints ──────────────────────────────────────────────────────────


def fetch_meta_data(curated: bool = False, accession_tags=None) -> list[dict]:
    return _request(
        "POST",
        "/api/data/meta",
        json={
            "curated": curated,
            "accession_tags": accession_tags,
        },
    )


def fetch_paper_data() -> list[dict]:
    return _request("GET", "/api/data/papers")


def fetch_ships(
    accession_tags=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> list[dict]:
    return _request(
        "POST",
        "/api/data/ships",
        json={
            "accession_tags": accession_tags,
            "curated": curated,
            "dereplicate": dereplicate,
            "with_sequence": with_sequence,
        },
    )


def fetch_ship_table(
    curated: bool = True,
    with_sequence: bool = False,
    with_gff_entries: bool = False,
) -> list[dict]:
    return _request(
        "POST",
        "/api/data/ship-table",
        json={
            "curated": curated,
            "with_sequence": with_sequence,
            "with_gff_entries": with_gff_entries,
        },
    )


def fetch_accession_ship(accession_tag: str) -> dict:
    return _request("GET", f"/api/data/accession/{accession_tag}/ship")


def fetch_captains(
    accession_tags=None,
    curated: bool = False,
    dereplicate: bool = True,
    with_sequence: bool = False,
) -> list[dict]:
    return _request(
        "POST",
        "/api/data/captains",
        json={
            "accession_tags": accession_tags,
            "curated": curated,
            "dereplicate": dereplicate,
            "with_sequence": with_sequence,
        },
    )


def fetch_captain_tree() -> str:
    result = _request("GET", "/api/data/captain-tree")
    return result["newick"]


def fetch_sf_data() -> list[dict]:
    return _request("GET", "/api/data/sf-data")


def get_database_version() -> str:
    result = _request("GET", "/api/data/database-version")
    return result["semantic_version"]


def set_database_version(
    semantic_version: str, description: str = "", created_by: str = "manual"
) -> bool:
    _request(
        "POST",
        "/api/data/database-version",
        json={
            "semantic_version": semantic_version,
            "description": description,
            "created_by": created_by,
        },
    )
    return True


def get_alembic_schema_version() -> str:
    result = _request("GET", "/api/data/alembic-version")
    return result["revision"]


def get_database_stats() -> dict:
    return _request("GET", "/api/data/stats")


def add_quality_tag(
    joined_ship_id: int,
    tag_type: str,
    tag_value: str | None = None,
    created_by: str = "api",
) -> int:
    result = _request(
        "POST",
        "/api/data/quality-tags",
        json={
            "joined_ship_id": joined_ship_id,
            "tag_type": tag_type,
            "tag_value": tag_value,
            "created_by": created_by,
        },
    )
    return result["id"]


def remove_quality_tag(joined_ship_id: int, tag_type: str) -> bool:
    result = _request(
        "POST",
        "/api/data/quality-tags/remove",
        json={
            "joined_ship_id": joined_ship_id,
            "tag_type": tag_type,
        },
    )
    return result["removed"]


def get_quality_tags(joined_ship_id: int) -> list[dict]:
    return _request("GET", f"/api/data/joined-ships/{joined_ship_id}/quality-tags")


def set_ship_deleted(joined_ship_id: int, deleted: bool = True) -> bool:
    result = _request(
        "POST",
        "/api/data/joined-ships/set-deleted",
        json={
            "joined_ship_id": joined_ship_id,
            "deleted": deleted,
        },
    )
    return result["updated"]


# ── Admin grid endpoints ─────────────────────────────────────────────────────


def fetch_admin_table(table_key: str) -> list[dict]:
    return _request("GET", f"/api/admin/table/{table_key}")


def admin_insert(table_key: str, col_values: dict) -> dict:
    return _request(
        "POST",
        "/api/admin/insert",
        json={"table_key": table_key, "col_values": col_values},
    )


def admin_update(table_key: str, row_id, col_id: str, new_value) -> dict:
    return _request(
        "POST",
        "/api/admin/update",
        json={
            "table_key": table_key,
            "row_id": row_id,
            "col_id": col_id,
            "new_value": new_value,
        },
    )


def run_admin_job(job_key: str) -> dict:
    return _request("POST", f"/api/admin/jobs/{job_key}/run", timeout=_BLAST_TIMEOUT)


def process_admin_submissions(sub_ids: list) -> list:
    return _request(
        "POST",
        "/api/admin/submissions/process",
        json={"sub_ids": sub_ids},
        timeout=_BLAST_TIMEOUT,
    )


def promote_admin_submission(sub_id: int) -> dict:
    return _request(
        "POST",
        "/api/admin/submissions/promote",
        json={"sub_id": sub_id},
        timeout=_BLAST_TIMEOUT,
    )


# ── Curation / annotation review endpoints ──────────────────────────────────


def fetch_annotation_queue(
    flag: int = None, assigned_to: str = None, limit: int = 50
) -> list:
    return _request(
        "POST",
        "/api/curation/queue",
        json={"flag": flag, "assigned_to": assigned_to, "limit": limit},
    )


def fetch_annotation(annotation_id: int) -> dict:
    return _request("GET", f"/api/curation/annotations/{annotation_id}")


def update_annotation(annotation_id: int, changes: dict, changed_by: str) -> dict:
    return _request(
        "POST",
        f"/api/curation/annotations/{annotation_id}/update",
        json={"changes": changes, "changed_by": changed_by},
    )


def fetch_ships_overview() -> list:
    return _request("GET", "/api/curation/ships/overview")


def fetch_ship_gene_features(joined_ship_id: int) -> dict:
    return _request("GET", f"/api/curation/ships/{joined_ship_id}/gene-features")


# ── Starfish pipeline run management ────────────────────────────────────────


def list_starfish_runs(status: str = None, limit: int = 50) -> list:
    return _request(
        "POST", "/api/starfish/runs/list", json={"status": status, "limit": limit}
    )


def get_starfish_run(run_id: int) -> dict:
    return _request("GET", f"/api/starfish/runs/{run_id}")


def create_starfish_run(run_name: str, genomes: list, **kwargs) -> dict:
    return _request(
        "POST",
        "/api/starfish/runs",
        json={"run_name": run_name, "genomes": genomes, **kwargs},
    )


def start_starfish_run(run_id: int) -> dict:
    return _request("POST", f"/api/starfish/runs/{run_id}/start")


def cancel_starfish_run(run_id: int) -> dict:
    return _request("POST", f"/api/starfish/runs/{run_id}/cancel")


def rerun_starfish_run(run_id: int) -> dict:
    return _request("POST", f"/api/starfish/runs/{run_id}/rerun")


def resume_starfish_run(run_id: int) -> dict:
    return _request("POST", f"/api/starfish/runs/{run_id}/resume")


def import_starfish_element(
    element_id: int, uploader: str = "starfish-pipeline"
) -> dict:
    return _request(
        "POST",
        f"/api/starfish/elements/{element_id}/import",
        json={"uploader": uploader},
    )


def get_starfish_run_log(run_id: int) -> str:
    result = _request("GET", f"/api/starfish/runs/{run_id}/log")
    return result["log"]


# ── BLAST / HMMER endpoints ─────────────────────────────────────────────────


def blast_search(
    query_header: str,
    query_seq: str,
    query_type: str = "nucl",
    eval_threshold: float = 0.01,
    curated: bool | None = None,
) -> dict | None:
    result = _request(
        "POST",
        "/api/blast/search",
        json={
            "query_header": query_header,
            "query_seq": query_seq,
            "query_type": query_type,
            "eval_threshold": eval_threshold,
            "curated": curated,
        },
        timeout=_BLAST_TIMEOUT,
    )
    if not result.get("ok"):
        logger.warning("Backend BLAST search failed: %s", result.get("error"))
        return None
    return result.get("result")


def hmmer_search(
    query_header: str,
    query_seq: str,
    query_type: str = "nucl",
    eval_threshold: float = 0.01,
) -> dict | None:
    result = _request(
        "POST",
        "/api/blast/hmmer",
        json={
            "query_header": query_header,
            "query_seq": query_seq,
            "query_type": query_type,
            "eval_threshold": eval_threshold,
        },
        timeout=_BLAST_TIMEOUT,
    )
    if not result.get("ok"):
        logger.warning("Backend HMMER search failed: %s", result.get("error"))
        return None
    return result.get("result")


# ── Classification workflow endpoints ───────────────────────────────────────


def _json_safe(value: Any) -> Any:
    """Recursively convert values to strict-JSON-safe equivalents.

    NaN/Inf → None; pandas DataFrames → list of record dicts (the workflow
    accepts both forms).
    """
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return None
        return value
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    to_dict = getattr(value, "to_dict", None)
    if callable(to_dict):
        try:
            return _json_safe(to_dict("records"))
        except Exception:
            return str(value)
    return value


def _inline_fasta(data: Any) -> Any:
    """Inline the local temp-file path in `fasta_file` as sequence content.

    The backend machine cannot see frontend /tmp paths, so the sequence must
    travel in the payload. The workflow accepts {"content": ...} dicts in
    fasta_file and materializes them locally.
    """
    if not isinstance(data, dict):
        return data
    fasta_file = data.get("fasta_file")
    if isinstance(fasta_file, str):
        if os.path.exists(fasta_file):
            with open(fasta_file, "r") as f:
                data["fasta_file"] = {"content": f.read()}
        elif data.get("sequence"):
            data["fasta_file"] = {"content": f">query\n{data['sequence']}\n"}
    return data


def _workflow_error_result(message: str) -> dict:
    from src.utils.blast_data import WorkflowState

    return WorkflowState.error_result(message)


def classification_workflow(
    workflow_state: dict,
    blast_data: dict | None = None,
    classification_data: dict | None = None,
    meta_dict: list | None = None,
) -> dict:
    """Run the classification workflow on the compute backend."""
    payload = _json_safe(
        {
            "workflow_state": workflow_state,
            "blast_data": _inline_fasta(blast_data) if blast_data else None,
            "classification_data": _inline_fasta(classification_data)
            if classification_data
            else None,
            "meta_dict": meta_dict,
        }
    )
    result = _request(
        "POST",
        "/api/classification/workflow",
        json=payload,
        timeout=_CLASSIFICATION_TIMEOUT,
    )
    if not result.get("ok"):
        logger.warning(
            "Backend classification workflow failed: %s", result.get("error")
        )
        return _workflow_error_result(str(result.get("error")))
    return result.get("result")

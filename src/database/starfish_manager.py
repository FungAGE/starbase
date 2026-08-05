"""starfish_manager — dual-mode facade for Starfish pipeline run management.

When BACKEND_API_URL is set: HTTP to compute backend.
When unset: direct call via starfish_manager_impl (local DB debug / monolith).

Mirrors src/database/admin_manager.py's split. Slice 4a: create/list/get
only -- start/cancel/rerun/resume land in slice 4b with actual execution.
"""

from __future__ import annotations

from src.config import backend_client

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.starfish_manager_impl as impl

        _impl = impl
    return _impl


def _via_backend(http_fn, local_fn, *args, **kwargs):
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def list_runs(status: str = None, limit: int = 50) -> list:
    def _http():
        return backend_client.list_starfish_runs(status=status, limit=limit)

    def _local():
        return _get_impl().list_runs(status=status, limit=limit)

    return _via_backend(_http, _local)


def get_run(run_id: int) -> dict:
    def _http():
        return backend_client.get_starfish_run(run_id)

    def _local():
        return _get_impl().get_run(run_id)

    return _via_backend(_http, _local)


def create_run(run_name: str, genomes: list, **kwargs) -> dict:
    def _http():
        return backend_client.create_starfish_run(run_name, genomes, **kwargs)

    def _local():
        return _get_impl().create_run(run_name, genomes, **kwargs)

    return _via_backend(_http, _local)

"""admin_jobs_manager — dual-mode facade for admin consistency/cleanup jobs.

When BACKEND_API_URL is set: HTTP to compute backend.
When unset: direct call via admin_jobs_manager_impl (local DB debug / monolith).

Mirrors src/database/admin_manager.py's split.
"""

from __future__ import annotations

from src.config import backend_client

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.admin_jobs_manager_impl as impl

        _impl = impl
    return _impl


def _via_backend(http_fn, local_fn, *args, **kwargs):
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def run_job(job_key: str) -> dict:
    def _http():
        return backend_client.run_admin_job(job_key)

    def _local():
        return _get_impl().run_job(job_key)

    return _via_backend(_http, _local)

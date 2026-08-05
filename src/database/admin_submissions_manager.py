"""admin_submissions_manager — dual-mode facade for the process/promote workflow.

When BACKEND_API_URL is set: HTTP to compute backend.
When unset: direct call via admin_submissions_manager_impl (local DB debug / monolith).

Mirrors src/database/admin_manager.py's split.
"""

from __future__ import annotations

from src.config import backend_client

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.admin_submissions_manager_impl as impl

        _impl = impl
    return _impl


def _via_backend(http_fn, local_fn, *args, **kwargs):
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def process_submissions(sub_ids: list) -> list:
    def _http():
        return backend_client.process_admin_submissions(sub_ids)

    def _local():
        return _get_impl().process_submissions(sub_ids)

    return _via_backend(_http, _local)


def promote_submission(sub_id: int):
    def _http():
        result = backend_client.promote_admin_submission(sub_id)
        return result["success"], result["accession"], result["error"]

    def _local():
        return _get_impl().promote_submission(sub_id)

    return _via_backend(_http, _local)

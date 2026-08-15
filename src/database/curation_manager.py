"""curation_manager — dual-mode facade for the curation/annotation-review workflow.

When BACKEND_API_URL is set: HTTP to compute backend.
When unset: direct call via curation_manager_impl (local DB debug / monolith).

Mirrors src/database/admin_manager.py's split.
"""

from __future__ import annotations

from src.config import backend_client

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.curation_manager_impl as impl

        _impl = impl
    return _impl


def _via_backend(http_fn, local_fn, *args, **kwargs):
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def fetch_annotation_queue(
    flag: int = None, assigned_to: str = None, limit: int = 50
) -> list:
    def _http():
        return backend_client.fetch_annotation_queue(
            flag=flag, assigned_to=assigned_to, limit=limit
        )

    def _local():
        return _get_impl().fetch_annotation_queue(
            flag=flag, assigned_to=assigned_to, limit=limit
        )

    return _via_backend(_http, _local)


def fetch_annotation(annotation_id: int) -> dict:
    def _http():
        return backend_client.fetch_annotation(annotation_id)

    def _local():
        return _get_impl().fetch_annotation(annotation_id)

    return _via_backend(_http, _local)


def update_annotation(annotation_id: int, changes: dict, changed_by: str) -> dict:
    def _http():
        return backend_client.update_annotation(annotation_id, changes, changed_by)

    def _local():
        return _get_impl().update_annotation(annotation_id, changes, changed_by)

    return _via_backend(_http, _local)


def fetch_ships_overview() -> list:
    def _http():
        return backend_client.fetch_ships_overview()

    def _local():
        return _get_impl().fetch_ships_overview()

    return _via_backend(_http, _local)


def fetch_ship_gene_features(joined_ship_id: int) -> dict:
    def _http():
        return backend_client.fetch_ship_gene_features(joined_ship_id)

    def _local():
        return _get_impl().fetch_ship_gene_features(joined_ship_id)

    return _via_backend(_http, _local)

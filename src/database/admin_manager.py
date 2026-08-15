"""admin_manager — dual-mode facade for admin grid read/write operations.

When BACKEND_API_URL is set: HTTP to compute backend.
When unset: direct SQLAlchemy via admin_manager_impl (local DB debug / monolith).

Mirrors src/database/sql_manager.py's split.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from src.config import backend_client

_impl = None


def _get_impl():
    global _impl
    if _impl is None:
        import src.database.admin_manager_impl as impl

        _impl = impl
    return _impl


def _via_backend(http_fn, local_fn, *args, **kwargs):
    if backend_client.is_configured():
        return http_fn(*args, **kwargs)
    return local_fn(*args, **kwargs)


def fetch_admin_table(table_key: str) -> pd.DataFrame:
    def _http():
        records = backend_client.fetch_admin_table(table_key)
        return pd.DataFrame.from_records(records) if records else pd.DataFrame()

    def _local():
        return _get_impl().fetch_admin_table(table_key)

    return _via_backend(_http, _local)


def do_insert(table_key: str, col_values: dict[str, Any]):
    def _http():
        result = backend_client.admin_insert(table_key, col_values)
        return result["success"], result["error"], result["new_id"]

    def _local():
        return _get_impl().do_insert(table_key, col_values)

    return _via_backend(_http, _local)


def do_update(table_key: str, row_id, col_id: str, new_value):
    def _http():
        result = backend_client.admin_update(table_key, row_id, col_id, new_value)
        return result["success"], result["error"]

    def _local():
        return _get_impl().do_update(table_key, row_id, col_id, new_value)

    return _via_backend(_http, _local)

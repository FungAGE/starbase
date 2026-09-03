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


def start_run(run_id: int) -> dict:
    def _http():
        return backend_client.start_starfish_run(run_id)

    def _local():
        return _get_impl().start_run(run_id)

    return _via_backend(_http, _local)


def cancel_run(run_id: int) -> dict:
    def _http():
        return backend_client.cancel_starfish_run(run_id)

    def _local():
        return _get_impl().cancel_run(run_id)

    return _via_backend(_http, _local)


def rerun_run(run_id: int) -> dict:
    def _http():
        return backend_client.rerun_starfish_run(run_id)

    def _local():
        return _get_impl().rerun_run(run_id)

    return _via_backend(_http, _local)


def resume_run(run_id: int) -> dict:
    def _http():
        return backend_client.resume_starfish_run(run_id)

    def _local():
        return _get_impl().resume_run(run_id)

    return _via_backend(_http, _local)


def import_element_to_submission(
    element_id: int, uploader: str = "starfish-pipeline"
) -> dict:
    def _http():
        return backend_client.import_starfish_element(element_id, uploader=uploader)

    def _local():
        return _get_impl().import_element_to_submission(element_id, uploader=uploader)

    return _via_backend(_http, _local)


def get_run_log(run_id: int) -> str:
    def _http():
        return backend_client.get_starfish_run_log(run_id)

    def _local():
        return _get_impl().get_run_log(run_id)

    return _via_backend(_http, _local)


def list_visualizations(run_id: int) -> dict:
    def _http():
        return backend_client.list_starfish_visualizations(run_id)

    def _local():
        return _get_impl().list_visualizations(run_id)

    return _via_backend(_http, _local)


def get_visualization_file(run_id: int, section: str, filename: str) -> tuple:
    def _http():
        return backend_client.get_starfish_visualization_file(run_id, section, filename)

    def _local():
        return _get_impl().get_visualization_file(run_id, section, filename)

    return _via_backend(_http, _local)


def delete_element(element_id: int) -> dict:
    def _http():
        return backend_client.delete_starfish_element(element_id)

    def _local():
        return _get_impl().delete_element(element_id)

    return _via_backend(_http, _local)


def update_element(
    element_id: int,
    family: str = None,
    navis: str = None,
    haplotype: str = None,
    confidence: str = None,
    notes: str = None,
) -> dict:
    def _http():
        return backend_client.update_starfish_element(
            element_id,
            family=family,
            navis=navis,
            haplotype=haplotype,
            confidence=confidence,
            notes=notes,
        )

    def _local():
        return _get_impl().update_element(
            element_id,
            family=family,
            navis=navis,
            haplotype=haplotype,
            confidence=confidence,
            notes=notes,
        )

    return _via_backend(_http, _local)

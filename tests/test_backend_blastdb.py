"""Tests for the scheduled BLAST DB rebuild (backend/tasks/blastdb.py)."""

import pytest
from fastapi.testclient import TestClient

from backend.main import app as backend_app
from backend.tasks import blastdb as blastdb_tasks
from src.config import celery_config

client = TestClient(backend_app)


@pytest.fixture(autouse=True)
def _no_api_key(monkeypatch):
    monkeypatch.delenv("BACKEND_API_KEY", raising=False)


def _allow_rebuild(monkeypatch):
    monkeypatch.setattr(blastdb_tasks, "IS_DEV", False)
    monkeypatch.setattr(
        blastdb_tasks.shutil, "which", lambda name: "/usr/bin/makeblastdb"
    )


# ── Core rebuild logic ──────────────────────────────────────────────────────


def test_rebuild_skips_in_dev(monkeypatch):
    monkeypatch.setattr(blastdb_tasks, "IS_DEV", True)
    result = blastdb_tasks.rebuild_blast_dbs()
    assert result == {"rebuilt": False, "skipped": "dev mode"}


def test_rebuild_skips_without_makeblastdb(monkeypatch):
    monkeypatch.setattr(blastdb_tasks, "IS_DEV", False)
    monkeypatch.setattr(blastdb_tasks.shutil, "which", lambda name: None)
    result = blastdb_tasks.rebuild_blast_dbs()
    assert result["rebuilt"] is False
    assert "makeblastdb" in result["skipped"]


def test_rebuild_success(monkeypatch):
    _allow_rebuild(monkeypatch)
    calls = {}

    def fake_create_dbs():
        calls["created"] = True

    monkeypatch.setattr("src.database.blastdb.create_dbs", fake_create_dbs)
    monkeypatch.setattr("src.database.blastdb.blast_db_exists", lambda path: True)

    result = blastdb_tasks.rebuild_blast_dbs()
    assert calls["created"] is True
    assert result["rebuilt"] is True
    assert result["databases"] == {
        "ships_all": True,
        "ships_curated": True,
        "captains": True,
    }
    assert result["duration_s"] >= 0


def test_rebuild_failure_returns_error_shape(monkeypatch):
    _allow_rebuild(monkeypatch)

    def boom():
        raise RuntimeError("makeblastdb exploded")

    monkeypatch.setattr("src.database.blastdb.create_dbs", boom)
    result = blastdb_tasks.rebuild_blast_dbs()
    assert result["rebuilt"] is False
    assert "makeblastdb exploded" in result["error"]


def test_rebuild_verification_failure_flagged(monkeypatch):
    _allow_rebuild(monkeypatch)
    monkeypatch.setattr("src.database.blastdb.create_dbs", lambda: None)
    monkeypatch.setattr("src.database.blastdb.blast_db_exists", lambda path: False)
    result = blastdb_tasks.rebuild_blast_dbs()
    assert result["rebuilt"] is False
    assert "verification failed" in result["error"]


# ── Celery wiring ───────────────────────────────────────────────────────────


def test_task_registered_with_celery():
    if not celery_config.CELERY_AVAILABLE:
        pytest.skip("celery not installed")
    assert "rebuild_blast_dbs" in celery_config.celery.tasks


def test_beat_schedule_has_daily_rebuild():
    if not celery_config.CELERY_AVAILABLE:
        pytest.skip("celery not installed")
    schedule = celery_config.celery.conf.beat_schedule
    assert "rebuild-blast-db-daily" in schedule
    assert schedule["rebuild-blast-db-daily"]["task"] == "rebuild_blast_dbs"


def test_worker_imports_include_backend_tasks():
    if not celery_config.CELERY_AVAILABLE:
        pytest.skip("celery not installed")
    imports = list(celery_config.celery.conf.imports)
    assert "backend.tasks.blastdb" in imports
    assert "backend.tasks.starfish" in imports


# ── Admin endpoints ─────────────────────────────────────────────────────────


def test_blastdb_status_endpoint(monkeypatch):
    monkeypatch.setattr("src.database.blastdb.blast_db_exists", lambda path: True)
    monkeypatch.setattr("src.database.blastdb.sourmash_sig_exists", lambda path: False)

    response = client.get("/api/admin/blastdb/status")
    assert response.status_code == 200
    body = response.json()
    for name in ("ships_all", "ships_curated", "captains"):
        assert body[name]["built"] is True
        assert body[name]["signatures"] is False
        assert "path" in body[name]
        assert "updated_at" in body[name]


def test_blastdb_rebuild_dispatches_task(monkeypatch):
    captured = {}

    def fake_dispatch():
        captured["called"] = True
        return {"dispatched": True, "task_id": "abc123"}

    monkeypatch.setattr("backend.tasks.blastdb.dispatch_blastdb_rebuild", fake_dispatch)

    response = client.post("/api/admin/blastdb/rebuild")
    assert response.status_code == 200
    assert response.json() == {"dispatched": True, "task_id": "abc123"}
    assert captured["called"] is True

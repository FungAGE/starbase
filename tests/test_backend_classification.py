import pytest
from fastapi.testclient import TestClient

import src.config.backend_client as backend_client
from backend.main import app as backend_app
from src.tasks import run_classification_workflow_sync

client = TestClient(backend_app)


@pytest.fixture(autouse=True)
def _no_api_key(monkeypatch):
    monkeypatch.delenv("BACKEND_API_KEY", raising=False)


# ── Backend endpoint ────────────────────────────────────────────────────────


class _FakeAsyncResult:
    def __init__(self, result):
        self._result = result

    def get(self, timeout=None):
        return self._result


def test_workflow_endpoint_direct_call(monkeypatch):
    import backend.routers.classification as classification_router

    captured = {}

    def fake_task(
        workflow_state, blast_data=None, classification_data=None, meta_dict=None
    ):
        captured.update(
            {
                "workflow_state": workflow_state,
                "blast_data": blast_data,
                "classification_data": classification_data,
                "meta_dict": meta_dict,
            }
        )
        return {"complete": True, "found_match": True, "match_stage": "exact"}

    monkeypatch.setattr(
        classification_router, "run_classification_workflow_task", fake_task
    )

    response = client.post(
        "/api/classification/workflow",
        json={
            "workflow_state": {"complete": False, "stages": {}},
            "blast_data": {"seq_type": "nucl", "fasta_file": {"content": ">q\nATGC"}},
            "classification_data": {"seq_type": "nucl"},
            "meta_dict": [{"accession_tag": "SSA002851"}],
        },
    )

    assert response.status_code == 200
    body = response.json()
    assert body["ok"] is True
    assert body["result"]["found_match"] is True
    assert captured["workflow_state"] == {"complete": False, "stages": {}}
    assert captured["blast_data"]["fasta_file"]["content"] == ">q\nATGC"
    assert captured["meta_dict"] == [{"accession_tag": "SSA002851"}]


def test_workflow_endpoint_celery_apply(monkeypatch):
    import backend.routers.classification as classification_router

    captured = {}

    class FakeCeleryTask:
        def apply(self, args=(), kwargs=None):
            captured["args"] = args
            captured["kwargs"] = kwargs
            return _FakeAsyncResult({"complete": True, "found_match": False})

    monkeypatch.setattr(
        classification_router, "run_classification_workflow_task", FakeCeleryTask()
    )

    response = client.post(
        "/api/classification/workflow",
        json={"workflow_state": {"complete": False}, "blast_data": None},
    )

    assert response.status_code == 200
    assert response.json()["ok"] is True
    assert captured["args"] == ({"complete": False},)
    assert captured["kwargs"]["blast_data"] is None


def test_workflow_endpoint_task_exception(monkeypatch):
    import backend.routers.classification as classification_router

    def fake_task(
        workflow_state, blast_data=None, classification_data=None, meta_dict=None
    ):
        raise RuntimeError("boom")

    monkeypatch.setattr(
        classification_router, "run_classification_workflow_task", fake_task
    )

    response = client.post(
        "/api/classification/workflow",
        json={"workflow_state": {"complete": False}},
    )

    assert response.status_code == 200
    body = response.json()
    assert body["ok"] is False
    assert "boom" in body["error"]
    assert body["result"] is None


def test_workflow_endpoint_requires_workflow_state():
    response = client.post("/api/classification/workflow", json={"blast_data": {}})
    assert response.status_code == 422


# ── Frontend dispatch ───────────────────────────────────────────────────────


def test_sync_dispatches_to_backend_when_configured(monkeypatch, tmp_path):
    fasta_file = tmp_path / "query.fa"
    fasta_file.write_text(">query\nATGCATGC\n")

    captured = {}

    def fake_request(method, path, *, json=None, timeout=None):
        captured.update(
            {"method": method, "path": path, "json": json, "timeout": timeout}
        )
        return {"ok": True, "result": {"complete": True, "found_match": True}}

    monkeypatch.setattr(backend_client, "_BACKEND_API_URL", "http://backend:8001")
    monkeypatch.setattr(backend_client, "_request", fake_request)

    result = run_classification_workflow_sync(
        workflow_state={"complete": False, "stages": {}},
        blast_data={
            "seq_type": "nucl",
            "fasta_file": str(fasta_file),
            "blast_df": [{"evalue": float("nan"), "pident": 100.0}],
            "sequence": None,
        },
        classification_data={"seq_type": "nucl", "fasta_file": str(fasta_file)},
        meta_dict=[{"accession_tag": "SSA002851", "familyName": None}],
    )

    assert result == {"complete": True, "found_match": True}
    assert captured["path"] == "/api/classification/workflow"
    assert captured["timeout"] == backend_client._CLASSIFICATION_TIMEOUT
    payload = captured["json"]
    # fasta path inlined as content
    assert payload["blast_data"]["fasta_file"] == {"content": ">query\nATGCATGC\n"}
    assert payload["classification_data"]["fasta_file"] == {
        "content": ">query\nATGCATGC\n"
    }
    # NaN sanitized so the payload is strict JSON
    assert payload["blast_data"]["blast_df"][0]["evalue"] is None


def test_sync_uses_local_path_when_not_configured(monkeypatch):
    import src.utils.classification_utils as classification_utils

    monkeypatch.setattr(backend_client, "_BACKEND_API_URL", "")

    captured = {}

    def fake_workflow(workflow_state, blast_data, classification_data, meta_dict=None):
        captured["called"] = True
        return {"complete": True, "found_match": False, "local": True}

    monkeypatch.setattr(
        classification_utils, "run_classification_workflow", fake_workflow
    )

    result = run_classification_workflow_sync(
        workflow_state={"complete": False, "stages": {}},
        blast_data={"seq_type": "nucl", "fasta_file": ">q\nATGC"},
    )

    assert captured["called"] is True
    assert result["local"] is True


def test_backend_failure_returns_error_shape(monkeypatch):
    monkeypatch.setattr(backend_client, "_BACKEND_API_URL", "http://backend:8001")
    monkeypatch.setattr(
        backend_client,
        "_request",
        lambda method, path, *, json=None, timeout=None: {
            "ok": False,
            "error": "Workflow returned no result",
        },
    )

    result = run_classification_workflow_sync(
        workflow_state={"complete": False},
        blast_data={"seq_type": "nucl", "sequence": "ATGC"},
    )

    assert result["complete"] is True
    assert result["status"] == "failed"
    assert "Workflow returned no result" in result["error"]
    assert result["found_match"] is False


def test_json_safe_handles_nested_and_tuples():
    nan = float("nan")
    inf = float("inf")
    assert backend_client._json_safe(nan) is None
    assert backend_client._json_safe(inf) is None
    assert backend_client._json_safe(1.5) == 1.5
    assert backend_client._json_safe({"a": [nan, (inf, None)], "b": 2}) == {
        "a": [None, [None, None]],
        "b": 2,
    }


def test_json_safe_converts_dataframe():
    import pandas as pd

    df = pd.DataFrame({"a": [1.0, float("nan")], "b": ["x", "y"]})
    assert backend_client._json_safe(df) == [
        {"a": 1.0, "b": "x"},
        {"a": None, "b": "y"},
    ]

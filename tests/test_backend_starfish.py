"""Tests for starfish run/element/visualization backend (impl + router)."""

from contextlib import contextmanager

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import src.database.starfish_manager_impl as impl
from backend.main import app as backend_app
from src.database.models.schema import (
    Base,
    StarfishElement,
    StarfishRun,
    StarfishRunGenome,
)

client = TestClient(backend_app)


@pytest.fixture(autouse=True)
def _no_api_key(monkeypatch):
    monkeypatch.delenv("BACKEND_API_KEY", raising=False)


@pytest.fixture()
def starfish_db(monkeypatch, tmp_path):
    """File-backed SQLite (shared across TestClient threads) + tmp runs dir,
    swapped into the impl module."""
    engine = create_engine(f"sqlite:///{tmp_path / 'test_starfish.db'}")
    Base.metadata.create_all(engine)
    TestingSession = sessionmaker(bind=engine)

    @contextmanager
    def fake_session():
        session = TestingSession()
        try:
            yield session
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

    monkeypatch.setattr(impl, "get_starbase_session", fake_session)
    monkeypatch.setattr(impl, "STARFISH_RUNS_DIR", str(tmp_path / "starfish_runs"))

    def add_element(run_id, genome_id, element_id="elem-1", **kw):
        with fake_session() as s:
            e = StarfishElement(
                element_id=element_id,
                run_id=run_id,
                genome_id=genome_id,
                contig_id=kw.get("contig_id", "ctg1"),
                start=kw.get("start", 100),
                end=kw.get("end", 200),
                strand=kw.get("strand", "+"),
                family=kw.get("family"),
                navis=kw.get("navis"),
                haplotype=kw.get("haplotype"),
                confidence=kw.get("confidence"),
                notes=kw.get("notes"),
                imported_submission_id=kw.get("imported_submission_id"),
            )
            s.add(e)
            s.flush()
            return e.id

    def set_counters(run_id, num_elements_found=1, num_elements=1):
        with fake_session() as s:
            r = s.query(StarfishRun).filter_by(id=run_id).first()
            r.num_elements_found = num_elements_found
            g = s.query(StarfishRunGenome).filter_by(run_id=run_id).first()
            g.num_elements = num_elements

    def get_genome_id(run_id):
        with fake_session() as s:
            g = s.query(StarfishRunGenome).filter_by(run_id=run_id).first()
            return g.id if g else None

    def get_run_counters(run_id):
        with fake_session() as s:
            r = s.query(StarfishRun).filter_by(id=run_id).first()
            g = s.query(StarfishRunGenome).filter_by(run_id=run_id).first()
            return (
                r.num_elements_found if r else None,
                g.num_elements if g else None,
            )

    def get_element(element_id):
        with fake_session() as s:
            e = s.query(StarfishElement).filter_by(id=element_id).first()
            if e is None:
                return None
            return {
                "id": e.id,
                "element_id": e.element_id,
                "family": e.family,
                "navis": e.navis,
                "haplotype": e.haplotype,
                "confidence": e.confidence,
                "notes": e.notes,
                "imported_submission_id": e.imported_submission_id,
            }

    yield {
        "add_element": add_element,
        "set_counters": set_counters,
        "get_genome_id": get_genome_id,
        "get_run_counters": get_run_counters,
        "get_element": get_element,
    }


def _make_run(name="vizrun"):
    return impl.create_run(
        name,
        [
            {
                "genome_id": "g1",
                "fna_path": "/data/g1.fna",
                "gff3_path": "/data/g1.gff3",
            }
        ],
    )


def _viz_dir(run, section):
    import os

    return os.path.join(run["output_dir"], run["run_name"], section)


# ── impl: list_visualizations ───────────────────────────────────────────────


def test_list_visualizations_empty(starfish_db):
    run = _make_run()
    result = impl.list_visualizations(run["id"])
    assert result == {"locusViz_files": [], "pairViz_files": []}


def test_list_visualizations_files(starfish_db):
    import os

    run = _make_run()
    locus = _viz_dir(run, "locusViz")
    pair = _viz_dir(run, "pairViz")
    os.makedirs(locus)
    os.makedirs(pair)
    for name in ("b.png", "a.png"):
        with open(os.path.join(locus, name), "wb") as f:
            f.write(b"png-bytes")
    with open(os.path.join(pair, "pair1.png"), "wb") as f:
        f.write(b"png-bytes")
    os.makedirs(os.path.join(locus, "subdir"))

    result = impl.list_visualizations(run["id"])
    assert result["locusViz_files"] == ["a.png", "b.png"]
    assert result["pairViz_files"] == ["pair1.png"]


def test_list_visualizations_unknown_run(starfish_db):
    with pytest.raises(ValueError):
        impl.list_visualizations(9999)


# ── impl: get_visualization_file ────────────────────────────────────────────


def test_get_visualization_file(starfish_db):
    import os

    run = _make_run()
    locus = _viz_dir(run, "locusViz")
    os.makedirs(locus)
    payload = b"\x89PNG-fake"
    with open(os.path.join(locus, "x.png"), "wb") as f:
        f.write(payload)

    data, media_type = impl.get_visualization_file(run["id"], "locusViz", "x.png")
    assert data == payload
    assert media_type == "image/png"

    with open(os.path.join(locus, "notes.txt"), "w") as f:
        f.write("hello")
    data, media_type = impl.get_visualization_file(run["id"], "locusViz", "notes.txt")
    assert data == b"hello"
    assert media_type == "application/octet-stream"


def test_get_visualization_file_rejects_traversal(starfish_db):
    run = _make_run()
    for bad in ("../evil.png", "..", "", "a/b.png"):
        with pytest.raises(ValueError):
            impl.get_visualization_file(run["id"], "locusViz", bad)


def test_get_visualization_file_bad_section(starfish_db):
    run = _make_run()
    with pytest.raises(ValueError):
        impl.get_visualization_file(run["id"], "bogus", "x.png")


def test_get_visualization_file_missing(starfish_db):
    run = _make_run()
    with pytest.raises(ValueError):
        impl.get_visualization_file(run["id"], "locusViz", "nope.png")


def test_get_visualization_file_unknown_run(starfish_db):
    with pytest.raises(ValueError):
        impl.get_visualization_file(9999, "locusViz", "x.png")


# ── impl: delete_element ────────────────────────────────────────────────────


def test_delete_element(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    starfish_db["set_counters"](run["id"])
    element_id = starfish_db["add_element"](run["id"], genome_id, "elem-del")

    result = impl.delete_element(element_id)
    assert result["deleted"] == "elem-del"
    assert starfish_db["get_element"](element_id) is None
    assert starfish_db["get_run_counters"](run["id"]) == (0, 0)


def test_delete_element_imported_blocked(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    element_id = starfish_db["add_element"](
        run["id"], genome_id, "elem-imp", imported_submission_id=42
    )
    with pytest.raises(ValueError, match="already imported"):
        impl.delete_element(element_id)
    assert starfish_db["get_element"](element_id) is not None


def test_delete_element_unknown(starfish_db):
    with pytest.raises(ValueError):
        impl.delete_element(9999)


# ── impl: update_element ────────────────────────────────────────────────────


def test_update_element(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    element_id = starfish_db["add_element"](run["id"], genome_id, "elem-edit")

    updated = impl.update_element(
        element_id,
        family="Prometheus",
        navis="Prometheus",
        haplotype="2",
        confidence="High",
        notes="looks good",
    )
    assert updated["family"] == "Prometheus"
    assert updated["navis"] == "Prometheus"
    assert updated["haplotype"] == "2"
    assert updated["confidence"] == "High"
    assert updated["notes"] == "looks good"

    cleared = impl.update_element(
        element_id, family=None, navis=None, haplotype=None, confidence=None, notes=None
    )
    assert cleared["family"] is None
    assert cleared["notes"] is None


def test_update_element_unknown(starfish_db):
    with pytest.raises(ValueError):
        impl.update_element(9999, family="x")


# ── router endpoints ────────────────────────────────────────────────────────


def test_visualizations_endpoint(starfish_db):
    import os

    run = _make_run()
    locus = _viz_dir(run, "locusViz")
    os.makedirs(locus)
    with open(os.path.join(locus, "x.png"), "wb") as f:
        f.write(b"png")

    response = client.get(f"/api/starfish/runs/{run['id']}/visualizations")
    assert response.status_code == 200
    assert response.json() == {"locusViz_files": ["x.png"], "pairViz_files": []}


def test_visualizations_endpoint_unknown_run(starfish_db):
    response = client.get("/api/starfish/runs/9999/visualizations")
    assert response.status_code == 404


def test_viz_file_endpoint(starfish_db):
    import os

    run = _make_run()
    locus = _viz_dir(run, "locusViz")
    os.makedirs(locus)
    payload = b"png-bytes-here"
    with open(os.path.join(locus, "x.png"), "wb") as f:
        f.write(payload)

    response = client.get(f"/api/starfish/runs/{run['id']}/viz-files/locusViz/x.png")
    assert response.status_code == 200
    assert response.content == payload
    assert response.headers["content-type"].startswith("image/png")


def test_viz_file_endpoint_404s(starfish_db):
    run = _make_run()
    assert (
        client.get(
            f"/api/starfish/runs/{run['id']}/viz-files/locusViz/nope.png"
        ).status_code
        == 404
    )
    assert (
        client.get(f"/api/starfish/runs/{run['id']}/viz-files/bogus/x.png").status_code
        == 404
    )
    assert (
        client.get("/api/starfish/runs/9999/viz-files/locusViz/x.png").status_code
        == 404
    )


def test_delete_element_endpoint(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    element_id = starfish_db["add_element"](run["id"], genome_id, "elem-api")

    response = client.delete(f"/api/starfish/elements/{element_id}")
    assert response.status_code == 200
    assert response.json() == {"deleted": "elem-api"}

    response = client.delete(f"/api/starfish/elements/{element_id}")
    assert response.status_code == 400


def test_delete_element_endpoint_imported_blocked(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    element_id = starfish_db["add_element"](
        run["id"], genome_id, "elem-api-imp", imported_submission_id=7
    )
    response = client.delete(f"/api/starfish/elements/{element_id}")
    assert response.status_code == 400


def test_update_element_endpoint(starfish_db):
    run = _make_run()
    genome_id = starfish_db["get_genome_id"](run["id"])
    element_id = starfish_db["add_element"](run["id"], genome_id, "elem-api-edit")

    response = client.patch(
        f"/api/starfish/elements/{element_id}",
        json={
            "family": "F",
            "navis": "N",
            "haplotype": "H",
            "confidence": "C",
            "notes": "M",
        },
    )
    assert response.status_code == 200
    body = response.json()
    assert body["family"] == "F"
    assert body["navis"] == "N"
    assert body["haplotype"] == "H"
    assert body["confidence"] == "C"
    assert body["notes"] == "M"

    response = client.patch("/api/starfish/elements/9999", json={"family": "F"})
    assert response.status_code == 400

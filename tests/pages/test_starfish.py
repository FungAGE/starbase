"""Layout smoke tests for the Starfish Dash page (no live app needed)."""

from unittest.mock import patch

from dash.development.base_component import Component

with patch("dash.register_page"):
    import src.pages.starfish as starfish_page


def _collect_ids(component):
    ids = []
    cid = getattr(component, "id", None)
    if cid is not None:
        ids.append(cid)
    children = getattr(component, "children", None)
    if children is None:
        return ids
    items = children if isinstance(children, (list, tuple)) else [children]
    for child in items:
        if isinstance(child, Component):
            ids.extend(_collect_ids(child))
    return ids


def _element_row(**overrides):
    row = {
        "id": 1,
        "element_id": "e1",
        "genome_id": "g1",
        "contig_id": "c1",
        "start": 1,
        "end": 5,
        "strand": "+",
        "length": 5,
        "family": "F",
        "navis": None,
        "haplotype": None,
        "confidence": "ref",
        "notes": None,
        "imported_submission_id": None,
    }
    row.update(overrides)
    return row


def test_elements_grid_shape():
    grid = starfish_page._elements_grid([])
    assert grid.id == "starfish-elements-grid"

    grid = starfish_page._elements_grid([_element_row()])
    row = grid.rowData[0]
    assert row["position_display"] == "1-5"
    assert row["import_display"] == "not imported"

    grid = starfish_page._elements_grid([_element_row(imported_submission_id=9)])
    assert grid.rowData[0]["import_display"] == "imported"


def test_element_actions():
    assert isinstance(starfish_page._element_actions(False), Component)
    ids = _collect_ids(starfish_page._element_actions(True))
    assert "starfish-element-import-btn" in ids
    assert "starfish-element-edit-btn" in ids
    assert "starfish-element-delete-btn" in ids


def test_viz_section_buttons():
    section = starfish_page._viz_section(
        {"locusViz_files": ["a.png", "b.png"], "pairViz_files": []}
    )
    ids = _collect_ids(section)
    assert {"type": "starfish-viz-open", "section": "locusViz", "index": 0} in ids
    assert {"type": "starfish-viz-open", "section": "locusViz", "index": 1} in ids
    assert not any(isinstance(i, dict) and i.get("section") == "pairViz" for i in ids)


def test_viz_section_empty():
    section = starfish_page._viz_section({"locusViz_files": [], "pairViz_files": []})
    assert isinstance(section, Component)


def test_detail_section_structure():
    ids = _collect_ids(starfish_page._build_detail_section())
    assert "starfish-detail-accordion" in ids
    assert "starfish-detail-viz" in ids
    assert "starfish-detail-elements-hint" in ids
    assert "starfish-detail-elements-actions" in ids
    assert "starfish-detail-elements" in ids


def test_modals_present_and_unique():
    ids = _collect_ids(starfish_page._viz_modal())
    assert "starfish-viz-modal" in ids
    assert "starfish-viz-download-btn" in ids
    ids = _collect_ids(starfish_page._edit_element_modal())
    assert "starfish-edit-modal" in ids
    assert "starfish-edit-save-btn" in ids
    ids = _collect_ids(starfish_page._delete_element_modal())
    assert "starfish-delete-modal" in ids
    assert "starfish-delete-confirm-btn" in ids


def test_no_duplicate_ids(monkeypatch):
    monkeypatch.setattr(starfish_page, "ADMIN_TOKEN", "test-token")
    children, authed = starfish_page.render_starfish_content("?token=test-token")
    assert authed is True
    all_ids = _collect_ids(starfish_page.layout) + _collect_ids(children)
    dupes = {i for i in all_ids if all_ids.count(i) > 1}
    assert not dupes, f"duplicate component ids: {dupes}"


def test_unauthenticated():
    _, authed = starfish_page.render_starfish_content("?token=wrong")
    assert authed is False

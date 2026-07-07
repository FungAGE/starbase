"""Unit tests for SubmissionProcessor._insert_to_database branch logic."""

from unittest.mock import MagicMock, patch

import pytest

from src.utils.submission_utils import DuplicateInfo, SubmissionProcessor


def _minimal_data():
    return {
        "sequence": "ATGCATGCATGCATGC",
        "starshipID": "test_ship",
        "evidence": "manual",
        "source": "web_submission_promoted",
        "genus": "Escherichia",
        "species": "coli",
    }


@pytest.fixture
def processor():
    return SubmissionProcessor()


def test_insert_creates_new_ship_when_not_duplicate(processor):
    """New sequences must create a ship; must not dereference ship=None."""
    duplicate_info = DuplicateInfo(is_duplicate=False)
    mock_session = MagicMock()
    added = []

    def track_add(obj):
        added.append(obj)
        if hasattr(obj, "id") and obj.id is None:
            obj.id = len(added)

    mock_session.add.side_effect = track_add
    mock_session.query.return_value.filter.return_value.first.return_value = None

    with (
        patch("src.utils.submission_utils.StarbaseSession", return_value=mock_session),
        patch("src.utils.submission_utils.fetch_ships", return_value=[]),
        patch(
            "src.utils.submission_utils.assign_accession",
            return_value=("SSA000099", False),
        ),
        patch("src.utils.submission_utils.ensure_ship_has_ssb", return_value=None),
        patch.object(processor, "_get_or_create_taxonomy", return_value=1),
        patch.object(processor, "_get_or_create_family", return_value=None),
        patch.object(processor, "_get_or_create_navis", return_value=None),
        patch.object(processor, "_get_or_create_haplotype", return_value=None),
        patch.object(processor, "_get_or_create_genome", return_value=None),
        patch.object(processor, "_add_gff_entries"),
        patch.object(processor, "_create_starship_features"),
    ):
        ship_id, accession_tag = processor._insert_to_database(
            _minimal_data(), duplicate_info
        )

    assert accession_tag == "SSA000099"
    assert ship_id is not None
    mock_session.commit.assert_called_once()


def test_insert_uses_staging_accession_when_trusted(processor):
    """Promote path should reuse staging accession_tag instead of assign_accession."""
    duplicate_info = DuplicateInfo(is_duplicate=False)
    data = {
        **_minimal_data(),
        "trust_staging": True,
        "staging_accession_tag": "SSA000042",
    }
    mock_session = MagicMock()
    added = []

    def track_add(obj):
        added.append(obj)
        if hasattr(obj, "id") and obj.id is None:
            obj.id = len(added)

    mock_session.add.side_effect = track_add
    mock_session.query.return_value.filter.return_value.first.return_value = None

    with (
        patch("src.utils.submission_utils.StarbaseSession", return_value=mock_session),
        patch(
            "src.utils.submission_utils.assign_accession",
        ) as mock_assign,
        patch("src.utils.submission_utils.ensure_ship_has_ssb", return_value=None),
        patch.object(processor, "_get_or_create_taxonomy", return_value=1),
        patch.object(processor, "_get_or_create_family", return_value=None),
        patch.object(processor, "_get_or_create_navis", return_value=None),
        patch.object(processor, "_get_or_create_haplotype", return_value=None),
        patch.object(processor, "_get_or_create_genome", return_value=None),
        patch.object(processor, "_add_gff_entries"),
        patch.object(processor, "_create_starship_features"),
    ):
        ship_id, accession_tag = processor._insert_to_database(data, duplicate_info)

    mock_assign.assert_not_called()
    assert accession_tag == "SSA000042"
    assert ship_id is not None
    ship_accession_adds = [
        obj for obj in added if obj.__class__.__name__ == "ShipAccessions"
    ]
    assert not ship_accession_adds


def test_process_submission_trust_staging_skips_duplicate_check(processor):
    with (
        patch(
            "src.utils.submission_utils.validate_submission",
        ) as mock_validate,
        patch(
            "src.utils.submission_utils.check_sequence_duplicate",
        ) as mock_duplicate,
        patch.object(
            processor,
            "_insert_to_database",
            return_value=(99, "SSA000042"),
        ),
    ):
        result = processor.process_submission(
            {
                **_minimal_data(),
                "trust_staging": True,
                "staging_accession_tag": "SSA000042",
            }
        )

    mock_validate.assert_not_called()
    mock_duplicate.assert_not_called()
    assert result["success"] is True
    assert result["accession"] == "SSA000042"


def test_create_ship_reuses_existing_for_same_accession():
    """Do not create duplicate ships when accession already has a matching ship."""
    from src.utils.submission_utils import _create_ship_with_accession_tag

    mock_session = MagicMock()
    mock_accession = MagicMock()
    mock_accession.id = 5
    mock_accession.accession_tag = "SSA000042"
    mock_ship = MagicMock()
    mock_ship.id = 99

    mock_session.query.return_value.filter.return_value.first.side_effect = [
        mock_accession,
        mock_ship,
    ]

    with patch("src.utils.submission_utils.sequence_matches_ship", return_value=True):
        ship, accession, ship_id, tag = _create_ship_with_accession_tag(
            mock_session, "ATGCATGCATGCATGC", "SSA000042"
        )

    assert ship_id == 99
    assert tag == "SSA000042"
    mock_session.add.assert_not_called()

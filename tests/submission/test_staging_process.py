"""Tests for admin staging submission processing helpers."""

from unittest.mock import MagicMock

from src.utils.web_submission_adapter import (
    _classification_from_workflow_result,
    _submission_has_classification,
    summarize_staging_process_results,
)


def test_submission_has_classification():
    assert _submission_has_classification({"classification_family": "Prometheus"})
    assert not _submission_has_classification({"classification_family": ""})
    assert not _submission_has_classification({})


def test_classification_from_workflow_result():
    wf = {
        "complete": True,
        "found_match": True,
        "match_stage": "family",
        "classification_data": {
            "source": "family",
            "family": "Prometheus",
            "confidence": "Medium",
        },
    }
    cls = _classification_from_workflow_result(wf)
    assert cls is not None
    assert cls["family"] == "Prometheus"
    assert cls["source"] == "family"


def test_summarize_staging_process_results():
    summary = summarize_staging_process_results(
        [
            {"sub_id": 1, "success": True, "already_processed": True},
            {
                "sub_id": 2,
                "success": True,
                "accession": "SSA001.1",
                "classified": True,
                "classification": {"family": "Arwing"},
            },
            {
                "sub_id": 3,
                "success": True,
                "accession": "SSA002.1",
                "skipped_accession": True,
                "skipped_classification": True,
            },
        ]
    )
    assert "#1: skipped" in summary
    assert "Arwing" in summary
    assert "unchanged" in summary


def test_parse_gff_contents_plain_text():
    from src.utils.web_submission_adapter import _parse_gff_contents

    gff = "##gff-version 3\nchr1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n"
    entries = _parse_gff_contents(gff, "test.gff")
    assert len(entries) == 1
    assert entries[0]["seqid"] == "chr1"
    assert entries[0]["start"] == 1


def test_add_gff_entries_maps_to_gff_model():
    from src.utils.submission_utils import SubmissionProcessor

    processor = SubmissionProcessor()
    mock_session = MagicMock()
    entries = [
        {
            "seqid": "chr1",
            "source": "starfish",
            "type": "gene",
            "start": 1,
            "end": 100,
            "score": ".",
            "strand": "+",
            "phase": "0",
            "attributes": "ID=g1",
        }
    ]

    processor._add_gff_entries(mock_session, entries, accession_id=1, ship_id=2)

    added = mock_session.add.call_args[0][0]
    assert added.__class__.__name__ == "Gff"
    assert added.source == "starfish"
    assert added.phase == 0
    assert added.attributes == "seqid=chr1;ID=g1"
    assert added.accession_id == 1
    assert added.ship_id == 2

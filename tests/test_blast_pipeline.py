"""Unit tests for BLAST/classification pipeline helpers."""

from unittest.mock import patch

from src.utils.blast_utils import resolve_blast_result_status

with patch("dash.register_page"):
    from src.pages.blast import legacy_per_sequence_result


class TestResolveBlastResultStatus:
    def test_failed_when_error_present(self):
        status, error = resolve_blast_result_status(
            {"error": "BLAST timed out"}, top_level_error=None
        )
        assert status == "failed"
        assert error == "BLAST timed out"

    def test_failed_uses_top_level_error(self):
        status, error = resolve_blast_result_status(
            {}, top_level_error="Pipeline error"
        )
        assert status == "failed"
        assert error == "Pipeline error"

    def test_success_when_blast_content_present(self):
        status, error = resolve_blast_result_status(
            {"blast_content": "<xml/>", "processed": True}
        )
        assert status == "success"
        assert error is None

    def test_no_hits_when_processed_without_content(self):
        status, error = resolve_blast_result_status({"processed": True})
        assert status == "no_hits"
        assert error is None

    def test_loading_when_not_processed(self):
        status, error = resolve_blast_result_status({})
        assert status == "loading"
        assert error is None


class TestLegacyPerSequenceResult:
    def test_extracts_inner_sequence_result(self):
        converted = {
            "blast_content": "xml",
            "sequence_results": {
                "0": {
                    "blast_content": "xml",
                    "error": None,
                    "processed": True,
                }
            },
        }
        result = legacy_per_sequence_result(converted)
        assert result["blast_content"] == "xml"
        assert result["processed"] is True

    def test_falls_back_to_top_level_dict(self):
        converted = {"blast_content": "xml", "error": "fail"}
        result = legacy_per_sequence_result(converted)
        assert result["blast_content"] == "xml"
        assert result["error"] == "fail"

    def test_empty_input(self):
        assert legacy_per_sequence_result(None) == {}

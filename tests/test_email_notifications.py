"""Tests for submission email notifications."""

from unittest.mock import patch

import pytest

from src.utils import email_notifications as email_mod
from src.utils.email_notifications import (
    SubmissionEmailResult,
    email_configuration_status,
    notify_submission_received,
)


@pytest.fixture
def enabled_email_config(monkeypatch):
    monkeypatch.setattr(email_mod, "EMAIL_ENABLED", True)
    monkeypatch.setattr(email_mod, "SMTP_FROM_EMAIL", "noreply@test.org")
    monkeypatch.setattr(email_mod, "CURATOR_EMAILS", ["curator@test.org"])
    monkeypatch.setattr(email_mod, "SMTP_HOST", "localhost")
    monkeypatch.setattr(email_mod, "SMTP_PORT", 1025)
    monkeypatch.setattr(email_mod, "SMTP_USE_TLS", False)
    monkeypatch.setattr(email_mod, "SMTP_USER", "")
    monkeypatch.setattr(email_mod, "SMTP_PASSWORD", "")


def test_email_configuration_status_disabled(monkeypatch):
    monkeypatch.setattr(email_mod, "EMAIL_ENABLED", False)
    status = email_configuration_status()
    assert status["enabled"] is False
    assert status["configured"] is False
    assert "EMAIL_ENABLED" in status["missing"]


def test_notify_submission_skipped_when_disabled(monkeypatch):
    monkeypatch.setattr(email_mod, "EMAIL_ENABLED", False)
    result = notify_submission_received(
        "group-1",
        [{"seq_id": "s1", "genus": "A", "species": "a"}],
        uploader="user@test.org",
    )
    assert result.skipped is True
    assert result.curator_sent is False
    assert result.confirmation_sent is False


def test_notify_submission_sends_both(enabled_email_config):
    entries = [
        {
            "seq_id": "ship1",
            "genus": "Fusarium",
            "species": "oxysporum",
            "hostchr": "chr1",
            "shipstart": 100,
            "shipend": 2000,
        }
    ]
    with patch.object(email_mod, "_send_email") as mock_send:
        result = notify_submission_received(
            "group-abc",
            entries,
            uploader="submitter@test.org",
            evidence="manual",
            comment="test",
        )

    assert isinstance(result, SubmissionEmailResult)
    assert result.curator_sent is True
    assert result.confirmation_sent is True
    assert mock_send.call_count == 2


def test_notify_submission_records_smtp_failure(enabled_email_config):
    with patch.object(email_mod, "_send_email", side_effect=OSError("SMTP down")):
        result = notify_submission_received(
            "group-abc",
            [{"seq_id": "s1"}],
            uploader="submitter@test.org",
        )

    assert result.curator_sent is False
    assert result.confirmation_sent is False
    assert result.errors

"""Tests for submission email notifications."""

from unittest.mock import MagicMock, patch

import pytest

from src.utils import email_notifications as email_mod
from src.utils.email_notifications import (
    _format_ship_summary,
    send_curator_notification,
    send_submission_confirmation,
)


@pytest.fixture
def smtp_config(monkeypatch):
    monkeypatch.setattr(email_mod, "SMTP_HOST", "localhost")
    monkeypatch.setattr(email_mod, "SMTP_PORT", 25)
    monkeypatch.setattr(email_mod, "SMTP_USER", "smtp-user")
    monkeypatch.setattr(email_mod, "SMTP_PASSWORD", "smtp-pass")
    monkeypatch.setattr(email_mod, "SMTP_FROM_EMAIL", "noreply@test.org")
    monkeypatch.setattr(email_mod, "CURATOR_EMAILS", ["curator@test.org"])


@pytest.fixture
def mock_smtp():
    server = MagicMock()
    server.__enter__.return_value = server
    server.__exit__.return_value = False
    with patch.object(email_mod.smtplib, "SMTP", return_value=server):
        yield server


def test_curator_notification_skips_when_unconfigured(monkeypatch):
    monkeypatch.setattr(email_mod, "SMTP_USER", "")
    monkeypatch.setattr(email_mod, "CURATOR_EMAILS", [])

    with patch.object(email_mod.smtplib, "SMTP") as smtp:
        result = send_curator_notification(
            "group-1", [{"seq_id": "s1"}], uploader="user@test.org"
        )

    assert result is False
    smtp.assert_not_called()


def test_curator_notification_skips_without_curators(monkeypatch):
    monkeypatch.setattr(email_mod, "SMTP_USER", "smtp-user")
    monkeypatch.setattr(email_mod, "CURATOR_EMAILS", [""])

    with patch.object(email_mod.smtplib, "SMTP") as smtp:
        result = send_curator_notification(
            "group-1", [{"seq_id": "s1"}], uploader="user@test.org"
        )

    assert result is False
    smtp.assert_not_called()


def test_curator_notification_sends_when_configured(smtp_config, mock_smtp):
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

    result = send_curator_notification(
        "group-abc",
        entries,
        uploader="submitter@test.org",
        evidence="manual",
        comment="test",
    )

    assert result is True
    mock_smtp.sendmail.assert_called_once()
    from_addr, to_addrs, msg = mock_smtp.sendmail.call_args[0]
    assert from_addr == "noreply@test.org"
    assert to_addrs == ["curator@test.org"]
    assert "group-abc" in msg
    assert "ship1" in msg
    assert "submitter@test.org" in msg


def test_curator_notification_records_smtp_failure(smtp_config):
    with patch.object(email_mod.smtplib, "SMTP", side_effect=OSError("SMTP down")):
        result = send_curator_notification(
            "group-abc",
            [{"seq_id": "s1"}],
            uploader="submitter@test.org",
        )

    assert result is False


def test_confirmation_skips_when_unconfigured(monkeypatch):
    monkeypatch.setattr(email_mod, "SMTP_USER", "")

    with patch.object(email_mod.smtplib, "SMTP") as smtp:
        result = send_submission_confirmation("submitter@test.org", "group-1", 1)

    assert result is False
    smtp.assert_not_called()


def test_confirmation_sends_when_configured(smtp_config, mock_smtp):
    result = send_submission_confirmation("submitter@test.org", "group-abc", 3)

    assert result is True
    mock_smtp.sendmail.assert_called_once()
    from_addr, to_addrs, msg = mock_smtp.sendmail.call_args[0]
    assert from_addr == "noreply@test.org"
    assert to_addrs == ["submitter@test.org"]
    assert "group-abc" in msg
    assert "3 Starships" in msg


def test_confirmation_records_smtp_failure(smtp_config):
    with patch.object(email_mod.smtplib, "SMTP", side_effect=OSError("SMTP down")):
        result = send_submission_confirmation("submitter@test.org", "group-abc", 1)

    assert result is False


def test_format_ship_summary():
    entries = [
        {
            "seq_id": "ship1",
            "genus": "Fusarium",
            "species": "oxysporum",
            "hostchr": "chr1",
            "shipstart": 100,
            "shipend": 2000,
        },
    ]

    summary = _format_ship_summary(entries)

    assert "ship1" in summary
    assert "Fusarium oxysporum" in summary
    assert "chr1" in summary
    assert "100 - 2000" in summary

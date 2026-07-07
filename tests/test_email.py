#!/usr/bin/env python3
"""Manual test script for email configuration."""

from dotenv import load_dotenv

load_dotenv()

from src.utils.email_notifications import (
    email_configuration_status,
    notify_submission_received,
)

test_entry = {
    "seq_id": "test_ship",
    "genus": "Fusarium",
    "species": "oxysporum",
    "hostchr": "chr1",
    "shipstart": 1000,
    "shipend": 20000,
    "evidence": "manual annotation",
    "comment": "Test submission",
}

print("Email configuration:", email_configuration_status())
result = notify_submission_received(
    submission_group_id="test-group-123",
    entries=[test_entry],
    uploader="test@example.com",
    evidence=test_entry["evidence"],
    comment=test_entry["comment"],
)

if result.curator_sent and result.confirmation_sent:
    print("✓ Curator and confirmation emails sent.")
elif result.skipped:
    print(f"✗ Skipped: {result.skip_reason}")
else:
    print(f"✗ Partial/failed: {result.errors}")

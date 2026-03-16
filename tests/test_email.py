#!/usr/bin/env python3
"""Test email configuration."""

from dotenv import load_dotenv
from src.utils.email_notifications import send_curator_notification

load_dotenv()  # Load .env file

# Test submission data
test_data = {
    "uploader": "test@example.com",
    "evidence": "manual annotation",
    "genus": "Fusarium",
    "species": "oxysporum",
    "hostchr": "chr1",
    "shipstart": 1000,
    "shipend": 20000,
    "strand_radio": 1,
    "seq_filename": "test.fasta",
    "comment": "Test submission",
}

success = send_curator_notification(
    submission_id="test-123", submission_data=test_data, accession="SSA000001"
)

if success:
    print("✓ Email sent successfully!")
else:
    print("✗ Email failed. Check logs for details.")

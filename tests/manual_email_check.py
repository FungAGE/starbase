#!/usr/bin/env python3
"""Manual test script for email configuration.

Run directly (not part of the automated test suite):
    python tests/manual_email_check.py
"""

from dotenv import load_dotenv

load_dotenv()

from src.utils import email_notifications as email_mod

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


def configuration_status():
    curators = [e for e in email_mod.CURATOR_EMAILS if e]
    return {
        "smtp_host": email_mod.SMTP_HOST,
        "smtp_port": email_mod.SMTP_PORT,
        "smtp_user_configured": bool(email_mod.SMTP_USER),
        "from_email": email_mod.SMTP_FROM_EMAIL,
        "curator_emails": curators,
        "configured": bool(email_mod.SMTP_USER) and bool(curators),
    }


def main():
    status = configuration_status()
    print("Email configuration:", status)

    if not status["configured"]:
        print("Skipped: email not configured (set SMTP_USER and CURATOR_EMAILS).")
        return

    curator_sent = email_mod.send_curator_notification(
        "test-group-123",
        [test_entry],
        uploader="test@example.com",
        evidence=test_entry["evidence"],
        comment=test_entry["comment"],
    )
    confirmation_sent = email_mod.send_submission_confirmation(
        "test@example.com", "test-group-123", 1
    )

    if curator_sent and confirmation_sent:
        print("Curator and confirmation emails sent.")
    else:
        print(f"Failed: curator={curator_sent}, confirmation={confirmation_sent}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Email notifications for Starbase submissions.
"""

import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, Any, Optional, List
from datetime import datetime

from src.config.logging import get_logger

logger = get_logger(__name__)


# Configuration - these should be set in environment variables
SMTP_HOST = os.getenv("SMTP_HOST", "localhost")
SMTP_PORT = int(os.getenv("SMTP_PORT", "587"))
SMTP_USER = os.getenv("SMTP_USER", "")
SMTP_PASSWORD = os.getenv("SMTP_PASSWORD", "")
SMTP_FROM_EMAIL = os.getenv("SMTP_FROM_EMAIL", "noreply@starbase.org")
CURATOR_EMAILS = os.getenv("CURATOR_EMAILS", "").split(",")


def _format_ship_summary(entries: List[Dict[str, Any]]) -> str:
    lines = []
    for index, entry in enumerate(entries, start=1):
        label = entry.get("seq_id") or entry.get("filename") or f"Ship {index}"
        organism = f"{entry.get('genus', '')} {entry.get('species', '')}".strip()
        host = entry.get("hostchr") or "N/A"
        coords = f"{entry.get('shipstart', 'N/A')} - {entry.get('shipend', 'N/A')}"
        lines.append(f"  {index}. {label} ({organism}) @ {host} {coords}")
    return "\n".join(lines)


def send_curator_notification(
    submission_group_id: str,
    entries: List[Dict[str, Any]],
    uploader: Optional[str] = None,
    evidence: Optional[str] = None,
    comment: Optional[str] = None,
) -> bool:
    """
    Send email notification to curators about a new grouped submission.

    Args:
        submission_group_id: Shared submission group identifier
        entries: List of ship entry dicts from the upload form
        uploader: Submitter email
        evidence: Annotation method
        comment: Optional user comment

    Returns:
        True if email sent successfully, False otherwise
    """
    if not SMTP_USER or not CURATOR_EMAILS or not CURATOR_EMAILS[0]:
        logger.warning("Email notifications not configured - skipping")
        return False

    first = entries[0] if entries else {}
    submission_data = {
        "uploader": uploader or first.get("uploader"),
        "evidence": evidence or first.get("evidence"),
        "comment": comment or first.get("comment"),
        "ship_count": len(entries),
        "ship_summary": _format_ship_summary(entries),
    }

    try:
        subject = f"New Starship Submission ({len(entries)} ship{'s' if len(entries) != 1 else ''}): {submission_group_id[:8]}"
        body_html = _build_submission_email_html(submission_group_id, submission_data)
        body_text = _build_submission_email_text(submission_group_id, submission_data)

        msg = MIMEMultipart("alternative")
        msg["Subject"] = subject
        msg["From"] = SMTP_FROM_EMAIL
        msg["To"] = ", ".join(CURATOR_EMAILS)

        msg.attach(MIMEText(body_text, "plain"))
        msg.attach(MIMEText(body_html, "html"))

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            if SMTP_PORT == 587:
                server.starttls()
            if SMTP_USER and SMTP_PASSWORD:
                server.login(SMTP_USER, SMTP_PASSWORD)
            server.sendmail(SMTP_FROM_EMAIL, CURATOR_EMAILS, msg.as_string())

        logger.info(
            f"Curator notification sent for submission group {submission_group_id}"
        )
        return True

    except Exception as e:
        logger.error(f"Failed to send curator notification: {str(e)}")
        return False


def _build_submission_email_html(
    submission_group_id: str, submission_data: Dict[str, Any]
) -> str:
    """Build HTML email body."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    ship_summary = submission_data.get("ship_summary", "").replace("\n", "<br>")

    html = f"""
    <html>
      <head>
        <style>
          body {{ font-family: Arial, sans-serif; line-height: 1.6; color: #333; }}
          .container {{ max-width: 600px; margin: 0 auto; padding: 20px; }}
          .header {{ background-color: #4A5568; color: white; padding: 20px; border-radius: 5px; }}
          .content {{ background-color: #f7fafc; padding: 20px; border-radius: 5px; margin-top: 20px; }}
          .field {{ margin-bottom: 15px; }}
          .label {{ font-weight: bold; color: #2D3748; }}
          .value {{ color: #4A5568; margin-left: 10px; }}
          .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #E2E8F0; color: #718096; font-size: 12px; }}
        </style>
      </head>
      <body>
        <div class="container">
          <div class="header">
            <h2>New Starship Submission Received</h2>
            <p>Submission ID: {submission_group_id}</p>
          </div>
          
          <div class="content">
            <div class="field">
              <span class="label">Submitted:</span>
              <span class="value">{timestamp}</span>
            </div>
            
            <div class="field">
              <span class="label">Ships in batch:</span>
              <span class="value">{submission_data.get("ship_count", 1)}</span>
            </div>
            
            <h3>Submission Details</h3>
            
            <div class="field">
              <span class="label">Submitter Email:</span>
              <span class="value">{submission_data.get("uploader", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Evidence/Method:</span>
              <span class="value">{submission_data.get("evidence", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Ships:</span>
              <div class="value">{ship_summary}</div>
            </div>
            
            {f'<div class="field"><span class="label">Comments:</span><span class="value">{submission_data.get("comment", "None")}</span></div>' if submission_data.get("comment") else ""}
          </div>
          
          <div class="footer">
            <p>This is an automated notification from the Starbase submission system.</p>
            <p>Please review this submission in the curator dashboard.</p>
          </div>
        </div>
      </body>
    </html>
    """
    return html


def _build_submission_email_text(
    submission_group_id: str, submission_data: Dict[str, Any]
) -> str:
    """Build plain text email body."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    text = f"""
NEW STARSHIP SUBMISSION RECEIVED

Submission ID: {submission_group_id}
Submitted: {timestamp}
Ships in batch: {submission_data.get("ship_count", 1)}

SUBMISSION DETAILS
------------------
Submitter Email: {submission_data.get("uploader", "N/A")}
Evidence/Method: {submission_data.get("evidence", "N/A")}
Ships:
{submission_data.get("ship_summary", "")}
{f"Comments: {submission_data.get('comment', 'None')}" if submission_data.get("comment") else ""}

---
This is an automated notification from the Starbase submission system.
Please review this submission in the curator dashboard.
    """
    return text.strip()


def send_submission_confirmation(
    recipient_email: str,
    submission_group_id: str,
    ship_count: int = 1,
) -> bool:
    """
    Send confirmation email to submitter.

    Args:
        recipient_email: Submitter's email address
        submission_group_id: Shared submission group identifier
        ship_count: Number of ships in the batch

    Returns:
        True if email sent successfully, False otherwise
    """
    if not SMTP_USER:
        logger.warning("Email notifications not configured - skipping confirmation")
        return False

    try:
        subject = "Starship Submission Received"
        ship_label = f"{ship_count} Starship{'s' if ship_count != 1 else ''}"

        body_html = f"""
        <html>
          <body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
            <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
              <h2 style="color: #4A5568;">Thank you for your submission!</h2>
              
              <p>Your {ship_label} have been successfully uploaded to the submission portal.</p>
              
              <div style="background-color: #f7fafc; padding: 15px; border-radius: 5px; margin: 20px 0;">
                <p><strong>Submission ID:</strong> {submission_group_id}</p>
                <p><strong>Ships uploaded:</strong> {ship_count}</p>
              </div>
              
              <p>Your submission will be reviewed by our curation team. Once approved, it will be included in the next database release.</p>
              
              <p>Thank you for contributing to the Starship community!</p>
              
              <hr style="border: none; border-top: 1px solid #E2E8F0; margin: 20px 0;">
              <p style="color: #718096; font-size: 12px;">
                This is an automated message from Starbase. Please do not reply to this email.
              </p>
            </div>
          </body>
        </html>
        """

        msg = MIMEMultipart("alternative")
        msg["Subject"] = subject
        msg["From"] = SMTP_FROM_EMAIL
        msg["To"] = recipient_email

        msg.attach(MIMEText(body_html, "html"))

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            if SMTP_PORT == 587:
                server.starttls()
            if SMTP_USER and SMTP_PASSWORD:
                server.login(SMTP_USER, SMTP_PASSWORD)
            server.sendmail(SMTP_FROM_EMAIL, [recipient_email], msg.as_string())

        logger.info(f"Confirmation email sent to {recipient_email}")
        return True

    except Exception as e:
        logger.error(f"Failed to send confirmation email: {str(e)}")
        return False

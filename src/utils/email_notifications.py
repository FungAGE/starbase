#!/usr/bin/env python3
"""
Email notifications for Starbase submissions.
"""

import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import Dict, Any, Optional
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


def send_curator_notification(
    submission_id: str,
    submission_data: Dict[str, Any],
    accession: Optional[str] = None,
) -> bool:
    """
    Send email notification to curators about new submission.

    Args:
        submission_id: Unique submission identifier
        submission_data: Dict containing submission details
        accession: Assigned accession number (if available)

    Returns:
        True if email sent successfully, False otherwise
    """
    # Skip if no SMTP configuration
    if not SMTP_USER or not CURATOR_EMAILS or not CURATOR_EMAILS[0]:
        logger.warning("Email notifications not configured - skipping")
        return False

    try:
        # Prepare email content
        subject = f"New Starship Submission: {submission_id[:8]}"

        # Build email body
        body_html = _build_submission_email_html(
            submission_id, submission_data, accession
        )
        body_text = _build_submission_email_text(
            submission_id, submission_data, accession
        )

        # Create message
        msg = MIMEMultipart("alternative")
        msg["Subject"] = subject
        msg["From"] = SMTP_FROM_EMAIL
        msg["To"] = ", ".join(CURATOR_EMAILS)

        # Attach both plain text and HTML versions
        part1 = MIMEText(body_text, "plain")
        part2 = MIMEText(body_html, "html")
        msg.attach(part1)
        msg.attach(part2)

        # Send email
        with smtplib.SMTP(SMTP_HOST, SMTP_PORT) as server:
            if SMTP_PORT == 587:
                server.starttls()
            if SMTP_USER and SMTP_PASSWORD:
                server.login(SMTP_USER, SMTP_PASSWORD)
            server.sendmail(SMTP_FROM_EMAIL, CURATOR_EMAILS, msg.as_string())

        logger.info(f"Curator notification sent for submission {submission_id}")
        return True

    except Exception as e:
        logger.error(f"Failed to send curator notification: {str(e)}")
        return False


def _build_submission_email_html(
    submission_id: str, submission_data: Dict[str, Any], accession: Optional[str]
) -> str:
    """Build HTML email body."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

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
          .accession {{ background-color: #48BB78; color: white; padding: 10px; border-radius: 5px; margin: 20px 0; }}
          .footer {{ margin-top: 20px; padding-top: 20px; border-top: 1px solid #E2E8F0; color: #718096; font-size: 12px; }}
        </style>
      </head>
      <body>
        <div class="container">
          <div class="header">
            <h2>🚀 New Starship Submission Received</h2>
            <p>Submission ID: {submission_id}</p>
          </div>
          
          <div class="content">
            <div class="field">
              <span class="label">Submitted:</span>
              <span class="value">{timestamp}</span>
            </div>
            
            {f'<div class="accession"><strong>Accession Assigned:</strong> {accession}</div>' if accession else ""}
            
            <h3>Submission Details</h3>
            
            <div class="field">
              <span class="label">Curator Email:</span>
              <span class="value">{submission_data.get("uploader", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Evidence/Method:</span>
              <span class="value">{submission_data.get("evidence", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Organism:</span>
              <span class="value">{submission_data.get("genus", "N/A")} {submission_data.get("species", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Host Contig:</span>
              <span class="value">{submission_data.get("hostchr", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Coordinates:</span>
              <span class="value">{submission_data.get("shipstart", "N/A")} - {submission_data.get("shipend", "N/A")}</span>
            </div>
            
            <div class="field">
              <span class="label">Strand:</span>
              <span class="value">{"+" if submission_data.get("strand_radio") == 1 else "-"}</span>
            </div>
            
            <div class="field">
              <span class="label">Sequence File:</span>
              <span class="value">{submission_data.get("seq_filename", "N/A")}</span>
            </div>
            
            {f'<div class="field"><span class="label">Annotation File:</span><span class="value">{submission_data.get("anno_filename", "N/A")}</span></div>' if submission_data.get("anno_filename") else ""}
            
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
    submission_id: str, submission_data: Dict[str, Any], accession: Optional[str]
) -> str:
    """Build plain text email body."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    text = f"""
NEW STARSHIP SUBMISSION RECEIVED

Submission ID: {submission_id}
Submitted: {timestamp}
{f"Accession Assigned: {accession}" if accession else ""}

SUBMISSION DETAILS
------------------
Curator Email: {submission_data.get("uploader", "N/A")}
Evidence/Method: {submission_data.get("evidence", "N/A")}
Organism: {submission_data.get("genus", "N/A")} {submission_data.get("species", "N/A")}
Host Contig: {submission_data.get("hostchr", "N/A")}
Coordinates: {submission_data.get("shipstart", "N/A")} - {submission_data.get("shipend", "N/A")}
Strand: {"+" if submission_data.get("strand_radio") == 1 else "-"}
Sequence File: {submission_data.get("seq_filename", "N/A")}
{f"Annotation File: {submission_data.get('anno_filename', 'N/A')}" if submission_data.get("anno_filename") else ""}
{f"Comments: {submission_data.get('comment', 'None')}" if submission_data.get("comment") else ""}

---
This is an automated notification from the Starbase submission system.
Please review this submission in the curator dashboard.
    """
    return text.strip()


def send_submission_confirmation(
    recipient_email: str, submission_id: str, accession: Optional[str] = None
) -> bool:
    """
    Send confirmation email to submitter.

    Args:
        recipient_email: Submitter's email address
        submission_id: Unique submission identifier
        accession: Assigned accession number (if available)

    Returns:
        True if email sent successfully, False otherwise
    """
    # Skip if no SMTP configuration
    if not SMTP_USER:
        logger.warning("Email notifications not configured - skipping confirmation")
        return False

    try:
        subject = "Starship Submission Received"

        body_html = f"""
        <html>
          <body style="font-family: Arial, sans-serif; line-height: 1.6; color: #333;">
            <div style="max-width: 600px; margin: 0 auto; padding: 20px;">
              <h2 style="color: #4A5568;">Thank you for your submission!</h2>
              
              <p>Your Starship sequence has been successfully submitted to the database.</p>
              
              <div style="background-color: #f7fafc; padding: 15px; border-radius: 5px; margin: 20px 0;">
                <p><strong>Submission ID:</strong> {submission_id}</p>
                {f"<p><strong>Accession Number:</strong> {accession}</p>" if accession else ""}
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

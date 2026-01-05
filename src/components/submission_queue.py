#!/usr/bin/env python3
"""
Public submission queue component for Starbase.

Shows all pending/open submissions with anonymized user information.
"""

import dash_mantine_components as dmc
from sqlalchemy import text
from typing import List, Dict, Any

from src.config.logging import get_logger
from src.database.sql_engine import get_submissions_session

logger = get_logger(__name__)


def anonymize_email(email: str) -> str:
    """
    Anonymize email address for public display.

    Args:
        email: Full email address

    Returns:
        Anonymized email (e.g., j***@example.com)
    """
    if not email or "@" not in email:
        return "Anonymous"

    local, domain = email.split("@", 1)
    if len(local) <= 1:
        anonymized_local = "*"
    else:
        anonymized_local = f"{local[0]}{'*' * (len(local) - 1)}"

    return f"{anonymized_local}@{domain}"


def get_pending_submissions() -> List[Dict[str, Any]]:
    """
    Get all pending submissions from database.

    Returns:
        List of pending submission dicts
    """
    try:
        session = get_submissions_session()

        # Query for submissions that need review or are pending
        query = text("""
            SELECT 
                id,
                seq_filename,
                seq_date,
                uploader,
                evidence,
                genus,
                species,
                accession_tag,
                needs_review,
                created_at
            FROM submissions
            WHERE needs_review = TRUE OR needs_review IS NULL
            ORDER BY created_at DESC
            LIMIT 50
        """)

        result = session.execute(query)
        submissions = []

        for row in result:
            submissions.append(
                {
                    "id": row.id,
                    "filename": row.seq_filename,
                    "date": row.seq_date,
                    "submitter": anonymize_email(row.uploader),
                    "evidence": row.evidence,
                    "genus": row.genus,
                    "species": row.species,
                    "accession": row.accession_tag,
                    "needs_review": row.needs_review,
                    "created_at": row.created_at,
                }
            )

        session.close()
        return submissions

    except Exception as e:
        logger.error(f"Error fetching pending submissions: {str(e)}")
        return []


def create_submission_queue_card(submission: Dict[str, Any]) -> dmc.Card:
    """
    Create a card for a single submission.

    Args:
        submission: Submission data dict

    Returns:
        Card component
    """
    # Format date
    if isinstance(submission["date"], str):
        date_str = submission["date"]
    else:
        date_str = (
            submission["date"].strftime("%Y-%m-%d") if submission["date"] else "Unknown"
        )

    return dmc.Card(
        children=[
            dmc.Group(
                [
                    dmc.Stack(
                        [
                            dmc.Group(
                                [
                                    dmc.Text(
                                        submission["filename"] or "Unknown",
                                        fw=600,
                                        size="md",
                                    ),
                                    dmc.Badge(
                                        submission["accession"] or "Pending",
                                        color="blue"
                                        if submission["accession"]
                                        else "gray",
                                        variant="light",
                                    ),
                                ],
                                gap="sm",
                            ),
                            dmc.Group(
                                [
                                    dmc.Text(
                                        f"🧬 {submission['genus']} {submission['species']}",
                                        size="sm",
                                        c="dimmed",
                                    ),
                                    dmc.Text("•", size="sm", c="dimmed"),
                                    dmc.Text(
                                        f"📧 {submission['submitter']}",
                                        size="sm",
                                        c="dimmed",
                                    ),
                                ],
                                gap="xs",
                            ),
                            dmc.Group(
                                [
                                    dmc.Text(
                                        f"Method: {submission['evidence']}",
                                        size="xs",
                                        c="dimmed",
                                    ),
                                    dmc.Text("•", size="xs", c="dimmed"),
                                    dmc.Text(
                                        f"Submitted: {date_str}",
                                        size="xs",
                                        c="dimmed",
                                    ),
                                ],
                                gap="xs",
                            ),
                        ],
                        gap="xs",
                    ),
                ],
            ),
        ],
        withBorder=True,
        p="md",
        radius="md",
        mb="sm",
        shadow="xs",
    )


def create_submission_queue(max_items: int = 20) -> dmc.Container:
    """
    Create submission queue component.

    Args:
        max_items: Maximum number of submissions to display

    Returns:
        Container with submission queue
    """
    submissions = get_pending_submissions()[:max_items]

    if not submissions:
        return dmc.Container(
            children=[
                dmc.Paper(
                    children=[
                        dmc.Title("📋 Pending Submissions", order=3, mb="md"),
                        dmc.Alert(
                            "No pending submissions at this time.",
                            title="All Clear!",
                            color="green",
                            variant="light",
                        ),
                    ],
                    p="xl",
                    radius="md",
                    withBorder=True,
                )
            ],
            size="lg",
        )

    # Create submission cards
    cards = [create_submission_queue_card(sub) for sub in submissions]

    return dmc.Container(
        children=[
            dmc.Paper(
                children=[
                    dmc.Group(
                        [
                            dmc.Title("📋 Pending Submissions", order=3),
                            dmc.Badge(
                                f"{len(submissions)} pending",
                                color="orange",
                                variant="light",
                                size="lg",
                            ),
                        ],
                        justify="space-between",
                        mb="md",
                    ),
                    dmc.Text(
                        "All submissions are reviewed by our curation team before being added to the public database.",
                        size="sm",
                        c="dimmed",
                        mb="lg",
                    ),
                    dmc.Stack(cards, gap="sm"),
                ],
                p="xl",
                radius="md",
                withBorder=True,
            )
        ],
        size="lg",
    )


def create_compact_submission_queue(max_items: int = 5) -> dmc.Card:
    """
    Create compact submission queue for dashboard.

    Args:
        max_items: Maximum number of submissions to show

    Returns:
        Compact card with recent submissions
    """
    submissions = get_pending_submissions()[:max_items]

    if not submissions:
        return dmc.Card(
            children=[
                dmc.Text("📋 Recent Submissions", fw=700, mb="xs"),
                dmc.Text("No pending submissions", size="sm", c="dimmed"),
            ],
            withBorder=True,
            p="md",
        )

    items = []
    for sub in submissions:
        items.append(
            dmc.Group(
                [
                    dmc.Stack(
                        [
                            dmc.Text(
                                sub["filename"] or "Unknown",
                                size="sm",
                                fw=500,
                                truncate=True,
                            ),
                            dmc.Text(
                                f"{sub['genus']} {sub['species']}",
                                size="xs",
                                c="dimmed",
                                truncate=True,
                            ),
                        ],
                        gap=0,
                        style={"flex": 1},
                    ),
                    dmc.Badge(
                        sub["accession"] or "Pending",
                        size="sm",
                        variant="light",
                    ),
                ],
                mb="sm",
            )
        )

    return dmc.Card(
        children=[
            dmc.Text("📋 Recent Submissions", fw=700, mb="md"),
            dmc.Stack(items, gap="xs"),
            dmc.Text(
                f"Showing {len(submissions)} of {len(get_pending_submissions())} pending",
                size="xs",
                c="dimmed",
                ta="center",
                mt="md",
            ),
        ],
        withBorder=True,
        p="md",
        shadow="sm",
    )

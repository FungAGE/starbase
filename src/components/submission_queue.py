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
    Get pending submission groups from database.

    Returns one row per submission group (most recent ship row per group).
    """
    try:
        with get_submissions_session() as session:
            query = text("""
                SELECT
                    COALESCE(submission_group_id, CAST(id AS TEXT)) AS submission_group_id,
                    MIN(id) AS id,
                    COUNT(*) AS ship_count,
                    MIN(seq_filename) AS seq_filename,
                    MIN(seq_date) AS seq_date,
                    MIN(uploader) AS uploader,
                    MIN(evidence) AS evidence,
                    MIN(genus) AS genus,
                    MIN(species) AS species,
                    MIN(processing_status) AS processing_status
                FROM submissions
                WHERE processing_status IN ('pending', 'processed')
                  AND (needs_review = TRUE OR needs_review IS NULL)
                GROUP BY COALESCE(submission_group_id, CAST(id AS TEXT))
                ORDER BY MIN(id) DESC
                LIMIT 50
            """)

            result = session.execute(query)
            submissions = []

            for row in result:
                group_id = row.submission_group_id or str(row.id)
                status = row.processing_status or "pending"
                submissions.append(
                    {
                        "id": row.id,
                        "group_id": group_id,
                        "ship_count": row.ship_count,
                        "filename": row.seq_filename,
                        "date": row.seq_date,
                        "submitter": anonymize_email(row.uploader),
                        "evidence": row.evidence,
                        "genus": row.genus,
                        "species": row.species,
                        "processing_status": status,
                        "created_at": row.id,
                    }
                )

            return submissions

    except Exception as e:
        logger.error(f"Error fetching pending submissions: {str(e)}")
        return []


def create_submission_queue_card(submission: Dict[str, Any]) -> dmc.Card:
    """
    Create a card for a single submission group.

    Args:
        submission: Submission data dict

    Returns:
        Card component
    """
    if isinstance(submission["date"], str):
        date_str = submission["date"]
    else:
        date_str = (
            submission["date"].strftime("%Y-%m-%d") if submission["date"] else "Unknown"
        )

    status = submission.get("processing_status") or "pending"
    badge_label = "Under review" if status == "processed" else "Pending review"
    badge_color = "blue" if status == "processed" else "gray"
    group_label = submission.get("group_id") or f"#{submission['id']}"
    ship_count = submission.get("ship_count") or 1

    return dmc.Card(
        children=[
            dmc.Group(
                [
                    dmc.Stack(
                        [
                            dmc.Group(
                                [
                                    dmc.Text(
                                        group_label[:13]
                                        + ("…" if len(group_label) > 13 else ""),
                                        fw=600,
                                        size="md",
                                    ),
                                    dmc.Badge(
                                        badge_label,
                                        color=badge_color,
                                        variant="light",
                                    ),
                                ],
                                gap="sm",
                            ),
                            dmc.Group(
                                [
                                    dmc.Text(
                                        f"{submission['genus']} {submission['species']}",
                                        size="sm",
                                        c="dimmed",
                                    ),
                                    dmc.Text("•", size="sm", c="dimmed"),
                                    dmc.Text(
                                        f"{submission['submitter']}",
                                        size="sm",
                                        c="dimmed",
                                    ),
                                ],
                                gap="xs",
                            ),
                            dmc.Group(
                                [
                                    dmc.Text(
                                        f"{ship_count} ship{'s' if ship_count != 1 else ''}",
                                        size="xs",
                                        c="dimmed",
                                    ),
                                    dmc.Text("•", size="xs", c="dimmed"),
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
                        dmc.Title("Pending Submissions", order=3, mb="md"),
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
                    mb="xl",
                    style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
                )
            ],
            size="lg",
        )

    cards = [create_submission_queue_card(sub) for sub in submissions]

    output = dmc.Container(
        children=[
            dmc.Paper(
                children=[
                    dmc.Group(
                        [
                            dmc.Title("Pending Submissions", order=3),
                            dmc.Badge(
                                f"{len(submissions)} pending",
                                color="var(--mantine-color-orange-6)",
                                variant="light",
                                size="lg",
                            ),
                        ],
                        justify="space-between",
                        mb="md",
                    ),
                    dmc.Stack(cards, gap="sm"),
                ],
                p="xl",
                radius="md",
                withBorder=True,
                style={"borderLeft": "4px solid var(--mantine-color-indigo-5)"},
            )
        ],
        size="lg",
    )

    return output


def create_submission_queue_banner(max_items: int = 5):
    """Compact pending queue for the top of the submit page."""
    all_pending = get_pending_submissions()
    submissions = all_pending[:max_items]
    total_pending = len(all_pending)

    if not submissions:
        return dmc.Paper(
            children=[
                dmc.Group(
                    [
                        dmc.Text("Pending Submissions", fw=600, size="sm"),
                        dmc.Badge(
                            "All clear", color="green", variant="light", size="sm"
                        ),
                    ],
                    justify="space-between",
                ),
            ],
            p="md",
            radius="md",
            withBorder=True,
            mb="md",
        )

    items = []
    for sub in submissions:
        status = sub.get("processing_status") or "pending"
        badge_label = "Review" if status == "processed" else "Pending"
        group_label = sub.get("group_id") or f"#{sub['id']}"
        if len(group_label) > 16:
            group_label = group_label[:16] + "…"
        ship_count = sub.get("ship_count") or 1
        items.append(
            dmc.Group(
                [
                    dmc.Text(
                        group_label, size="sm", fw=500, style={"minWidth": "100px"}
                    ),
                    dmc.Text(
                        f"{sub.get('genus', '')} {sub.get('species', '')}".strip()
                        or "Unknown",
                        size="sm",
                        c="dimmed",
                        style={"flex": 1},
                        truncate=True,
                    ),
                    dmc.Text(
                        f"{ship_count} ship{'s' if ship_count != 1 else ''}",
                        size="xs",
                        c="dimmed",
                    ),
                    dmc.Badge(badge_label, size="sm", variant="light"),
                ],
                gap="sm",
                wrap="nowrap",
                style={"width": "100%"},
            )
        )

    return dmc.Paper(
        children=[
            dmc.Group(
                [
                    dmc.Text("Pending Submissions", fw=600, size="sm"),
                    dmc.Badge(
                        f"{total_pending} pending",
                        color="var(--mantine-color-orange-6)",
                        variant="light",
                        size="sm",
                    ),
                ],
                justify="space-between",
                mb="sm",
            ),
            dmc.Stack(items, gap="xs"),
            dmc.Text(
                f"Showing {len(submissions)} most recent"
                + (f" of {total_pending}" if total_pending > len(submissions) else ""),
                size="xs",
                c="dimmed",
                ta="right",
                mt="xs",
            ),
        ],
        p="md",
        radius="md",
        withBorder=True,
        mb="md",
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
                dmc.Text("Recent Submissions", fw=700, mb="xs"),
                dmc.Text("No pending submissions", size="sm", c="dimmed"),
            ],
            withBorder=True,
            p="md",
        )

    items = []
    for sub in submissions:
        status = sub.get("processing_status") or "pending"
        badge_label = "Review" if status == "processed" else "Pending"
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
                        badge_label,
                        size="sm",
                        variant="light",
                    ),
                ],
                mb="sm",
            )
        )

    return dmc.Card(
        children=[
            dmc.Text("Recent Submissions", fw=700, mb="md"),
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

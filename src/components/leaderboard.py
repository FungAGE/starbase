#!/usr/bin/env python3
"""
Contributor leaderboard component for Starbase.

Shows top contributors based on number of submissions.
"""

import dash_mantine_components as dmc
from typing import List, Dict, Any
from sqlalchemy import text

from src.config.logging import get_logger
from src.database.sql_engine import get_submissions_session

logger = get_logger(__name__)


def get_top_contributors(limit: int = 10) -> List[Dict[str, Any]]:
    """
    Get top contributors from submissions database.

    Args:
        limit: Number of top contributors to return

    Returns:
        List of dicts with contributor info
    """
    try:
        session = get_submissions_session()

        # Query to get top contributors (anonymized)
        query = text("""
            SELECT 
                uploader,
                COUNT(*) as submission_count,
                MIN(seq_date) as first_submission,
                MAX(seq_date) as last_submission
            FROM submissions
            WHERE uploader IS NOT NULL
            GROUP BY uploader
            ORDER BY submission_count DESC
            LIMIT :limit
        """)

        result = session.execute(query, {"limit": limit})
        contributors = []

        for row in result:
            # Anonymize email (show first letter + domain)
            email = row.uploader
            if "@" in email:
                local, domain = email.split("@", 1)
                anonymized = f"{local[0]}***@{domain}"
            else:
                anonymized = "Anonymous"

            contributors.append(
                {
                    "anonymized_email": anonymized,
                    "submission_count": row.submission_count,
                    "first_submission": row.first_submission,
                    "last_submission": row.last_submission,
                }
            )

        session.close()
        return contributors

    except Exception as e:
        logger.error(f"Error fetching top contributors: {str(e)}")
        return []


def create_leaderboard(limit: int = 10, title: str = "Top Contributors") -> dmc.Paper:
    """
    Create leaderboard component.

    Args:
        limit: Number of contributors to show
        title: Title for the leaderboard

    Returns:
        Dash Mantine Component Paper with leaderboard
    """
    contributors = get_top_contributors(limit)

    if not contributors:
        return dmc.Paper(
            children=[
                dmc.Title(title, order=3, mb="md"),
                dmc.Text(
                    "No submissions yet. Be the first to contribute!",
                    c="dimmed",
                    ta="center",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
        )

    # Create leaderboard rows
    leaderboard_items = []

    for idx, contributor in enumerate(contributors, 1):
        # Medal emoji for top 3
        if idx == 1:
            medal = "🥇"
        elif idx == 2:
            medal = "🥈"
        elif idx == 3:
            medal = "🥉"
        else:
            medal = f"#{idx}"

        leaderboard_items.append(
            dmc.Group(
                [
                    dmc.Text(medal, size="xl", fw=700, style={"width": "50px"}),
                    dmc.Stack(
                        [
                            dmc.Text(
                                contributor["anonymized_email"],
                                fw=500,
                                size="sm",
                            ),
                            dmc.Text(
                                f"{contributor['submission_count']} submission{'s' if contributor['submission_count'] > 1 else ''}",
                                size="xs",
                                c="dimmed",
                            ),
                        ],
                        gap=0,
                        style={"flex": 1},
                    ),
                    dmc.Badge(
                        f"{contributor['submission_count']}",
                        color="blue",
                        variant="light",
                        size="lg",
                    ),
                ],
                justify="space-between",
                align="center",
                p="md",
                style={
                    "borderBottom": "1px solid #E2E8F0"
                    if idx < len(contributors)
                    else "none"
                },
            )
        )

    return dmc.Paper(
        children=[
            dmc.Group(
                [
                    dmc.Title(title, order=3),
                    dmc.Badge(
                        f"Total: {sum(c['submission_count'] for c in contributors)} submissions",
                        color="green",
                        variant="light",
                    ),
                ],
                justify="space-between",
                mb="md",
            ),
            dmc.Stack(leaderboard_items, gap=0),
            dmc.Text(
                "📊 Rankings based on number of submitted sequences",
                size="xs",
                c="dimmed",
                ta="center",
                mt="md",
            ),
        ],
        p="xl",
        radius="md",
        withBorder=True,
        shadow="sm",
    )


def create_compact_leaderboard(limit: int = 5) -> dmc.Card:
    """
    Create compact leaderboard for sidebar or smaller spaces.

    Args:
        limit: Number of contributors to show

    Returns:
        Compact leaderboard card
    """
    contributors = get_top_contributors(limit)

    if not contributors:
        return dmc.Card(
            children=[
                dmc.Text("🏆 Top Contributors", fw=700, mb="xs"),
                dmc.Text("No submissions yet", size="sm", c="dimmed"),
            ],
            withBorder=True,
            p="md",
        )

    items = []
    for idx, contributor in enumerate(contributors, 1):
        items.append(
            dmc.Group(
                [
                    dmc.Text(f"{idx}.", size="sm", c="dimmed"),
                    dmc.Text(
                        contributor["anonymized_email"],
                        size="sm",
                        style={"flex": 1},
                        truncate=True,
                    ),
                    dmc.Badge(
                        str(contributor["submission_count"]),
                        size="sm",
                        variant="light",
                    ),
                ],
                gap="xs",
                mb="xs",
            )
        )

    return dmc.Card(
        children=[
            dmc.Text("🏆 Top Contributors", fw=700, mb="md"),
            dmc.Stack(items, gap="xs"),
        ],
        withBorder=True,
        p="md",
        shadow="sm",
    )

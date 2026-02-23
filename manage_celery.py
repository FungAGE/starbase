#!/usr/bin/env python3
"""
Management script for Celery processes in the Starbase application.
Usage:
    python manage_celery.py worker     # Start only the worker
    python manage_celery.py beat       # Start only the beat scheduler
    python manage_celery.py both       # Start both worker and beat
"""

import os
import sys
import subprocess
import argparse

# Add the project root to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config.logging import get_logger

logger = get_logger(__name__)


def start_worker():
    """Start Celery worker."""
    logger.info("Starting Celery worker...")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "celery",
            "-A",
            "src.config.celery_config",
            "worker",
            "--loglevel=info",
            "--concurrency=2",
        ]
    )


def start_beat():
    """Start Celery beat scheduler."""
    logger.info("Starting Celery beat scheduler...")
    subprocess.run(
        [
            sys.executable,
            "-m",
            "celery",
            "-A",
            "src.config.celery_config",
            "beat",
            "--loglevel=info",
        ]
    )


def start_both():
    """Start both worker and beat scheduler using the start_celery.py script."""
    logger.info("Starting both Celery worker and beat scheduler...")
    subprocess.run([sys.executable, "start_celery.py"])


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(description="Manage Celery processes for Starbase")
    parser.add_argument(
        "command",
        choices=["worker", "beat", "both"],
        help="Which Celery process to start",
    )

    args = parser.parse_args()

    if args.command == "worker":
        start_worker()
    elif args.command == "beat":
        start_beat()
    elif args.command == "both":
        start_both()
    else:
        logger.error(f"Unknown command: {args.command}")
        sys.exit(1)


if __name__ == "__main__":
    main()

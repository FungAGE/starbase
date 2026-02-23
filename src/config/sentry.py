import logging
import os

from dotenv import load_dotenv
import sentry_sdk
from sentry_sdk.integrations.celery import CeleryIntegration
from sentry_sdk.integrations.flask import FlaskIntegration
from sentry_sdk.integrations.logging import LoggingIntegration

from src.config.logging import get_logger


load_dotenv()
SENTRY_DSN = os.getenv("SENTRY_DSN")

logger = get_logger(__name__)


def init_sentry():
    """
    Initialize Sentry
    Note:
        - The default ThreadingIntegration wraps Thread.start, which conflicts with
        - Flask-Limiter's in-memory storage timer and can raise:
          - RuntimeError: threads can only be started once when the same Timer is started multiple times.

        - Disable default integrations and specify only the needed ones explicitly (Flask, Celery, Logging) to avoid the conflict.
    """
    if not SENTRY_DSN:
        return

    sentry_logging = LoggingIntegration(
        level=logging.INFO,
        event_level=logging.ERROR,
    )

    sentry_sdk.init(
        dsn=SENTRY_DSN,
        default_integrations=False,
        integrations=[
            FlaskIntegration(),
            CeleryIntegration(),
            sentry_logging,
        ],
        send_default_pii=True,
        enable_logs=True,
        traces_sample_rate=1.0,
        profile_session_sample_rate=1.0,
        profile_lifecycle="trace",
    )

    logger.info("Sentry initialized")

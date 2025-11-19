from contextlib import contextmanager
from src.config.database import StarbaseSession, SubmissionsSession, TelemetrySession

from src.config.logging import get_logger

logger = get_logger(__name__)


@contextmanager
def get_starbase_session():
    """Context manager for starbase database sessions"""
    session = StarbaseSession()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        logger.error(f"Database error: {e}")
        raise
    finally:
        session.close()


@contextmanager
def get_submissions_session():
    """Context manager for submissions database sessions"""
    session = SubmissionsSession()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        logger.error(f"Database error: {e}")
        raise
    finally:
        session.close()


@contextmanager
def get_telemetry_session():
    """Context manager for telemetry database sessions"""
    session = TelemetrySession()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        logger.error(f"Database error: {e}")
        raise
    finally:
        session.close()

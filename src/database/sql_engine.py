import logging
from contextlib import contextmanager
from src.config.database import (
    StarbaseSession,
    SubmissionsSession,
    TelemetrySession
)

logger = logging.getLogger(__name__)

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

def check_database_connection():
    """Check if database connections are working"""
    try:
        with get_starbase_session() as session:
            session.execute("SELECT 1")
        return True
    except Exception as e:
        logger.error(f"Database connection check failed: {e}")
        return False

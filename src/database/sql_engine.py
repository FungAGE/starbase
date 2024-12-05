import os
import logging
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError
from src.config.database import StarbaseSession, SubmissionsSession
from contextlib import contextmanager
from src.config.database import starbase_engine, submissions_engine, telemetry_engine

logger = logging.getLogger(__name__)

sql_connected = False


def connect_to_database(db_path):
    if not os.path.exists(db_path):
        logger.error(f"Database file not found at {db_path}")
        return None, None

    connection_str = f"sqlite:///{db_path}"
    try:
        engine = create_engine(connection_str)
        with engine.connect():
            logger.info(f"Successfully connected to {db_path}")
            Session = sessionmaker(bind=engine)
            sql_connected = True
            return engine, Session, sql_connected
    except OperationalError as e:
        logger.error(f"Operational error connecting to {db_path}: {e}")
    except Exception as e:
        logger.exception(f"Unexpected error connecting to {db_path}: {e}")
    return None, None


# Connect to databases
starbase_path = (
    "src/database/db/starbase.sqlite"
    if os.path.exists("src/database/db/starbase.sqlite")
    else "db/starbase.sqlite"
)
submissions_path = (
    "src/database/db/submissions.sqlite"
    if os.path.exists("src/database/db/submissions.sqlite")
    else "db/submissions.sqlite"
)

telemetry_path = (
    "src/database/db/telemetry.sqlite"
    if os.path.exists("src/database/db/telemetry.sqlite")
    else "db/telemetry.sqlite"
)

db_dir = os.path.dirname(starbase_path)

if "starbase_engine" in globals():
    starbase_engine.dispose()  # Properly closes the connection pool
    del starbase_engine

if "submissions_engine" in globals():
    submissions_engine.dispose()  # Properly closes the connection pool
    del submissions_engine


if "telemetry_engine" in globals():
    telemetry_engine.dispose()  # Properly closes the connection pool
    del telemetry_engine

starbase_engine, starbase_session_factory, sql_connected = connect_to_database(
    starbase_path
)
submissions_engine, submissions_session_factory, submissions_connected = (
    connect_to_database(submissions_path)
)

telemetry_engine, telemetry_session_factory, telemetry_connected = connect_to_database(
    telemetry_path
)

@contextmanager
def get_starbase_session():
    """Context manager for starbase database sessions"""
    session = StarbaseSession()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
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
    except:
        session.rollback()
        raise
    finally:
        session.close()

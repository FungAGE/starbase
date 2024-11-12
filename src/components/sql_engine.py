import os
import logging
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError

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
    "src/data/db/starbase.sqlite"
    if os.path.exists("src/data/db/starbase.sqlite")
    else "db/starbase.sqlite"
)
submissions_path = (
    "src/data/db/submissions.sqlite"
    if os.path.exists("src/data/db/submissions.sqlite")
    else "db/submissions.sqlite"
)

if "starbase_engine" in globals():
    starbase_engine.dispose()  # Properly closes the connection pool
    del starbase_engine

if "submissions_engine" in globals():
    submissions_engine.dispose()  # Properly closes the connection pool
    del submissions_engine


starbase_engine, starbase_session_factory, sql_connected = connect_to_database(
    starbase_path
)
submissions_engine, submissions_session_factory, submissions_connected = (
    connect_to_database(submissions_path)
)

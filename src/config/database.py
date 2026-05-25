from src.config.settings import DB_PATHS, IS_DEV
from sqlalchemy import create_engine
from sqlalchemy.pool import QueuePool
from sqlalchemy.orm import sessionmaker
import sqlalchemy.exc

from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

DATABASE_URLS = {
    "starbase": f"sqlite:///{DB_PATHS['starbase']}",
    "submissions": f"sqlite:///{DB_PATHS['submissions']}",
    "telemetry": f"sqlite:///{DB_PATHS['telemetry']}",
}

# Create engines with proper connection pooling
# For SQLite, we need special handling due to threading restrictions
engines = {}
for name, url in DATABASE_URLS.items():
    engines[name] = create_engine(
        url,
        poolclass=QueuePool,
        pool_size=10,
        max_overflow=20,
        pool_timeout=30,
        pool_recycle=1800,
        echo=IS_DEV,
        connect_args={"check_same_thread": False},
    )

# Extract individual engines for backward compatibility
starbase_engine = engines["starbase"]
submissions_engine = engines["submissions"]
telemetry_engine = engines["telemetry"]

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionSession = sessionmaker(bind=submissions_engine)
TelemetrySession = sessionmaker(bind=telemetry_engine)


# Define a common retry decorator for database operations
def db_retry_decorator(additional_retry_exceptions=()):
    """
    Create a retry decorator for database operations
    Args:
        additional_retry_exceptions: Tuple of additional exceptions to retry on
    """
    retry_exceptions = (sqlalchemy.exc.OperationalError,) + additional_retry_exceptions

    return retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        retry=retry_if_exception_type(retry_exceptions),
        before_sleep=lambda retry_state: logger.warning(
            f"Retrying database operation after error: {retry_state.outcome.exception()}"
        ),
    )


def check_database_connection(database_name="starbase"):
    """Verify the submissions database is accessible and properly configured."""
    from sqlalchemy import text
    from src.database.sql_engine import (
        get_starbase_session,
        get_submissions_session,
        get_telemetry_session,
    )

    try:
        if database_name == "starbase":
            with get_starbase_session() as session:
                result = session.execute(
                    text("SELECT COUNT(*) FROM joined_ships")
                ).scalar()
                logger.info(f"Starbase database check passed. Current ships: {result}")
                return True
        elif database_name == "submissions":
            with get_submissions_session() as session:
                result = session.execute(
                    text("SELECT COUNT(*) FROM submissions")
                ).scalar()
                logger.info(
                    f"Submissions database check passed. Current submissions: {result}"
                )
                return True
        elif database_name == "telemetry":
            with get_telemetry_session() as session:
                result = session.execute(
                    text("SELECT COUNT(*) FROM telemetry")
                ).scalar()
                logger.info(
                    f"Telemetry database check passed. Current telemetry: {result}"
                )
                return True
    except Exception as e:
        logger.error(f"Database check failed: {str(e)}")
        return False

import logging
import warnings
from src.config.settings import IS_DEV

DEFAULT_LOG_LEVEL = logging.DEBUG if IS_DEV else logging.INFO

logging.basicConfig(
    level=DEFAULT_LOG_LEVEL,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
    ],
)

# Control ALL loggers by default instead of just the root logger
if not IS_DEV:
    for logger_name in logging.root.manager.loggerDict:
        logger = logging.getLogger(logger_name)
        if logger_name.startswith(
            (
                "uvicorn",
                "werkzeug",
                "urllib3",
                "requests",
                "dash",
                "sqlalchemy",
            )
        ):
            logger.setLevel(logging.WARNING)
        else:
            # For application loggers, use INFO
            logger.setLevel(logging.INFO)


def configure_third_party_loggers():
    """
    Configure third-party loggers to be less verbose in production
    """
    for logger_name in ["uvicorn", "uvicorn.error", "uvicorn.access"]:
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG if IS_DEV else logging.WARNING)

    for logger_name in [
        "sqlalchemy",
        "sqlalchemy.engine",
        "sqlalchemy.pool",
        "sqlalchemy.orm",
    ]:
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG if IS_DEV else logging.WARNING)
        if not IS_DEV:
            logger.propagate = False

    for logger_name in ["werkzeug", "urllib3", "requests", "dash", "matplotlib"]:
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.INFO if IS_DEV else logging.WARNING)


configure_third_party_loggers()

warnings.filterwarnings("ignore")


def get_logger(name: str, level=None):
    """Get a configured logger for a specific module.

    Args:
        name: Logger name
        level: Optional level override (defaults to environment-based level)

    Returns:
        A logger that uses the root logger's handlers to avoid duplication
    """
    logger = logging.getLogger(name)
    logger.setLevel(level or DEFAULT_LOG_LEVEL)

    # DON'T add handlers - let messages propagate to root logger
    # This prevents duplicate log messages
    # The root logger already has the appropriate handler configured

    return logger

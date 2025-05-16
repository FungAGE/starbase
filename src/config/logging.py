import logging
import warnings
import os
from celery.signals import setup_logging

ENV = os.getenv("ENVIRONMENT", "development")
IS_DEV = ENV.lower() == "development"
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
                "celery",
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


@setup_logging.connect
def config_loggers(*args, **kwargs):
    """Configure celery logging to use our custom logger"""
    return True


def configure_third_party_loggers():
    """
    Configure third-party loggers to be less verbose in production
    """
    for logger_name in ["uvicorn", "uvicorn.error", "uvicorn.access"]:
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG if IS_DEV else logging.WARNING)

    for logger_name in [
        "celery",
        "celery.task",
        "celery.worker",
        "celery.worker.strategy",
        "celery.app.trace",
    ]:
        logger = logging.getLogger(logger_name)
        logger.setLevel(logging.DEBUG if IS_DEV else logging.WARNING)
        if not IS_DEV:
            logger.propagate = False

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
    """
    logger = logging.getLogger(name)
    logger.setLevel(level or DEFAULT_LOG_LEVEL)

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level or DEFAULT_LOG_LEVEL)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    return logger

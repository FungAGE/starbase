import os
import logging

from src.config.settings import IS_DEV

# Check if Celery should be enabled (default to False for single-pod deployments)
CELERY_ENABLED = os.getenv("CELERY_ENABLED", "false").lower() == "true"

# Try to import Celery - if not available, disable it
CELERY_AVAILABLE = False
celery = None

try:
    from celery import Celery
    from celery.schedules import crontab
    CELERY_AVAILABLE = True
except ImportError:
    # Celery not installed - will run tasks directly
    CELERY_AVAILABLE = False
    CELERY_ENABLED = False
    logger = logging.getLogger(__name__)
    logger.info("Celery not available - tasks will run directly")

# Only proceed with Celery setup if it's available
if CELERY_AVAILABLE:
    # Redis configuration
    CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
    CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

    # Configure Celery loggers immediately - before app initialization
    if not IS_DEV:
        for logger_name in [
            "celery",
            "celery.task",
            "celery.worker",
            "celery.worker.strategy",
            "celery.app.trace",
        ]:
            logger = logging.getLogger(logger_name)
            logger.setLevel(logging.WARNING)

    # Initialize Celery
    celery = Celery("starships")
else:
    # Create a dummy crontab for when Celery isn't available
    def crontab(*args, **kwargs):
        """Dummy crontab function when Celery is not available"""
        return None
    
    CELERY_BROKER_URL = None
    CELERY_RESULT_BACKEND = None

# Configure Celery if available
if CELERY_AVAILABLE and celery:
    # Configure using settings
    celery.conf.update(
        broker_url=CELERY_BROKER_URL,
        result_backend=CELERY_RESULT_BACKEND,
        task_serializer="json",
        result_serializer="json",
        accept_content=["json"],
        imports=[
            "src.tasks",
            "src.telemetry.tasks",
        ],
        worker_prefetch_multiplier=1,  # Disable prefetching for more predictable behavior
        task_time_limit=300,  # 5 minute timeout
        task_soft_time_limit=90,  # Soft timeout of 1.5 minutes
        task_track_started=True,
        timezone="UTC",
        enable_utc=True,
        broker_connection_retry_on_startup=True,
        worker_max_tasks_per_child=1000,
    )

    # Configure logging
    celery.conf.worker_hijack_root_logger = False
    celery.conf.worker_log_color = IS_DEV

    # In production, only log warnings and errors for task events
    if not IS_DEV:
        # Disable events in production to reduce logging
        celery.conf.worker_send_task_events = False
        celery.conf.task_send_sent_event = False
        celery.conf.task_track_started = False
        celery.conf.task_store_errors_even_if_ignored = True

        # Configure log formats for less verbosity
        celery.conf.worker_task_log_format = (
            "%(asctime)s: %(levelname)s - %(task_name)s[%(task_id)s] - %(message)s"
        )
        celery.conf.worker_log_format = (
            "%(asctime)s: %(levelname)s - %(name)s - %(message)s"
        )

    # Configure periodic tasks
    celery.conf.beat_schedule = {
        "health-check-every-5min": {
            "task": "src.telemetry.routes.health",
            "schedule": crontab(hour=0),  # Every 1 hour at midnight
        },
        "update-ip-locations-daily": {
            "task": "src.telemetry.tasks.update_ip_locations_task",
            "schedule": crontab(hour=0),  # Every day at midnight
        },
    }

    # This ensures tasks are properly registered
    celery.autodiscover_tasks(["src.tasks", "src.telemetry"])

    # Final adjustment to ensure loggers are properly set
    # This needs to run after autodiscover_tasks
    if not IS_DEV:
        for logger_name in [
            "celery",
            "celery.task",
            "celery.worker",
            "celery.worker.strategy",
            "celery.app.trace",
        ]:
            logger = logging.getLogger(logger_name)
            logger.setLevel(logging.WARNING)


def run_task(task_func, *args, **kwargs):
    """
    Run a task either via Celery (if enabled) or directly (if disabled).
    
    For single-pod deployments, set CELERY_ENABLED=false to run tasks directly.
    For distributed deployments with Redis, set CELERY_ENABLED=true.
    If Celery fails, log error and run directly as fallback
    Args:
        task_func: The Celery task function to run
        *args: Positional arguments to pass to the task
        **kwargs: Keyword arguments to pass to the task
    
    Returns:
        AsyncResult if Celery is enabled, None otherwise
    """
    if CELERY_ENABLED:
        try:
            return task_func.apply_async(args=args, kwargs=kwargs)
        except Exception as e:
            logger = logging.getLogger(__name__)
            logger.warning(f"Celery task failed, running directly: {str(e)}")
            return task_func(*args, **kwargs)
    else:
        return task_func(*args, **kwargs)

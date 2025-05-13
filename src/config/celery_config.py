from celery import Celery
from celery.schedules import crontab
import os

# Redis configuration
CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", "redis://localhost:6379/0")
CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", "redis://localhost:6379/0")

# Initialize Celery
celery = Celery("starships")

# Configure using settings
celery.conf.update(
    broker_url=CELERY_BROKER_URL,
    result_backend=CELERY_RESULT_BACKEND,
    task_serializer="json",
    result_serializer="json",
    accept_content=["json"],
    imports=[
        "src.tasks",
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

# Configure periodic tasks
celery.conf.beat_schedule = {
    "update-ip-locations-hourly": {
        "task": "src.utils.telemetry.update_ip_locations",
        "schedule": crontab(minute=0),  # Every hour at minute 0
    },
    "refresh-telemetry-every-15min": {
        "task": "src.utils.telemetry.refresh_telemetry",
        "schedule": crontab(minute="*/15"),  # Every 15 minutes
    },
    "check-cache-status-every-5min": {
        "task": "src.utils.telemetry.check_cache_status",
        "schedule": crontab(minute="*/5"),  # Every 5 minutes
    },
}

# This ensures tasks are properly registered
celery.autodiscover_tasks(["src.tasks"])

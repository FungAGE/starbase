from celery import Celery
import os
from celery.schedules import crontab

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
        "src.tasks",  # Make sure this is included
    ],
    worker_prefetch_multiplier=1,  # Disable prefetching for more predictable behavior
    task_time_limit=300,  # 5 minute timeout
    task_soft_time_limit=90,  # Soft timeout of 1.5 minutes
    task_track_started=True,
    timezone="UTC",
    enable_utc=True,
    beat_schedule={
        "daily-cache-cleanup": {
            "task": "cleanup_cache",
            "schedule": crontab(hour=3, minute=0),  # Run at 3 AM UTC
            "args": (48,),  # Keep files newer than 48 hours
        },
    },
)

# This ensures tasks are properly registered
celery.autodiscover_tasks(["src.tasks"])

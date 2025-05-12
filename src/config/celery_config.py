from celery import Celery
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
    task_serializer="pickle",
    result_serializer="pickle",
    accept_content=["json", "pickle"],
    imports=[
        "src.tasks",  # Make sure this is included
    ],
    worker_prefetch_multiplier=1,  # Disable prefetching for more predictable behavior
    task_time_limit=300,  # 5 minute timeout
    task_soft_time_limit=90,  # Soft timeout of 1.5 minutes
    task_track_started=True,
    timezone="UTC",
    enable_utc=True,
)

# This ensures tasks are properly registered
celery.autodiscover_tasks(["src.tasks"])

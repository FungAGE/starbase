from celery import Celery
from src.config.settings import DATA_DIR
import os

# Redis configuration
CELERY_BROKER_URL = os.getenv('CELERY_BROKER_URL', 'redis://localhost:6379/0')
CELERY_RESULT_BACKEND = os.getenv('CELERY_RESULT_BACKEND', 'redis://localhost:6379/0')

# Initialize Celery
celery = Celery(
    'starbase',
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND,
    include=['src.tasks']  # This will include all tasks defined in src/tasks/
)

# Celery configuration
celery.conf.update(
    worker_prefetch_multiplier=1,  # Disable prefetching for more predictable behavior
    task_time_limit=3600,  # 1 hour timeout
    task_soft_time_limit=3000,  # Soft timeout of 50 minutes
    task_track_started=True,
    task_serializer='json',
    result_serializer='json',
    accept_content=['json'],
    timezone='UTC',
    enable_utc=True,
) 
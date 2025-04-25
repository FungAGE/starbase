#!/bin/bash

# Start Redis server in the background if not using external Redis
redis-server --daemonize yes

# Start cron in the background using supercronic
supercronic $HOME/cron/crontab &

restart_celery() {
    # Kill existing Celery worker if it exists
    if [ -f /tmp/celery.pid ]; then
        kill $(cat /tmp/celery.pid) 2>/dev/null || true
        rm /tmp/celery.pid
    fi

    # Start new Celery worker
    celery -A src.config.celery_config:celery worker \
        --loglevel=DEBUG \
        --concurrency=4 \
        --max-tasks-per-child=1000 \
        --max-memory-per-child=1024000 \
        --pidfile=/tmp/celery.pid \
        &
}

# Start Celery initially
restart_celery

# Wait a moment for Redis and Celery to initialize
sleep 3

# Start the application with optimized settings
gunicorn --bind=0.0.0.0:8000 \
    --reload \
    --workers=4 \
    --threads=4 \
    --worker-class=gthread \
    --worker-tmp-dir=/dev/shm \
    --timeout=300 \
    --graceful-timeout=60 \
    --keep-alive=5 \
    --max-requests=1000 \
    --max-requests-jitter=50 \
    --forwarded-allow-ips='*' \
    --access-logfile - \
    app:server
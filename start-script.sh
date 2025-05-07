#!/bin/bash
set -e  # Exit on error

# Ensure script is executable
if [ ! -x "$0" ]; then
    echo "Error: Script is not executable"
    exit 1
fi

# Create and set permissions for cache directory
mkdir -p $HOME/src/database/db/cache/tmp
chmod -R 777 $HOME/src/database/db/cache

# Start Redis server in the background if not using external Redis
redis-server --daemonize yes

# Start cron in the background using supercronic
supercronic $HOME/cron/crontab &

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --dev)
            DEV_MODE=true
            shift # Remove --dev from processing
            ;;
        *)
            # Unknown option
            echo "Unknown option: $1"
            echo "Usage: $0 [--dev]"
            exit 1
            ;;
    esac
done

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

restart_celery_beat() {
    # Kill existing Celery beat if it exists
    if [ -f /tmp/celerybeat.pid ]; then
        kill $(cat /tmp/celerybeat.pid) 2>/dev/null || true
        rm /tmp/celerybeat.pid
    fi

    # Start new Celery beat scheduler
    celery -A src.config.celery_config:celery beat \
        --loglevel=INFO \
        --pidfile=/tmp/celerybeat.pid \
        --schedule=/tmp/celerybeat-schedule \
        &
}

# Start Celery worker and beat
restart_celery
restart_celery_beat

# Wait a moment for Redis and Celery to initialize
sleep 3

# Check if --dev flag was provided
if [ "$DEV_MODE" = "true" ]; then
    # Development mode with reload
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --reload \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --access-log \
        app:server
else
    # Production mode with workers
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --workers=4 \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --access-log \
        app:server
fi

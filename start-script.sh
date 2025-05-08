#!/bin/bash
set -e  # Exit on error

# Activate conda environment if not already activated
if [[ -z "$CONDA_PREFIX" || "$CONDA_PREFIX" != *"/starbase" ]]; then
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate starbase
fi

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

# Ensure script is executable
if [ ! -x "$0" ]; then
    echo "Error: Script is not executable"
    exit 1
fi

# Ensure RAM disk directory exists and has proper permissions
mkdir -p /dev/shm/starbase_cache/tmp
mkdir -p /dev/shm/starbase_cache/celery
chmod 777 /dev/shm/starbase_cache

# Start Redis server in the background if not using external Redis
redis-server --daemonize yes

# Start cron in the background using supercronic
supercronic $HOME/cron/crontab &

restart_celery() {
    # Kill existing Celery worker if it exists
    if [ -f /dev/shm/starbase_cache/celery/celery.pid ]; then
        kill $(cat /dev/shm/starbase_cache/celery/celery.pid) 2>/dev/null || true
        rm /dev/shm/starbase_cache/celery/celery.pid
    fi

    # Start new Celery worker
    celery -A src.config.celery_config:celery worker \
        --loglevel=DEBUG \
        --concurrency=4 \
        --max-tasks-per-child=1000 \
        --max-memory-per-child=1024000 \
        --pidfile=/dev/shm/starbase_cache/celery/celery.pid \
        &
}

restart_celery_beat() {
    # Kill existing Celery beat if it exists
    if [ -f /dev/shm/starbase_cache/celery/celerybeat.pid ]; then
        kill $(cat /dev/shm/starbase_cache/celery/celerybeat.pid) 2>/dev/null || true
        rm /dev/shm/starbase_cache/celery/celerybeat.pid
    fi

    # Start new Celery beat scheduler
    celery -A src.config.celery_config:celery beat \
        --loglevel=INFO \
        --pidfile=/dev/shm/starbase_cache/celery/celerybeat.pid \
        --schedule=/dev/shm/starbase_cache/celery/celerybeat-schedule \
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
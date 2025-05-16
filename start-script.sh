#!/bin/bash

# Start Redis server in the background if not using external Redis
redis-server --daemonize yes

# Default to production mode
DEV_MODE=false

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

# Activate conda environment - with fallback if conda not configured
if [ "$DEV_MODE" = "false" ]; then
    if command -v conda &> /dev/null; then
        # If conda exists but not in standard location
        eval "$(conda shell.bash hook)"
        conda activate starbase
    else
        source /opt/conda/etc/profile.d/conda.sh
        conda activate starbase 
    fi
fi

# Set environment variables based on mode
export PYTHONPATH=$(pwd)
if [ "$DEV_MODE" = "true" ]; then
    export ENVIRONMENT="development"
else
    export ENVIRONMENT="production"
    export PYTHONUNBUFFERED=1
    export SQLALCHEMY_WARN_20=1
    export SQLALCHEMY_SILENCE_UBER_WARNING=1
fi

restart_celery() {
    # Kill existing Celery worker if it exists
    if [ -f /tmp/celery.pid ]; then
        kill $(cat /tmp/celery.pid) 2>/dev/null || true
        rm /tmp/celery.pid
    fi

    # Kill existing Celery beat if it exists
    if [ -f /tmp/celerybeat.pid ]; then
        kill $(cat /tmp/celerybeat.pid) 2>/dev/null || true
        rm /tmp/celerybeat.pid
    fi

    # Set log level based on environment
    if [ "$DEV_MODE" = "true" ]; then
        CELERY_LOG_LEVEL=DEBUG
    else
        CELERY_LOG_LEVEL=WARNING
    fi

    # This environment variable will silence internal Celery logging
    if [ "$DEV_MODE" != "true" ]; then
        export CELERY_WORKER_HIJACK_ROOT_LOGGER=False
        export CELERY_WORKER_REDIRECT_STDOUTS=False
        
        # These env vars can help silence Celery logs in production
        export CELERY_WORKER_TASK_LOG_FORMAT="%(message)s"
        export CELERY_WORKER_REDIRECT_STDOUTS_LEVEL=WARNING
    fi

    # Start new Celery worker
    celery -A src.config.celery_config:celery worker \
        --loglevel=$CELERY_LOG_LEVEL \
        --concurrency=4 \
        --max-tasks-per-child=1000 \
        --max-memory-per-child=1024000 \
        --pidfile=/tmp/celery.pid \
        &
    
    # Start Celery beat for scheduled tasks
    celery -A src.config.celery_config:celery beat \
        --loglevel=$CELERY_LOG_LEVEL \
        --pidfile=/tmp/celerybeat.pid \
        &
}

# Start Celery initially
restart_celery

# Wait a moment for Redis and Celery to initialize
sleep 3

# Check if --dev flag was provided
if [ "$DEV_MODE" = "true" ]; then
    # Development mode with reload
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --reload \
        --log-level=debug \
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
        --log-level=warning \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --no-access-log \
        app:server
fi

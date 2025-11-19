#!/bin/bash

# get DEV_MODE from `.env` file
DEV_MODE=$(grep DEV_MODE .env | cut -d '=' -f 2)

if [ -z "$DEV_MODE" ]; then
    DEV_MODE=false
fi


# Parse command line arguments
while [ $# -gt 0 ]; do
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
    # we are provided 4cpu and 8gb of memory within a single pod
    # we can use 4 workers to maximize the utilization of the resources

    # key parameters to consider:
    # --worker-class=uvicorn.workers.UvicornWorker: Explicitly specifies the worker class for consistency.
    # --worker-timeout=300: Kills workers that have been running for more than 5 minutes, preventing infinite hangs while allowing your "few minutes" task duration.
    # --max-requests-jitter=50: Adds randomness to the 1000 request limit to prevent all workers from restarting simultaneously.
    # --timeout-graceful-shutdown=30: Reduces graceful shutdown time from 60 to 30 seconds for faster recovery.

    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --workers=4 \
        --worker-timeout=300 \
        --max-requests=1000 \
        --max-requests-jitter=50 \
        --log-level=warning \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=30 \
        --no-access-log \
        app:server
fi

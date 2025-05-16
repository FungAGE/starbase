#!/bin/bash

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

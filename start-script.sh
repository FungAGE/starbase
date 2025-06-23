#!/bin/bash

# Starbase Multi-Service Startup Script
# Supports: web, worker, beat, all

# Default values
DEV_MODE=false
SERVICE_ROLE="web"

# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        --dev)
            DEV_MODE=true
            shift
            ;;
        --role)
            SERVICE_ROLE="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --dev          Run in development mode"
            echo "  --role ROLE    Service role: web, worker, beat, all (default: web)"
            echo "  --help         Show this help message"
            echo ""
            echo "Service Roles:"
            echo "  web     - Start the web application (Dash + Flask)"
            echo "  worker  - Start Celery worker for task processing" 
            echo "  beat    - Start Celery beat for scheduled tasks"
            echo "  all     - Start all services (single-pod mode)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "🚀 Starting Starbase service: $SERVICE_ROLE"

# Activate conda environment - with fallback if conda not configured
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate starbase
else
    source /opt/conda/etc/profile.d/conda.sh
    conda activate starbase 
fi

# Set common environment variables
export PYTHONPATH=$(pwd):$(pwd)/src
if [ "$DEV_MODE" = "true" ]; then
    export ENVIRONMENT="development"
    echo "🔧 Development mode enabled"
else
    export ENVIRONMENT="production"
    export PYTHONUNBUFFERED=1
    export SQLALCHEMY_WARN_20=1
    export SQLALCHEMY_SILENCE_UBER_WARNING=1
    echo "🏭 Production mode enabled"
fi

# Function to wait for Redis (if Celery is enabled)
wait_for_redis() {
    if [ -n "$CELERY_BROKER_URL" ]; then
        echo "⏳ Waiting for Redis to be available..."
        max_attempts=30
        attempt=1
        
        while ! redis-cli -u "$CELERY_BROKER_URL" ping > /dev/null 2>&1; do
            if [ $attempt -ge $max_attempts ]; then
                echo "❌ Redis not available after $max_attempts attempts"
                exit 1
            fi
            echo "Waiting for Redis... (attempt $attempt/$max_attempts)"
            sleep 2
            attempt=$((attempt + 1))
        done
        echo "✅ Redis is available"
    fi
}

# Function to start web application
start_web() {
    echo "🌐 Starting web application..."
    
    if [ "$DEV_MODE" = "true" ]; then
        # Development mode with reload
        exec uvicorn --host=0.0.0.0 \
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
        exec uvicorn --host=0.0.0.0 \
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
}

# Function to start Celery worker
start_worker() {
    echo "👷 Starting Celery worker..."
    wait_for_redis
    
    if [ "$DEV_MODE" = "true" ]; then
        exec celery -A src.config.celery_config.celery worker \
            --loglevel=debug \
            --concurrency=${CELERY_CONCURRENCY:-1} \
            --max-tasks-per-child=100 \
            --time-limit=600 \
            --soft-time-limit=540
    else
        exec celery -A src.config.celery_config.celery worker \
            --loglevel=info \
            --concurrency=${CELERY_CONCURRENCY:-2} \
            --max-tasks-per-child=1000 \
            --time-limit=300 \
            --soft-time-limit=240
    fi
}

# Function to start Celery beat
start_beat() {
    echo "⏰ Starting Celery beat scheduler..."
    wait_for_redis
    
    # Remove old beat schedule file
    rm -f celerybeat-schedule /tmp/celerybeat-schedule
    
    if [ "$DEV_MODE" = "true" ]; then
        exec celery -A src.config.celery_config.celery beat \
            --loglevel=debug \
            --schedule=/tmp/celerybeat-schedule
    else
        exec celery -A src.config.celery_config.celery beat \
            --loglevel=info \
            --schedule=/tmp/celerybeat-schedule
    fi
}

# Function to start all services (single-pod deployment)
start_all() {
    echo "🚀 Starting all services in single-pod mode..."
    
    # Check if Redis is already running, if not start it
    if ! redis-cli ping > /dev/null 2>&1; then
        echo "🔴 Starting Redis server..."
        redis-server --daemonize yes --appendonly yes --dir /var/lib/redis --logfile /var/log/redis/redis.log
        sleep 2
        
        # Set Celery environment variables to use local Redis
        export CELERY_BROKER_URL="redis://localhost:6379/0"
        export CELERY_RESULT_BACKEND="redis://localhost:6379/0"
    fi
    
    wait_for_redis
    
    # Start Celery worker in background
    echo "👷 Starting Celery worker in background..."
    celery -A src.config.celery_config.celery worker \
        --loglevel=debug \
        --concurrency=1 \
        --max-tasks-per-child=100 &
    WORKER_PID=$!
    
    # Start Celery beat in background  
    echo "⏰ Starting Celery beat in background..."
    rm -f celerybeat-schedule /tmp/celerybeat-schedule
    celery -A src.config.celery_config.celery beat \
        --loglevel=debug \
        --schedule=/tmp/celerybeat-schedule &
    BEAT_PID=$!
    
    # Cleanup function
    cleanup() {
        echo "🛑 Shutting down services..."
        kill $WORKER_PID $BEAT_PID 2>/dev/null
        wait $WORKER_PID $BEAT_PID 2>/dev/null
        exit 0
    }
    
    # Set up signal handlers
    trap cleanup SIGTERM SIGINT
    
    # Start web application in foreground
    echo "🌐 Starting web application..."
    start_web
}

# Route to appropriate service based on role
case $SERVICE_ROLE in
    web)
        start_web
        ;;
    worker)
        start_worker
        ;;
    beat)
        start_beat
        ;;
    all)
        start_all
        ;;
    *)
        echo "❌ Unknown service role: $SERVICE_ROLE"
        echo "Valid roles: web, worker, beat, all"
        exit 1
        ;;
esac

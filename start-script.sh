#!/bin/bash

# Start cron in the background using supercronic
supercronic $HOME/cron/crontab &

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
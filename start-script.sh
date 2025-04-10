#!/bin/bash

# Start cron in the background using supercronic
supercronic $HOME/cron/crontab &

# Start the application with uvicorn in WSGI mode
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
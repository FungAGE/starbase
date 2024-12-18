#!/bin/bash

# Start cron in the background using supercronic (a cron for containers)
supercronic $HOME/cron/crontab &

# Start the application
gunicorn --bind=0.0.0.0:8000 --reload --workers=8 --thread=2 --worker-class=gthread --forwarded-allow-ips='*' --access-logfile - app:server
#!/bin/bash

gunicorn --bind=0.0.0.0:8000 --workers=8 --thread=2 --worker-class=gthread --forwarded-allow-ips='*' --access-logfile - app:server
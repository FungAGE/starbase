#!/bin/bash
gunicorn --bind=0.0.0.0:8000 --workers=2 --thread=4 --worker-class=gthread --forwarded-allow-ips='*' --access-logfile - src.app:server
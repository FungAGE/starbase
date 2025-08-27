"""
Telemetry package for STARBASE.
This package handles user analytics, IP geolocation, and request tracking.
"""

# Note: Tasks are imported in celery_config.py to avoid circular imports
# The tasks module is auto-discovered by Celery
from src.telemetry.utils import (
    log_request,
    get_client_ip,
    is_development_ip,
    blast_limit_decorator,
)

__all__ = [
    "log_request",
    "get_client_ip",
    "is_development_ip",
    "blast_limit_decorator",
]

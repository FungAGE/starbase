"""
Telemetry package for STARBASE.
This package handles user analytics, IP geolocation, and request tracking.
"""

# Note: Tasks are imported in celery_config.py to avoid circular imports
# The tasks module is auto-discovered by Celery
from src.telemetry.utils import (
    get_client_ip,
    is_development_ip,
)

__all__ = [
    "get_client_ip",
    "is_development_ip",
]

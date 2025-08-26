"""
Telemetry tasks module.
Contains tasks related to telemetry (formerly Celery tasks).
"""

import requests
from src.config.settings import IPSTACK_API_KEY
from src.telemetry.utils import update_ip_locations as _update_ip_locations
from src.config.logging import get_logger
from src.config.celery_config import celery

logger = get_logger(__name__)


@celery.task(name="update_ip_locations")
def update_ip_locations(api_key=None):
    """
    Update locations for any new IPs in request_logs that aren't in ip_locations.
    This was previously run hourly via Celery Beat.
    """
    if api_key is None:
        api_key = IPSTACK_API_KEY
    return _update_ip_locations(api_key)


@celery.task(name="refresh_telemetry")
def refresh_telemetry():
    """
    Refresh telemetry data (formerly a Celery task).
    This was previously run every 15 minutes via Celery Beat.
    """
    try:
        response = requests.post("http://localhost:8000/api/telemetry/refresh")
        return {"status": "success", "response": response.status_code}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "error": str(e)}


@celery.task(name="check_cache_status")
def check_cache_status():
    """
    Check cache status and refresh if needed.
    This was previously run every 5 minutes via Celery Beat.
    """
    try:
        # First check if cache is healthy
        status_response = requests.get("http://localhost:8000/api/cache/status")
        if status_response.status_code != 200:
            # Cache needs refresh
            refresh_response = requests.post("http://localhost:8000/api/cache/refresh")
            return {"status": "refreshed", "response": refresh_response.status_code}
        return {"status": "healthy"}
    except Exception as e:
        logger.error(f"Error checking cache status: {str(e)}")
        return {"status": "error", "error": str(e)}

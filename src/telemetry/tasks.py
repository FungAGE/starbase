"""
Telemetry tasks module.
Contains tasks related to telemetry (formerly Celery tasks).
"""

import requests
from datetime import datetime
from src.config.settings import IPSTACK_API_KEY
from src.telemetry.utils import update_ip_locations as _update_ip_locations
from src.config.logging import get_logger
from src.config.celery_config import celery, CELERY_AVAILABLE

logger = get_logger(__name__)


def _log_request_impl(ip_address, endpoint):
    """
    Implementation of log request task.
    Logs request details to telemetry database.
    """
    from src.telemetry.utils import is_development_ip, page_mapping
    from src.config.database import TelemetrySession
    from sqlalchemy import text
    
    try:
        # Skip development IPs
        if is_development_ip(ip_address):
            logger.debug(f"Skipping telemetry for development IP: {ip_address}")
            return
        
        # Only log valid endpoints
        if endpoint not in page_mapping:
            logger.debug(f"Skipping telemetry for non-mapped endpoint: {endpoint}")
            return
        
        session = TelemetrySession()
        try:
            # Check if this IP+endpoint was logged in the last hour
            check_query = """
            SELECT COUNT(*) FROM request_logs
            WHERE ip_address = :ip 
            AND endpoint = :endpoint
            AND datetime(timestamp) >= datetime('now', '-1 hour')
            """
            recent_count = session.execute(
                text(check_query),
                {"ip": ip_address, "endpoint": endpoint}
            ).scalar() or 0
            
            # Only log if not already logged in the last hour
            if recent_count == 0:
                insert_query = """
                INSERT INTO request_logs (ip_address, endpoint, timestamp)
                VALUES (:ip, :endpoint, :timestamp)
                """
                session.execute(
                    text(insert_query),
                    {"ip": ip_address, "endpoint": endpoint, "timestamp": datetime.now()},
                )
                session.commit()
                logger.debug(f"Logged request from {ip_address} to {endpoint}")
        finally:
            session.close()
            
    except Exception as e:
        logger.error(f"Error in log_request_task: {str(e)}")
        raise


def _update_ip_locations_impl(api_key=None):
    """
    Implementation of update IP locations task.
    Update locations for any new IPs in request_logs that aren't in ip_locations.
    """
    try:
        if api_key is None:
            api_key = IPSTACK_API_KEY
        result = _update_ip_locations(api_key)
        logger.info("IP locations updated successfully")
        return {"status": "success", "result": result}
    except Exception as e:
        logger.error(f"Error updating IP locations: {str(e)}")
        return {"status": "error", "error": str(e)}


def _refresh_telemetry_impl():
    """
    Implementation of refresh telemetry task.
    Refresh telemetry data (formerly a Celery task).
    """
    try:
        response = requests.post("http://localhost:8000/api/telemetry/refresh", timeout=30)
        logger.info(f"Telemetry refresh completed with status: {response.status_code}")
        return {"status": "success", "response": response.status_code}
    except requests.exceptions.Timeout:
        logger.error("Telemetry refresh timed out")
        return {"status": "error", "error": "Request timeout"}
    except requests.exceptions.ConnectionError:
        logger.error("Telemetry refresh connection error")
        return {"status": "error", "error": "Connection error"}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "error": str(e)}


def _check_cache_status_impl():
    """
    Implementation of check cache status task.
    Check cache status and refresh if needed.
    """
    try:
        # First check if cache is healthy
        status_response = requests.get("http://localhost:8000/api/cache/status", timeout=10)
        if status_response.status_code != 200:
            # Cache needs refresh
            logger.info("Cache status unhealthy, refreshing...")
            refresh_response = requests.post("http://localhost:8000/api/cache/refresh", timeout=30)
            logger.info(f"Cache refresh completed with status: {refresh_response.status_code}")
            return {"status": "refreshed", "response": refresh_response.status_code}
        
        logger.info("Cache status healthy")
        return {"status": "healthy"}
    except requests.exceptions.Timeout:
        logger.error("Cache status check timed out")
        return {"status": "error", "error": "Request timeout"}
    except requests.exceptions.ConnectionError:
        logger.error("Cache status check connection error")
        return {"status": "error", "error": "Connection error"}
    except Exception as e:
        logger.error(f"Error checking cache status: {str(e)}")
        return {"status": "error", "error": str(e)}


# Create task wrappers - either Celery tasks or plain functions
if CELERY_AVAILABLE and celery:
    # Wrap with Celery decorators
    @celery.task(name="src.telemetry.tasks.log_request", 
                 ignore_result=True, 
                 max_retries=3,
                 rate_limit='100/m')
    def log_request_task(ip_address, endpoint):
        """Celery task wrapper for log_request"""
        return _log_request_impl(ip_address, endpoint)

    @celery.task(name="src.telemetry.tasks.update_ip_locations")
    def update_ip_locations(api_key=None):
        """Celery task wrapper for update_ip_locations"""
        return _update_ip_locations_impl(api_key)

    @celery.task(name="src.telemetry.tasks.refresh_telemetry")
    def refresh_telemetry():
        """Celery task wrapper for refresh_telemetry"""
        return _refresh_telemetry_impl()

    @celery.task(name="src.telemetry.tasks.check_cache_status")
    def check_cache_status():
        """Celery task wrapper for check_cache_status"""
        return _check_cache_status_impl()
else:
    # Create plain function aliases
    def log_request_task(ip_address, endpoint):
        """Direct function call (no Celery)"""
        return _log_request_impl(ip_address, endpoint)

    def update_ip_locations(api_key=None):
        """Direct function call (no Celery)"""
        return _update_ip_locations_impl(api_key)

    def refresh_telemetry():
        """Direct function call (no Celery)"""
        return _refresh_telemetry_impl()

    def check_cache_status():
        """Direct function call (no Celery)"""
        return _check_cache_status_impl()

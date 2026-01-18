"""
Telemetry tasks module.
Contains tasks related to telemetry (formerly Celery tasks).
"""

import requests
from datetime import datetime
from src.config.logging import get_logger
from src.config.celery_config import run_task

logger = get_logger(__name__)


def _log_request_impl(ip_address, endpoint):
    """
    Implementation of log request task.
    Logs request details to telemetry database.
    """
    from src.telemetry.utils import is_development_ip
    from src.config.settings import PAGE_MAPPING
    from src.database.sql_engine import get_telemetry_session
    from sqlalchemy import text
    
    try:
        # Skip development IPs
        if is_development_ip(ip_address):
            logger.debug(f"Skipping telemetry for development IP: {ip_address}")
            return
        
        # Only log valid endpoints
        if endpoint not in PAGE_MAPPING:
            logger.debug(f"Skipping telemetry for non-mapped endpoint: {endpoint}")
            return
        
        with get_telemetry_session() as session:
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
            
    except Exception as e:
        logger.error(f"Error in log_request_task: {str(e)}")
        raise

def _update_ip_locations_impl(api_key=None):
    """Update IP geolocation data"""

    from sqlalchemy import text
    from src.config.settings import IPSTACK_API_KEY
    from src.telemetry.utils import GeolocatorService
    from src.database.sql_engine import get_telemetry_session
    
    if api_key is None:
        api_key = IPSTACK_API_KEY
        
    geolocator = GeolocatorService()
    
    # Find IPs that need lookup
    query = """
    SELECT DISTINCT r.ip_address 
    FROM request_logs r 
    LEFT JOIN ip_locations i ON r.ip_address = i.ip_address 
    WHERE i.ip_address IS NULL 
    """
    # Insert or update the location data
    upsert_query = """
    INSERT INTO ip_locations (
        ip_address, latitude, longitude, city, country, 
        last_updated, lookup_attempted
    ) VALUES (
        :ip, :lat, :lon, :city, :country, 
        :timestamp, TRUE
    )
    ON CONFLICT(ip_address) DO UPDATE SET
        latitude = :lat,
        longitude = :lon,
        city = :city,
        country = :country,
        last_updated = :timestamp,
        lookup_attempted = TRUE
    """

    with get_telemetry_session() as session:
        try:
            new_ips = session.execute(text(query)).fetchall()
            
            for (ip,) in new_ips:
                if geolocator.is_private_ip(ip):
                    continue
                    
                location = geolocator.get_location(ip, api_key)
                if location.lat == 0 and location.lon == 0:
                    continue

                session.execute(
                    text(upsert_query),
                    {
                        "ip": location.ip,
                        "lat": location.lat,
                        "lon": location.lon,
                        "city": location.city or "Unknown",
                        "country": location.country or "Unknown",
                        "timestamp": datetime.now(),
                    },
                )

                session.commit()
        except Exception as e:
            logger.error(f"Error updating location for IP {ip}: {str(e)}")
            session.rollback()

    return {"status": "success", "ips_processed": len(new_ips)}


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

def log_request_task(ip_address, endpoint):
    return run_task(_log_request_impl, ip_address, endpoint)
    
def update_ip_locations_task():
    return run_task(_update_ip_locations_impl)

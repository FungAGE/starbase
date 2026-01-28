"""
Telemetry utilities module.
Contains core telemetry functions used throughout the application.
"""

from pathlib import Path
from dotenv import load_dotenv
from typing import Optional
from dataclasses import dataclass
import plotly.graph_objects as go
from flask import request
import time
import requests
from functools import lru_cache
import ipaddress
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from sqlalchemy import text
from src.config.settings import PAGE_MAPPING
from typing import Dict, Any
from src.database.sql_engine import get_telemetry_session
from src.telemetry.visualize import create_time_series_figure, create_endpoints_figure, create_map_figure
from src.config.logging import get_logger

logger = get_logger(__name__)

# Load environment variables from .env file
env_path = Path(".") / ".env"
load_dotenv(dotenv_path=env_path)

@dataclass
class LocationInfo:
    ip: str
    lat: float
    lon: float
    city: Optional[str] = None
    country: Optional[str] = None

    @classmethod
    def create_empty(cls, ip: str):
        return cls(ip=ip, lat=0.0, lon=0.0)


def initialize_ip_locations_table():
    """Create the ip_locations table if it doesn't exist."""

    IP_LOCATIONS_TABLE_SQL = """
    CREATE TABLE IF NOT EXISTS ip_locations (
        ip_address TEXT PRIMARY KEY,
        latitude REAL,
        longitude REAL,
        city TEXT,
        country TEXT,
        last_updated TIMESTAMP,
        lookup_attempted BOOLEAN DEFAULT FALSE
    )
    """
    with get_telemetry_session() as session:
        try:
            session.execute(text(IP_LOCATIONS_TABLE_SQL))
            session.commit()
        except Exception as e:
            logger.error(f"Error creating ip_locations table: {str(e)}")
            session.rollback()


class GeolocatorService:
    def __init__(self):
        self.cache_timeout = 3600  # 1 hour cache
        self.session = self._create_retry_session()

    def _create_retry_session(self):
        """Create a requests session with retry logic"""
        session = requests.Session()
        retries = Retry(
            total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504]
        )
        session.mount("http://", HTTPAdapter(max_retries=retries))
        session.mount("https://", HTTPAdapter(max_retries=retries))
        return session

    def is_private_ip(self, ip: str) -> bool:
        try:
            ip_obj = ipaddress.ip_address(ip)
            return (
                ip_obj.is_private
                or ip_obj.is_loopback
                or ip_obj.is_link_local
                or ip_obj.is_multicast
            )
        except ValueError:
            return True  # Consider invalid IPs as private

    @lru_cache(maxsize=1000)
    def get_location_from_ipapi(self, ip: str) -> Optional[LocationInfo]:
        """Primary method using ip-api.com (free tier, limited to 45 req/minute)."""
        try:
            response = self.session.get(f"http://ip-api.com/json/{ip}", timeout=5)
            if response.status_code == 200:
                data = response.json()
                if data.get("status") == "success":
                    return LocationInfo(
                        ip=ip,
                        lat=float(data.get("lat", 0)),
                        lon=float(data.get("lon", 0)),
                        city=data.get("city"),
                        country=data.get("country"),
                    )
            return None
        except Exception as e:
            logger.warning(f"ip-api.com lookup failed for IP {ip}: {str(e)}")
            return None

    @lru_cache(maxsize=1000)
    def get_location_from_ipstack(
        self, ip: str, api_key: Optional[str] = None
    ) -> Optional[LocationInfo]:
        """Fallback using ipstack.com (requires API key)."""
        if not api_key:
            return None

        try:
            response = self.session.get(
                f"http://api.ipstack.com/{ip}?access_key={api_key}", timeout=5
            )
            if response.status_code == 200:
                data = response.json()
                return LocationInfo(
                    ip=ip,
                    lat=float(data.get("latitude", 0)),
                    lon=float(data.get("longitude", 0)),
                    city=data.get("city"),
                    country=data.get("country_name"),
                )
            return None
        except Exception as e:
            logger.warning(f"ipstack lookup failed for IP {ip}: {str(e)}")
            return None

    @lru_cache(maxsize=1000)
    def get_location_from_ipwhois(self, ip: str) -> Optional[LocationInfo]:
        """Second fallback using ipwhois.app (free tier, 10k requests per month)."""
        try:
            response = self.session.get(f"http://ipwho.is/{ip}", timeout=5)
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    return LocationInfo(
                        ip=ip,
                        lat=float(data.get("latitude", 0)),
                        lon=float(data.get("longitude", 0)),
                        city=data.get("city"),
                        country=data.get("country"),
                    )
            return None
        except Exception as e:
            logger.warning(f"ipwhois lookup failed for IP {ip}: {str(e)}")
            return None

    def get_location(
        self, ip: str, ipstack_api_key: Optional[str] = None
    ) -> LocationInfo:
        """
        Get location information for an IP address using multiple fallback services.
        Returns LocationInfo with zeroed coordinates if all services fail.
        """
        # Skip private/invalid IPs
        if self.is_private_ip(ip):
            logger.debug(f"Skipping private/invalid IP: {ip}")
            return LocationInfo.create_empty(ip)

        # Try primary service
        location = self.get_location_from_ipapi(ip)
        if location:
            time.sleep(1.5)  # Rate limiting for ip-api.com
            return location

        # Try first fallback
        location = self.get_location_from_ipwhois(ip)
        if location:
            return location

        # Try second fallback if API key is provided
        if ipstack_api_key:
            location = self.get_location_from_ipstack(ip, ipstack_api_key)
            if location:
                return location

        # Return empty location if all methods fail
        logger.warning(f"All geolocation methods failed for IP: {ip}")
        return LocationInfo.create_empty(ip)


def get_client_ip():
    """Get the client's IP address from the request."""
    if request.headers.get("X-Forwarded-For"):
        # If behind a proxy, get real IP
        return request.headers.get("X-Forwarded-For").split(",")[0]
    return request.remote_addr


def is_development_ip(ip_address):
    """Check if IP is a development/local IP"""
    development_ips = {
        "127.0.0.1",  # localhost
        "0.0.0.0",  # all interfaces
        "::1",  # IPv6 localhost
    }

    # Check exact matches first
    if ip_address in development_ips:
        return True

    # Check for local network ranges
    local_prefixes = (
        "192.168.",
        "10.",
        "172.16.",
        "172.17.",
        "172.18.",
        "172.19.",
        "172.20.",
        "172.21.",
        "172.22.",
        "172.23.",
        "172.24.",
        "172.25.",
        "172.26.",
        "172.27.",
        "172.28.",
        "172.29.",
        "172.30.",
        "172.31.",
    )

    return any(ip_address.startswith(prefix) for prefix in local_prefixes)


def get_private_ip_filter_sql():
    """
    Returns SQL WHERE clause conditions to filter out private/development IPs.
    Use this in queries to exclude local network traffic.
    """
    return """
        ip_address NOT LIKE '192.168.%'
        AND ip_address NOT LIKE '10.%'
        AND ip_address NOT LIKE '172.1_%'
        AND ip_address NOT LIKE '172.2_%'
        AND ip_address NOT LIKE '172.3_%'
        AND ip_address != '127.0.0.1'
        AND ip_address != 'localhost'
        AND ip_address != '::1'
    """.strip()


def get_telemetry_data():
    """Fetch telemetry data from the database."""
    
    # Get reusable IP filter
    ip_filter = get_private_ip_filter_sql()

    with get_telemetry_session() as session:
        try:
            # Get unique users count, excluding private IPs
            unique_users_query = f"""
            SELECT COUNT(DISTINCT r.ip_address) as count 
            FROM request_logs r
            WHERE {ip_filter}
            """
            unique_users = session.execute(text(unique_users_query)).scalar() or 0

            # Get time series data for unique daily visitors
            time_series_query = f"""
            SELECT DATE(timestamp) as date, COUNT(DISTINCT ip_address) as count 
            FROM request_logs 
            WHERE endpoint IN ({PAGE_MAPPING})
            AND {ip_filter}
            GROUP BY DATE(timestamp) 
            ORDER BY date
            """
            time_series_data = session.execute(text(time_series_query)).fetchall()

            # Get unique visitors per endpoint
            endpoints_query = f"""
            SELECT endpoint, COUNT(DISTINCT ip_address) as count 
            FROM request_logs 
            WHERE endpoint IN ({PAGE_MAPPING})
            AND {ip_filter}
            GROUP BY endpoint 
            ORDER BY count DESC
            """
            endpoints_data = session.execute(text(endpoints_query)).fetchall()

            return {
                "unique_users": unique_users,
                "time_series_data": time_series_data,
                "endpoints_data": endpoints_data,
            }

        except Exception as e:
            logger.error(f"Error fetching telemetry data: {str(e)}")
            return None


def get_ip_locations():
    """Get locations from the ip_locations table."""
    with get_telemetry_session() as session:
        try:
            # Only get successfully looked up locations
            query = """
            SELECT ip_address, latitude, longitude, city, country
            FROM ip_locations
            WHERE lookup_attempted = TRUE
            AND latitude != 0 AND longitude != 0
            """
            results = session.execute(text(query)).fetchall()

            return [
                {
                    "ip": row[0],
                    "lat": row[1],
                    "lon": row[2],
                    "city": row[3] or "Unknown",
                    "country": row[4] or "Unknown",
                }
                for row in results
            ]

        except Exception as e:
            logger.error(f"Error fetching cached locations: {str(e)}")
            return []


def count_blast_submissions(ip_address, hours=1):
    """Count BLAST submissions from an IP address in the last N hours."""
    with get_telemetry_session() as session:
        try:
            query = """
            SELECT COUNT(*) as count 
            FROM request_logs 
            WHERE ip_address = :ip 
            AND endpoint = '/api/blast-submit'
            AND datetime(timestamp) >= datetime('now', :hours_ago)
            """
            result = (
                session.execute(
                    text(query), {"ip": ip_address, "hours_ago": f"-{hours} hours"}
                ).scalar()
                or 0
            )
            return result
        except Exception as e:
            logger.error(f"Error counting BLAST submissions: {str(e)}")
            return 0


    return fig


def get_telemetry_visualizations() -> Dict[str, Any]:
    """Analyze telemetry data and create visualizations."""
    try:
        # Fetch basic telemetry data
        telemetry_data = get_telemetry_data()
        if telemetry_data is None:
            return {
                "unique_users": 0,
                "total_requests": 0,
                "time_series": go.Figure(),
                "endpoints": go.Figure(),
                "map": go.Figure(),
                "locations": [],
            }

        # Create modern visualizations
        time_series_fig = create_time_series_figure(telemetry_data["time_series_data"])

        # Create endpoints visualization
        endpoints_fig = create_endpoints_figure(telemetry_data["endpoints_data"])

        # Get cached locations and create map
        locations = get_ip_locations()
        map_fig = create_map_figure(locations)

        # Calculate total requests
        total_requests = (
            sum(row[1] for row in telemetry_data["endpoints_data"] if row[0] in PAGE_MAPPING)
            if telemetry_data["endpoints_data"]
            else 0
        )

        return {
            "unique_users": telemetry_data["unique_users"],
            "total_requests": total_requests,
            "time_series": time_series_fig,
            "endpoints": endpoints_fig,
            "map": map_fig,
            "locations": locations,
        }

    except Exception as e:
        logger.error(f"Error in analyze_telemetry: {str(e)}")
        return {
            "unique_users": 0,
            "total_requests": 0,
            "time_series": go.Figure(),
            "endpoints": go.Figure(),
            "map": go.Figure(),
            "locations": [],
        }



def get_blast_limit_info(ip_address):
    """Get rate limit info for an IP address."""
    try:
        submissions = count_blast_submissions(ip_address)
        limit = 10  # Configurable limit per hour
        remaining = max(0, limit - submissions)
        return {"remaining": remaining, "limit": limit, "submissions": submissions}
    except Exception as e:
        logger.error(f"Error getting BLAST limit info: {str(e)}")
        return {"remaining": 0, "limit": 10, "submissions": 10}


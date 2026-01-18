"""
Telemetry utilities module.
Contains core telemetry functions used throughout the application.
"""

from pathlib import Path
from dotenv import load_dotenv
from datetime import datetime
from typing import Optional, Dict, Any
from dataclasses import dataclass
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from flask import request
import time
import requests
from functools import lru_cache
import ipaddress
from functools import wraps
from dash.exceptions import PreventUpdate
from requests.adapters import HTTPAdapter
from urllib3.util import Retry
from sqlalchemy import text
from src.database.sql_engine import get_telemetry_session

from src.config.logging import get_logger

logger = get_logger(__name__)

# Load environment variables from .env file
env_path = Path(".") / ".env"
load_dotenv(dotenv_path=env_path)

# Define valid pages
page_mapping = {
    "/": "Home",
    "/download": "Download",
    "/pgv": "PGV",
    "/submit": "Submit",
    "/blast": "BLAST",
    "/wiki": "Wiki",
    "/metrics": "Metrics",
    "/starfish": "Starfish",
    "/about": "About",
}

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


def update_ip_locations(api_key=None):
    """
    Update locations for any new IPs in request_logs that aren't in ip_locations.
    This should be run periodically (e.g., hourly) rather than on every app launch.
    """
    geolocator = GeolocatorService()

    with get_telemetry_session() as session:
        try:
            # Find IPs that haven't been looked up yet
            query = """
            SELECT DISTINCT r.ip_address 
            FROM request_logs r 
            LEFT JOIN ip_locations i ON r.ip_address = i.ip_address 
            WHERE i.ip_address IS NULL 
            OR (i.lookup_attempted = FALSE AND i.last_updated < datetime('now', '-7 days'))
            """
            new_ips = session.execute(text(query)).fetchall()

            for (ip,) in new_ips:
                # Skip private IPs entirely
                if geolocator.is_private_ip(ip):
                    logger.debug(f"Skipping private IP: {ip}")
                    continue

                try:
                    # Get location data
                    location = geolocator.get_location(ip, api_key)

                    # Skip if we got an empty location
                    if location.lat == 0 and location.lon == 0:
                        logger.debug(f"Skipping IP with no location data: {ip}")
                        continue

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

        except Exception as e:
            logger.error(f"Error in update_ip_locations: {str(e)}")
            return {"status": "error", "error": str(e)}


def get_ip_locations():
    """
    Get cached location data from the ip_locations table.
    This is now a fast database query instead of making API calls.
    """
    with get_telemetry_session() as session:
        try:
            # Only get locations that were successfully looked up
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
            logger.error(f"Error fetching IP locations: {str(e)}")
            return []


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


def log_request(ip_address, endpoint, session):
    """Log request details to telemetry database."""
    if is_development_ip(ip_address):
        logger.debug(f"Skipping telemetry for development IP: {ip_address}")
        return

    # Only log valid endpoints
    if endpoint not in page_mapping:
        logger.debug(f"Skipping telemetry for non-mapped endpoint: {endpoint}")
        return

    try:
        query = """
        INSERT INTO request_logs (ip_address, endpoint, timestamp)
        VALUES (:ip, :endpoint, :timestamp)
        """
        session.execute(
            text(query),
            {"ip": ip_address, "endpoint": endpoint, "timestamp": datetime.now()},
        )
        session.commit()
        logger.debug(f"Logged request from {ip_address} to {endpoint}")
    except Exception as e:
        session.rollback()
        logger.error(f"Error logging request: {str(e)}")


def fetch_telemetry_data():
    """Fetch telemetry data from the database."""

    # Define valid pages as a comma-separated string of quoted values
    valid_endpoints = "'/','/download','/pgv','/submit','/blast','/wiki','/metrics','/starfish','/about'"

    with get_telemetry_session() as session:
        try:
            # Get unique users count, excluding private IPs
            unique_users_query = """
            SELECT COUNT(DISTINCT r.ip_address) as count 
            FROM request_logs r
            WHERE r.ip_address NOT LIKE '192.168.%'
            AND r.ip_address NOT LIKE '10.%'
            AND r.ip_address NOT LIKE '172.1_%'
            AND r.ip_address NOT LIKE '172.2_%'
            AND r.ip_address NOT LIKE '172.3_%'
            AND r.ip_address != '127.0.0.1'
            AND r.ip_address != 'localhost'
            AND r.ip_address != '::1'
            """
            unique_users = session.execute(text(unique_users_query)).scalar() or 0

            # Get time series data for unique daily visitors
            time_series_query = f"""
            SELECT DATE(timestamp) as date, COUNT(DISTINCT ip_address) as count 
            FROM request_logs 
            WHERE endpoint IN ({valid_endpoints})
            AND ip_address NOT LIKE '192.168.%'
            AND ip_address NOT LIKE '10.%'
            AND ip_address NOT LIKE '172.1_%'
            AND ip_address NOT LIKE '172.2_%'
            AND ip_address NOT LIKE '172.3_%'
            AND ip_address != '127.0.0.1'
            AND ip_address != 'localhost'
            AND ip_address != '::1'
            GROUP BY DATE(timestamp) 
            ORDER BY date
            """
            time_series_data = session.execute(text(time_series_query)).fetchall()

            # Get unique visitors per endpoint
            endpoints_query = f"""
            SELECT endpoint, COUNT(DISTINCT ip_address) as count 
            FROM request_logs 
            WHERE endpoint IN ({valid_endpoints})
            AND ip_address NOT LIKE '192.168.%'
            AND ip_address NOT LIKE '10.%'
            AND ip_address NOT LIKE '172.1_%'
            AND ip_address NOT LIKE '172.2_%'
            AND ip_address NOT LIKE '172.3_%'
            AND ip_address != '127.0.0.1'
            AND ip_address != 'localhost'
            AND ip_address != '::1'
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


def get_cached_locations():
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


def create_modern_time_series_figure(time_series_data):
    """Create a modern time series visualization with better styling."""
    if time_series_data:
        df = pd.DataFrame(time_series_data, columns=["date", "count"])
        df["date"] = pd.to_datetime(df["date"])

        # Calculate moving average for trend
        df["moving_avg"] = df["count"].rolling(window=3, center=True).mean()

        fig = go.Figure()

        # Add main line with gradient
        fig.add_trace(
            go.Scatter(
                x=df["date"],
                y=df["count"],
                mode="lines+markers",
                name="Daily Visitors",
                line=dict(color="#6366f1", width=3, shape="spline"),
                marker=dict(size=8, color="#6366f1", line=dict(width=2, color="white")),
                fill="tonexty",
                fillcolor="rgba(99, 102, 241, 0.1)",
            )
        )

        # Add trend line
        fig.add_trace(
            go.Scatter(
                x=df["date"],
                y=df["moving_avg"],
                mode="lines",
                name="Trend (3-day avg)",
                line=dict(color="#ef4444", width=2, dash="dash"),
                opacity=0.7,
            )
        )

        # Calculate growth rate
        if len(df) > 1:
            current = df["count"].iloc[-1]
            previous = df["count"].iloc[-2]
            growth_rate = ((current - previous) / previous * 100) if previous > 0 else 0
            growth_text = (
                f"📈 +{growth_rate:.1f}%"
                if growth_rate > 0
                else f"📉 {growth_rate:.1f}%"
            )
        else:
            growth_text = "📊 New data"

        fig.update_layout(
            title=dict(
                text=f"Visitor Trends {growth_text}",
                x=0.5,
                font=dict(size=18, color="#1f2937"),
            ),
            xaxis=dict(
                title="Date",
                showgrid=True,
                gridcolor="rgba(0,0,0,0.05)",
                zeroline=False,
                tickformat="%b %d",
            ),
            yaxis=dict(
                title="Unique Visitors",
                showgrid=True,
                gridcolor="rgba(0,0,0,0.05)",
                zeroline=False,
            ),
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
            margin=dict(l=50, r=20, t=60, b=50),
            hovermode="x unified",
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Inter"),
            showlegend=True,
            legend=dict(
                orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1
            ),
        )

        # Focus on last 14 days
        if not df.empty:
            end_date = df["date"].max()
            start_date = end_date - pd.Timedelta(days=14)
            fig.update_layout(xaxis_range=[start_date, end_date])

    else:
        fig = go.Figure()
        fig.update_layout(
            title="No visitor data available",
            xaxis_title="Date",
            yaxis_title="Unique Visitors",
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
        )

    return fig


def create_modern_endpoints_figure(endpoints_data):
    """Create a modern horizontal bar chart for page visits."""
    # Define page mapping locally
    page_mapping = {
        "/": "Home",
        "/download": "Download",
        "/pgv": "PGV",
        "/submit": "Submit",
        "/blast": "BLAST",
        "/wiki": "Wiki",
        "/metrics": "Metrics",
        "/starfish": "Starfish",
        "/about": "About",
    }

    if endpoints_data:
        filtered_data = [
            (page_mapping.get(row[0], row[0]), row[1])
            for row in endpoints_data
            if row[0] in page_mapping
        ]

        if filtered_data:
            endpoints = [row[0] for row in filtered_data]
            counts = [row[1] for row in filtered_data]

            # Create color gradient
            colors = px.colors.qualitative.Set3[: len(endpoints)]

            fig = go.Figure(
                go.Bar(
                    y=endpoints,
                    x=counts,
                    orientation="h",
                    marker=dict(color=colors, line=dict(color="white", width=1)),
                    text=counts,
                    textposition="auto",
                    texttemplate="%{text}",
                    textfont=dict(color="white", size=12),
                )
            )

            fig.update_layout(
                title=dict(
                    text="Page Popularity", x=0.5, font=dict(size=18, color="#1f2937")
                ),
                xaxis=dict(
                    title="Unique Visitors",
                    showgrid=True,
                    gridcolor="rgba(0,0,0,0.05)",
                    zeroline=False,
                ),
                yaxis=dict(showgrid=False, zeroline=False),
                plot_bgcolor="white",
                paper_bgcolor="white",
                height=400,
                margin=dict(l=50, r=20, t=60, b=50),
                bargap=0.3,
                showlegend=False,
            )

        else:
            fig = go.Figure()
            fig.update_layout(
                title="No page visit data available",
                plot_bgcolor="white",
                paper_bgcolor="white",
                height=400,
            )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="No page visit data available",
            plot_bgcolor="white",
            paper_bgcolor="white",
            height=400,
        )

    return fig


def create_modern_map_figure(locations):
    """Create a modern map visualization."""
    fig = go.Figure()

    fig.update_layout(
        mapbox=dict(
            style="carto-positron",
            zoom=1,
            center=dict(lat=20, lon=0),
        ),
        margin={"r": 0, "t": 30, "l": 0, "b": 0},
        height=400,
        title=dict(
            text="Global Visitor Distribution",
            x=0.5,
            font=dict(size=18, color="#1f2937"),
        ),
        paper_bgcolor="white",
    )

    if locations and len(locations) > 0:
        df = pd.DataFrame(locations)

        # Count visitors per location
        location_counts = (
            df.groupby(["lat", "lon", "city", "country"])
            .size()
            .reset_index(name="visits")
        )

        # Create size gradient
        min_visits = location_counts["visits"].min()
        max_visits = location_counts["visits"].max()
        location_counts["marker_size"] = 8 + (
            location_counts["visits"] - min_visits
        ) * (25 / (max_visits - min_visits if max_visits > min_visits else 1))

        fig.add_trace(
            go.Scattermapbox(
                lat=location_counts["lat"],
                lon=location_counts["lon"],
                mode="markers",
                marker=dict(
                    size=location_counts["marker_size"],
                    color=location_counts["visits"],
                    colorscale="Viridis",
                    opacity=0.8,
                    sizemode="diameter",
                    colorbar=dict(title="Visits", x=0.95, len=0.8),
                ),
                text=location_counts.apply(
                    lambda row: f"<b>{row['city']}, {row['country']}</b><br>👥 {int(row['visits'])} visitors",
                    axis=1,
                ),
                hoverinfo="text",
                hovertemplate="%{text}<extra></extra>",
            )
        )

        # Adjust center and zoom
        center_lat = location_counts["lat"].mean()
        center_lon = location_counts["lon"].mean()
        fig.update_layout(
            mapbox=dict(center=dict(lat=center_lat, lon=center_lon), zoom=1.5)
        )

    return fig


def analyze_telemetry() -> Dict[str, Any]:
    """Analyze telemetry data and create visualizations."""
    try:
        # Fetch basic telemetry data
        data = fetch_telemetry_data()
        if data is None:
            return {
                "unique_users": 0,
                "total_requests": 0,
                "time_series": go.Figure(),
                "endpoints": go.Figure(),
                "map": go.Figure(),
                "locations": [],
            }

        # Create modern visualizations
        time_series_fig = create_modern_time_series_figure(data["time_series_data"])

        # Create endpoints visualization
        endpoints_fig = create_modern_endpoints_figure(data["endpoints_data"])

        # Get cached locations and create map
        locations = get_cached_locations()
        map_fig = create_modern_map_figure(locations)

        # Calculate total requests
        total_requests = (
            sum(row[1] for row in data["endpoints_data"] if row[0] in page_mapping)
            if data["endpoints_data"]
            else 0
        )

        return {
            "unique_users": data["unique_users"],
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


def blast_limit_decorator(f):
    """Decorator to limit BLAST operations per user"""

    @wraps(f)
    def wrapped(*args, **kwargs):
        try:
            client_ip = get_client_ip()

            # Skip rate limiting for development IPs
            if is_development_ip(client_ip):
                return f(*args, **kwargs)

            # Only check limits for BLAST-related endpoints
            if request.path == "/api/blast-submit":
                limit_info = get_blast_limit_info(client_ip)

                if limit_info["remaining"] <= 0:
                    logger.warning(f"BLAST limit exceeded for IP: {client_ip}")
                    raise PreventUpdate("Hourly BLAST limit exceeded")

                # Log the BLAST request
                logger.debug(f"BLAST submission from IP: {client_ip}")
                log_request(client_ip, "/api/blast-submit")

            return f(*args, **kwargs)

        except PreventUpdate:
            raise
        except Exception as e:
            logger.error(f"Error in blast_limit_decorator: {str(e)}")
            raise

    return wrapped

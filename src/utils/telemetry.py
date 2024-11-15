from datetime import datetime
import warnings
import logging
import pandas as pd
from datetime import datetime, timedelta
import plotly.express as px
import plotly.graph_objects as go
from flask import request
from ip2geotools.databases.noncommercial import DbIpCity

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)


from src.components.sql_engine import telemetry_session_factory, telemetry_connected

def get_client_ip():
    """Get the client's IP address from the request."""
    if request.headers.get('X-Forwarded-For'):
        # If behind a proxy, get real IP
        return request.headers.get('X-Forwarded-For').split(',')[0]
    return request.remote_addr


# use the telemetry_engine to log request_logs
def log_request(ip_address, endpoint):
    """Log request details to telemetry database."""
    if not telemetry_connected:
        return

    session = telemetry_session_factory()
    
    try:
        query = """
        INSERT INTO request_logs (ip_address, endpoint, timestamp)
        VALUES (:ip, :endpoint, :timestamp)
        """
        
        session.execute(
            query,
            {
                "ip": ip_address,
                "endpoint": endpoint, 
                "timestamp": datetime.now()
            }
        )
        session.commit()
        
    except Exception as e:
        session.rollback()
    finally:
        session.close()

def fetch_telemetry_data():
    """Fetch telemetry data from the database."""
    session = telemetry_session_factory()
    try:
        # Get unique users count
        unique_users_query = "SELECT COUNT(DISTINCT ip_address) as count FROM request_logs"
        unique_users = session.execute(unique_users_query).scalar() or 0

        # Get time series data
        time_series_query = """
        SELECT DATE(timestamp) as date, COUNT(*) as count 
        FROM request_logs 
        GROUP BY DATE(timestamp) 
        ORDER BY date
        """
        time_series_data = session.execute(time_series_query).fetchall()
        
        # Get endpoint statistics
        endpoints_query = """
        SELECT endpoint, COUNT(*) as count 
        FROM request_logs 
        GROUP BY endpoint 
        ORDER BY count DESC
        """
        endpoints_data = session.execute(endpoints_query).fetchall()

        return {
            "unique_users": unique_users,
            "time_series_data": time_series_data,
            "endpoints_data": endpoints_data
        }

    except Exception as e:
        logger.error(f"Error fetching telemetry data: {str(e)}")
        return None
    finally:
        session.close()

def create_time_series_figure(time_series_data):
    """Create time series visualization."""
    if time_series_data:
        dates = [row[0] for row in time_series_data]
        counts = [row[1] for row in time_series_data]
        fig = px.line(
            x=dates, 
            y=counts,
            title="Daily Requests",
            labels={"x": "Date", "y": "Number of Requests"}
        )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="Daily Requests (No Data)",
            xaxis_title="Date",
            yaxis_title="Number of Requests"
        )
    return fig

def create_endpoints_figure(endpoints_data):
    """Create endpoints visualization."""
    if endpoints_data:
        endpoints = [row[0] for row in endpoints_data]
        counts = [row[1] for row in endpoints_data]
        fig = px.bar(
            x=endpoints,
            y=counts,
            title="Requests by Endpoint",
            labels={"x": "Endpoint", "y": "Number of Requests"}
        )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="Requests by Endpoint (No Data)",
            xaxis_title="Endpoint",
            yaxis_title="Number of Requests"
        )
    return fig

def get_ip_locations():
    """Get location data for IP addresses."""
    session = telemetry_session_factory()
    try:
        # Get unique IP addresses
        query = "SELECT DISTINCT ip_address FROM request_logs"
        ip_addresses = session.execute(query).fetchall()
        
        # Get location data for each IP
        locations = []
        for (ip,) in ip_addresses:
            try:
                # Skip localhost/private IPs
                if ip in ('127.0.0.1', 'localhost') or ip.startswith(('192.168.', '10.', '172.')):
                    continue
                    
                response = DbIpCity.get(ip)
                locations.append({
                    'ip': ip,
                    'lat': response.latitude,
                    'lon': response.longitude,
                    'city': response.city,
                    'country': response.country
                })
            except Exception as e:
                logger.warning(f"Could not get location for IP {ip}: {str(e)}")
                continue
                
        return locations
        
    except Exception as e:
        logger.error(f"Error fetching IP locations: {str(e)}")
        return []
    finally:
        session.close()

def create_map_figure(locations):
    """Create map visualization of user locations."""
    if locations:
        df = pd.DataFrame(locations)
        fig = px.scatter_mapbox(
            df,
            lat='lat',
            lon='lon',
            hover_data=['city', 'country'],
            zoom=1,
            title="User Locations"
        )
        fig.update_layout(
            mapbox_style="open-street-map",
            margin={"r":0,"t":30,"l":0,"b":0}
        )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="User Locations (No Data)",
        )
    return fig


def analyze_telemetry():
    """Analyze telemetry data and create visualizations."""
    # Fetch data
    data = fetch_telemetry_data()
    locations = get_ip_locations()
    
    if data is None:
        empty_fig = go.Figure()
        return {
            "unique_users": 0,
            "total_requests": 0,
            "time_series": empty_fig,
            "endpoints": empty_fig,
            "map": empty_fig
        }

    # Create visualizations
    time_series_fig = create_time_series_figure(data["time_series_data"])
    endpoints_fig = create_endpoints_figure(data["endpoints_data"])
    map_fig = create_map_figure(locations)
    
    # Calculate total requests
    total_requests = sum(row[1] for row in data["endpoints_data"]) if data["endpoints_data"] else 0

    return {
        "unique_users": data["unique_users"],
        "total_requests": total_requests,
        "time_series": time_series_fig,
        "endpoints": endpoints_fig,
        "map": map_fig
    }

def count_blast_submissions(ip_address, hours=1):
    """Count BLAST submissions from an IP address in the last N hours."""
    session = telemetry_session_factory()
    try:
        query = """
        SELECT COUNT(*) as count 
        FROM request_logs 
        WHERE ip_address = :ip 
        AND endpoint = '/api/blast-submit'
        AND datetime(timestamp) >= datetime('now', :hours_ago)
        """
        result = session.execute(
            query, 
            {
                "ip": ip_address,
                "hours_ago": f'-{hours} hours'
            }
        ).scalar() or 0
        return result
    except Exception as e:
        logger.error(f"Error counting BLAST submissions: {str(e)}")
        return 0
    finally:
        session.close()

def get_blast_limit_info(ip_address):
    """Get rate limit info for an IP address."""
    try:
        submissions = count_blast_submissions(ip_address)
        limit = 10  # Configurable limit per hour
        remaining = max(0, limit - submissions)
        return {
            "remaining": remaining,
            "limit": limit,
            "submissions": submissions
        }
    except Exception as e:
        logger.error(f"Error getting BLAST limit info: {str(e)}")
        return {
            "remaining": 0,
            "limit": 10,
            "submissions": 10
        }

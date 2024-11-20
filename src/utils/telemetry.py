from datetime import datetime
import warnings
import logging
import pandas as pd
from datetime import datetime, timedelta
import plotly.express as px
import plotly.graph_objects as go
from flask import request
from functools import wraps
from ip2geotools.databases.noncommercial import DbIpCity
from dash.exceptions import PreventUpdate

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

from src.components.sql_engine import telemetry_session_factory, telemetry_connected

# Define valid pages
page_mapping = {
    '/': 'Home',
    '/download': 'Download',
    '/pgv': 'PGV',
    '/submit': 'Submit',
    '/blast': 'BLAST',
    '/wiki': 'Wiki',
    '/metrics': 'Metrics',
    '/starfish': 'Starfish',
    '/about': 'About'
}


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
    
    # Define valid pages as a comma-separated string of quoted values
    valid_endpoints = "'/','/download','/pgv','/submit','/blast','/wiki','/metrics','/starfish','/about'"
    
    try:
        # Get unique users count (still count all unique users)
        unique_users_query = "SELECT COUNT(DISTINCT ip_address) as count FROM request_logs"
        unique_users = session.execute(unique_users_query).scalar() or 0

        # Get time series data for valid endpoints only
        time_series_query = f"""
        SELECT DATE(timestamp) as date, COUNT(*) as count 
        FROM request_logs 
        WHERE endpoint IN ({valid_endpoints})
        GROUP BY DATE(timestamp) 
        ORDER BY date
        """
        time_series_data = session.execute(time_series_query).fetchall()
        
        # Get endpoint statistics for valid endpoints only
        endpoints_query = f"""
        SELECT endpoint, COUNT(*) as count 
        FROM request_logs 
        WHERE endpoint IN ({valid_endpoints})
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
    """Create time series visualization focusing on recent activity."""
    if time_series_data:
        # Convert dates to pandas datetime objects
        df = pd.DataFrame(time_series_data, columns=['date', 'count'])
        df['date'] = pd.to_datetime(df['date'])  # Convert to datetime
        
        fig = px.line(
            df,
            x='date', 
            y='count',
            title="Daily Requests (Last 7 Days)",
            labels={"date": "Date", "count": "Number of Requests"}
        )
        
        # Focus on the last 7 days
        if not df.empty:
            end_date = df['date'].max()
            start_date = end_date - pd.Timedelta(days=7)
            
            fig.update_layout(
                xaxis_range=[start_date, end_date],
                xaxis_title="Date",
                yaxis_title="Number of Requests",
                showlegend=False
            )
            
            # Improve tick format
            fig.update_xaxes(
                dtick="D1",  # Show daily ticks
                tickformat="%Y-%m-%d",  # Date format
            )
            
            # Add markers to show actual data points
            fig.update_traces(mode='lines+markers')
            
    else:
        fig = go.Figure()
        fig.update_layout(
            title="Daily Requests (No Data)",
            xaxis_title="Date",
            yaxis_title="Number of Requests"
        )
    
    # General layout improvements
    fig.update_layout(
        height=400,
        margin=dict(l=50, r=20, t=40, b=50),
        plot_bgcolor='white',
        paper_bgcolor='white',
        hovermode='x unified'
    )
    
    # Grid lines
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGrey')
    
    return fig

def create_endpoints_figure(endpoints_data):
    """Create endpoints visualization for main app pages."""
    if endpoints_data:
        # Filter and transform the data
        filtered_data = [
            (page_mapping.get(row[0], row[0]), row[1]) 
            for row in endpoints_data 
            if row[0] in page_mapping
        ]
        
        if filtered_data:  # Check if we have any data after filtering
            endpoints = [row[0] for row in filtered_data]
            counts = [row[1] for row in filtered_data]
            
            fig = px.bar(
                x=endpoints,
                y=counts,
                title="Page Visits",
                labels={"x": "Page", "y": "Number of Visits"}
            )
            
            # Customize the layout
            fig.update_layout(
                xaxis_tickangle=-45,  # Angle the labels for better readability
                bargap=0.3,           # Adjust space between bars
            )
        else:
            fig = go.Figure()
            fig.update_layout(
                title="Page Visits (No Data)",
                xaxis_title="Page",
                yaxis_title="Number of Visits"
            )
    else:
        fig = go.Figure()
        fig.update_layout(
            title="Page Visits (No Data)",
            xaxis_title="Page",
            yaxis_title="Number of Visits"
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
    fig = go.Figure(go.Scattermapbox())
    
    # Base layout for both empty and populated maps
    fig.update_layout(
        mapbox=dict(
            style="carto-positron",
            zoom=1,
            center=dict(lat=20, lon=0),
        ),
        margin={"r":0,"t":30,"l":0,"b":0},
        height=400,
        title="User Locations"
    )
    
    if locations and len(locations) > 0:
        df = pd.DataFrame(locations)
        
        # Count visitors per location
        location_counts = df.groupby(['lat', 'lon', 'city', 'country']).size().reset_index(name='visits')
        
        # Calculate marker sizes (scale between 10 and 40 based on visit count)
        min_visits = location_counts['visits'].min()
        max_visits = location_counts['visits'].max()
        location_counts['marker_size'] = 10 + (location_counts['visits'] - min_visits) * (
            30 / (max_visits - min_visits if max_visits > min_visits else 1)
        )
        
        fig.add_trace(go.Scattermapbox(
            lat=location_counts['lat'],
            lon=location_counts['lon'],
            mode='markers',
            marker=dict(
                size=location_counts['marker_size'],
                color='indigo',
                opacity=0.7,
                sizemode='diameter'
            ),
            text=location_counts.apply(
                lambda row: f"{row['city']}, {row['country']}<br>{int(row['visits'])} visits", 
                axis=1
            ),
            hoverinfo='text'
        ))
        
        # Adjust center and zoom if we have points
        center_lat = location_counts['lat'].mean()
        center_lon = location_counts['lon'].mean()
        fig.update_layout(
            mapbox=dict(
                center=dict(lat=center_lat, lon=center_lon),
                zoom=1.5
            )
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
    
    # Calculate total requests only for valid pages
    total_requests = sum(
        row[1] for row in data["endpoints_data"] 
        if row[0] in page_mapping
    ) if data["endpoints_data"] else 0

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
    
def blast_limit_decorator(f):
    """Decorator to limit BLAST operations per user"""
    @wraps(f)
    def wrapped(*args, **kwargs):
        try:
            # Get client IP using existing function
            client_ip = get_client_ip()
            
            # Use existing limit checking function
            limit_info = get_blast_limit_info(client_ip)
            
            if limit_info["remaining"] <= 0:
                logger.warning(f"BLAST limit exceeded for IP: {client_ip}")
                raise PreventUpdate("Hourly BLAST limit exceeded")
            
            # Log the BLAST request
            log_request(client_ip, '/api/blast-submit')
            
            # Execute the callback
            return f(*args, **kwargs)
            
        except PreventUpdate:
            raise
        except Exception as e:
            logger.error(f"Error in blast_limit_decorator: {str(e)}")
            raise
    
    return wrapped
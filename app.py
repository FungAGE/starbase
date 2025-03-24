import warnings
import logging
from flask_compress import Compress
from werkzeug.middleware.proxy_fix import ProxyFix

import dash
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask, request, jsonify
from flask_limiter import Limiter
from sqlalchemy import text
import os
from sqlalchemy.pool import QueuePool
from sqlalchemy import create_engine

logging.basicConfig(
    level=logging.ERROR,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
    ]
)

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.ERROR)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

from src.components import navmenu
from src.components.callbacks import create_feedback_button
from src.utils.telemetry import log_request, get_client_ip, is_development_ip, maintain_ip_locations
from src.config.cache import cache, cache_dir, cleanup_old_cache
from src.config.database import TelemetrySession, SubmissionsSession
from src.config.settings import DB_PATHS

_dash_renderer._set_react_version("18.2.0")

external_stylesheets = [
    "https://unpkg.com/@mantine/core@7.11.0/styles.css",
    "https://unpkg.com/@mantine/dates@7.11.0/styles.css",
    "https://unpkg.com/@mantine/code-highlight@7.11.0/styles.css",
    "https://unpkg.com/@mantine/charts@7.11.0/styles.css",
    "https://unpkg.com/@mantine/carousel@7.11.0/styles.css",
    "https://unpkg.com/@mantine/notifications@7.11.0/styles.css",
    "https://unpkg.com/@mantine/nprogress@7.11.0/styles.css",
    dbc.icons.BOOTSTRAP,
    dbc.themes.BOOTSTRAP,
    "https://code.jquery.com/jquery-3.6.0.min.js",
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "/assets/styles.css",
    "https://unpkg.com/tabulator-tables@6.2.5/dist/css/tabulator.min.css",
    "https://cdn.jsdelivr.net/npm/ag-grid-community@30.0.0/styles/ag-grid.css",
    "https://cdn.jsdelivr.net/npm/ag-grid-community@30.0.0/styles/ag-theme-alpine.css",
]

external_scripts = [
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "https://unpkg.com/tabulator-tables@6.2.5/dist/js/tabulator.min.js",
    "https://unpkg.com/micromodal/dist/micromodal.min.js",
]

# Initialize Flask first
server = Flask(__name__)
server.wsgi_app = ProxyFix(server.wsgi_app, x_for=1, x_proto=1)
Compress(server)

server.config.update(
    MAX_CONTENT_LENGTH=10 * 1024 * 1024,  # 10MB limit
    CACHE_TYPE='SimpleCache',
    CACHE_DEFAULT_TIMEOUT=300,
    SEND_FILE_MAX_AGE_DEFAULT=0,
    COMPRESS_MIMETYPES=['text/html', 'text/css', 'application/javascript'],
    COMPRESS_LEVEL=6,
    COMPRESS_ALGORITHM=['gzip', 'br']
)

cache.init_app(server)
cleanup_old_cache()

# Initialize Dash app
app = Dash(
    __name__,
    server=server,
    use_pages=True,
    pages_folder="src/pages",
    suppress_callback_exceptions=True,
    title="starbase",
    external_stylesheets=external_stylesheets,
    external_scripts=external_scripts,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    update_title=None
)

limiter = Limiter(
    get_client_ip,
    app=server,
    storage_uri="memory://",
    default_limits=[]
)

def initialize_app():
    """Initialize all app components."""
    with server.app_context():
        from src.database.migrations import create_database_indexes
        create_database_indexes()
                
        cache.init_app(server)
        cleanup_old_cache()
                
    maintain_ip_locations()

# Create a context for cache operations
with server.app_context():
    cleanup_old_cache()
    initialize_app()

IPSTACK_API_KEY = os.environ.get('IPSTACK_API_KEY') or os.getenv('IPSTACK_API_KEY')

limiter = Limiter(
    get_client_ip,
    app=server,
    storage_uri="memory://",
    default_limits=["200 per minute"],
    strategy="fixed-window-elastic-expiry"  # More forgiving rate limiting
)

def serve_app_layout():
    return dmc.MantineProvider(
        html.Div([
            dmc.NotificationProvider(position="bottom-right"),
            html.Div(id="notifications-container"),
            dcc.Location(id="url", refresh=False),
            navmenu.navmenu(),
            html.Div(dash.page_container),
            create_feedback_button(),
        ])
    )

@limiter.request_filter
def log_rate_limit():
    remote_addr = get_client_ip()
    if is_development_ip(remote_addr):
        return True
    logging.info(f"Rate limit hit by IP: {remote_addr}")
    return False

@app.server.route('/api/blast-submit', methods=['POST'])
@limiter.limit("10 per hour")
def check_blast_limit():
    remote_addr = get_client_ip()
    logging.info(f"BLAST submission from IP: {remote_addr}")
    return {"allowed": True}

@app.server.before_request
def log_request_info():
    """Log every request to the telemetry database"""
    try:
        # Skip static files and API endpoints
        if request.path.startswith('/static/') or request.path.startswith('/_dash-'):
            return
            
        client_ip = get_client_ip()
        endpoint = request.path
        
        # Log the request
        log_request(client_ip, endpoint)
        
    except Exception as e:
        logger.error(f"Error in request logging middleware: {str(e)}")

@app.server.route('/health/telemetry')
def telemetry_health():
    """Check telemetry system health"""
    try:
        session = TelemetrySession()
        
        test_ip = "127.0.0.1"
        test_endpoint = "/"
        log_request(test_ip, test_endpoint)
        
        result = session.execute(text("SELECT COUNT(*) FROM request_logs")).scalar()
        
        return {
            "status": "healthy",
            "record_count": result
        }, 200
    except Exception as e:
        return {
            "status": "unhealthy",
            "error": str(e)
        }, 503
    finally:
        session.close()

@app.server.route('/api/refresh-telemetry', methods=['POST'])
def refresh_telemetry():
    """Endpoint to refresh telemetry data"""
    try:
        maintain_ip_locations(IPSTACK_API_KEY)
        cache.delete('telemetry_data')
        return {"status": "success"}, 200
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}, 500

def check_submissions_db():
    """Verify submissions database is accessible and properly configured"""
    try:
        session = SubmissionsSession()
        test_query = """
        SELECT COUNT(*) FROM submissions
        """
        result = session.execute(text(test_query)).scalar()
        logger.info(f"Submissions database check passed. Current submissions: {result}")
        return True
    except Exception as e:
        logger.error(f"Submissions database check failed: {str(e)}")
        return False

@app.server.route('/api/cache/status')
def cache_status():
    """Check cache status"""
    try:
        # Just return basic cache info
        return jsonify({
            "status": "success",
            "cache_dir": cache_dir,
            "cache_type": cache.config['CACHE_TYPE']
        }), 200
    except Exception as e:
        logger.error(f"Cache status check failed: {str(e)}")
        return jsonify({
            "status": "error",
            "error": str(e)
        }), 500

@app.server.route('/api/cache/refresh', methods=['POST'])
@limiter.limit("1 per minute")
def refresh_cache():
    """Force refresh of all cached data"""
    try:
        cache.clear()
        cleanup_old_cache()
        
        return jsonify({
            "status": "success",
            "message": "Cache cleared and old files cleaned up"
        }), 200
    except Exception as e:
        logger.error(f"Cache refresh failed: {str(e)}")
        return jsonify({
            "status": "error",
            "error": str(e)
        }), 500

@app.server.route('/api/cache/cleanup', methods=['POST'])
@limiter.limit("1 per minute")
def force_cache_cleanup():
    """Force cleanup of old cache files"""
    try:
        cleanup_old_cache()
        return jsonify({
            "status": "success",
            "message": "Cache cleanup completed"
        }), 200
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return jsonify({
            "status": "error",
            "error": str(e)
        }), 500

app.layout = serve_app_layout

# Add global error handlers
@server.errorhandler(500)
def handle_500(e):
    logger.error(f"Internal Server Error: {str(e)}")
    return {
        "error": "Internal Server Error",
        "message": "The server encountered an error. Please try again later."
    }, 500

@server.errorhandler(429)
def handle_429(e):
    return {
        "error": "Too Many Requests",
        "message": "Please wait before making more requests."
    }, 429

# Create database URLs
DATABASE_URLS = {
    'starbase': f"sqlite:///{DB_PATHS['starbase']}",
    'submissions': f"sqlite:///{DB_PATHS['submissions']}",
    'telemetry': f"sqlite:///{DB_PATHS['telemetry']}"
}

# Create engines with connection pooling
engines = {
    name: create_engine(
        url,
        poolclass=QueuePool,
        pool_size=10,
        max_overflow=20,
        pool_timeout=30,
        pool_recycle=1800
    ) for name, url in DATABASE_URLS.items()
}

if __name__ == "__main__":
    app.run_server(debug=False)

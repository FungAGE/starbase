import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
import dash
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask
from flask import request

from flask_limiter import Limiter
import logging

from functools import wraps
from flask import jsonify
import secrets
import os

from src.components import navmenu
from src.utils.telemetry import log_request, get_client_ip, is_development_ip, maintain_ip_locations
from src.components.sql_manager import refresh_cache


logging.basicConfig(level=logging.ERROR)

# from src.components.sql_manager import precompute_all
from src.components.cache import cache

# from src.utils.blastdb import create_dbs
from src.components.sql_engine import sql_connected


import warnings
import logging

warnings.filterwarnings("ignore")
if not logging.getLogger().hasHandlers():
    logging.basicConfig(level=logging.ERROR)
    logging.getLogger("matplotlib.font_manager").disabled = True

_dash_renderer._set_react_version("18.2.0")

external_stylesheets = [
    dmc.styles.ALL,
    dbc.icons.BOOTSTRAP,
    dbc.themes.BOOTSTRAP,
    "/assets/styles.css",
    "https://unpkg.com/tabulator-tables@6.2.5/dist/css/tabulator.min.css",
]

external_scripts = [
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "https://unpkg.com/tabulator-tables@6.2.5/dist/js/tabulator.min.js",
    "https://unpkg.com/micromodal/dist/micromodal.min.js",
]

server = Flask(__name__)

server.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024

# Initialize Dash app with the Flask server
app = Dash(
    __name__,
    server=server,
    use_pages=True,
    pages_folder="src/pages",
    suppress_callback_exceptions=True,
    title="starbase",
    external_stylesheets=external_stylesheets,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

# Set up cache with app
cache.init_app(server)

limiter = Limiter(
    get_client_ip,
    app=server,
    default_limits=["60 per day", "20 per hour"],
)

@limiter.request_filter
def log_rate_limit():
    remote_addr = get_client_ip()
    # Skip rate limit logging for development IPs
    if is_development_ip(remote_addr):
        return True  # Bypass rate limit for local IPs
    logging.info(f"Rate limit hit by IP: {remote_addr}")
    return False

@app.server.before_request
def before_request_func():
    # Log requests only if not a development IP
    if not is_development_ip(get_client_ip()):
        log_request(get_client_ip(), request.path)

@app.server.before_first_request
def initialize_telemetry():
    maintain_ip_locations()

# Create a secure token for the maintenance endpoint
MAINTENANCE_TOKEN = os.getenv('MAINTENANCE_TOKEN', secrets.token_urlsafe(32))

def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        auth_token = request.headers.get('Authorization')
        if not auth_token or auth_token != f'Bearer {MAINTENANCE_TOKEN}':
            return jsonify({"error": "Unauthorized"}), 401
        return f(*args, **kwargs)
    return decorated

@app.server.route('/maintenance/refresh-cache', methods=['POST'])
@requires_auth
def refresh_cache_endpoint():
    result = refresh_cache()
    return jsonify(result)
def serve_app_layout():
    return dmc.MantineProvider(
        html.Div(
            [
                dmc.NotificationProvider(position="top-center"),
                html.Div(id="notifications-container"),
                dcc.Location(id="url", refresh=False),
                navmenu.navmenu(buttons_disabled=not sql_connected),
                html.Div(dash.page_container),
            ]
        )
    )


app.layout = serve_app_layout

@app.server.route('/api/blast-submit', methods=['POST'])
@limiter.limit("10 per hour")  # Adjust limits as needed
def check_blast_limit():
    remote_addr = get_client_ip()
    logging.info(f"BLAST submission from IP: {remote_addr}")
    return {"allowed": True}

if __name__ == "__main__":
    app.run_server(debug=False)

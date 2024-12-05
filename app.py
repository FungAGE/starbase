import warnings
import logging

import dash
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask, request
from flask_limiter import Limiter
import pandas as pd

logging.basicConfig(level=logging.ERROR)
warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

from src.components import navmenu
from src.utils.telemetry import log_request, get_client_ip, is_development_ip, maintain_ip_locations
from src.components.cache import cache
from src.components.sql_engine import sql_connected
from src.components.sql_manager import fetch_meta_data

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

limiter = Limiter(
    get_client_ip,
    app=server,
    default_limits=["60 per day", "20 per hour"],
)

def initialize_cache(app):
    """Initialize cache with commonly used data during app startup."""
    with app.app_context():
        try:
            meta_data = fetch_meta_data()
            if meta_data is None or meta_data.empty:
                logger.error("Failed to fetch metadata for cache initialization")
                return

            cache.set("meta_data", meta_data)
            
            df = pd.DataFrame(meta_data)
            cache.set("taxonomy_options", sorted(df["genus"].dropna().unique()))
            cache.set("family_options", sorted(df["familyName"].dropna().unique()))
            cache.set("navis_options", sorted(df["starship_navis"].dropna().unique()))
            cache.set("haplotype_options", sorted(df["starship_haplotype"].dropna().unique()))
            
            logger.info("Cache initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing cache: {str(e)}", exc_info=True)


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
def before_request_func():
    if not is_development_ip(get_client_ip()):
        log_request(get_client_ip(), request.path)

def initialize_cache(app):
    """Initialize cache with commonly used data during app startup."""
    with app.app_context():
        try:
            meta_data = fetch_meta_data()
            if meta_data is None or meta_data.empty:
                logger.error("Failed to fetch metadata for cache initialization")
                return

            cache.set("meta_data", meta_data)
            
            df = pd.DataFrame(meta_data)
            cache.set("taxonomy_options", sorted(df["genus"].dropna().unique()))
            cache.set("family_options", sorted(df["familyName"].dropna().unique()))
            cache.set("navis_options", sorted(df["starship_navis"].dropna().unique()))
            cache.set("haplotype_options", sorted(df["starship_haplotype"].dropna().unique()))
            
            logger.info("Cache initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing cache: {str(e)}", exc_info=True)

def initialize_app():
    """Initialize all app components."""
    initialize_cache(app)
    maintain_ip_locations()

initialize_app()

def serve_app_layout():
    return dmc.MantineProvider(
        html.Div([
            dmc.NotificationProvider(position="top-center"),
            html.Div(id="notifications-container"),
            dcc.Location(id="url", refresh=False),
            navmenu.navmenu(buttons_disabled=not sql_connected),
            html.Div(dash.page_container),
        ])
    )

app.layout = serve_app_layout

if __name__ == "__main__":
    app.run_server(debug=False)

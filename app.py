from flask import Flask, request
from flask_compress import Compress
from werkzeug.middleware.proxy_fix import ProxyFix
from sqlalchemy import create_engine, text
from sqlalchemy.pool import QueuePool
import dash
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
from dash import Dash, html, dcc, _dash_renderer

from src.components import navmenu
from src.components.callbacks import create_feedback_button
from src.utils.telemetry import log_request, get_client_ip, maintain_ip_locations
from src.config.cache import cache, cleanup_old_cache
from src.config.database import SubmissionsSession
from src.config.settings import DB_PATHS

from src.config.logging import get_logger
logger = get_logger(__name__)

server = Flask(__name__)
server.wsgi_app = ProxyFix(server.wsgi_app, x_for=1, x_proto=1)
Compress(server)

server.config.update(
    MAX_CONTENT_LENGTH=10 * 1024 * 1024,  # 10MB limit
    CACHE_TYPE="SimpleCache",
    CACHE_DEFAULT_TIMEOUT=300,
    SEND_FILE_MAX_AGE_DEFAULT=0,
    COMPRESS_MIMETYPES=["text/html", "text/css", "application/javascript"],
    COMPRESS_LEVEL=6,
    COMPRESS_ALGORITHM=["gzip", "br"],
)

cache.init_app(server)
cleanup_old_cache()

_dash_renderer._set_react_version("18.2.0")


external_stylesheets = [
    "https://cdn.jsdelivr.net/npm/@mantine/core@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/dates@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/code-highlight@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/charts@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/carousel@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/notifications@7.11.0/styles.css",
    "https://cdn.jsdelivr.net/npm/@mantine/nprogress@7.11.0/styles.css",
    dbc.icons.BOOTSTRAP,
    dbc.themes.BOOTSTRAP,
    "https://code.jquery.com/jquery-3.6.0.min.js",
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "/assets/styles.css",
    "https://cdn.jsdelivr.net/npm/tabulator-tables@6.2.5/dist/css/tabulator.min.css",
    "https://cdn.jsdelivr.net/npm/ag-grid-community@30.0.0/styles/ag-grid.css",
    "https://cdn.jsdelivr.net/npm/ag-grid-community@30.0.0/styles/ag-theme-alpine.css",
]

external_scripts = [
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "https://cdn.jsdelivr.net/npm/tabulator-tables@6.2.5/dist/js/tabulator.min.js",
    "https://cdn.jsdelivr.net/npm/micromodal/dist/micromodal.min.js",
]

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
    update_title=None,
)

DATABASE_URLS = {
    "starbase": f"sqlite:///{DB_PATHS['starbase']}",
    "submissions": f"sqlite:///{DB_PATHS['submissions']}",
    "telemetry": f"sqlite:///{DB_PATHS['telemetry']}",
}

engines = {
    name: create_engine(
        url,
        poolclass=QueuePool,
        pool_size=10,
        max_overflow=20,
        pool_timeout=30,
        pool_recycle=1800,
    )
    for name, url in DATABASE_URLS.items()
}


def initialize_app():
    """Initialize app components and perform setup tasks."""
    with server.app_context():
        from src.database.migrations import create_database_indexes

        create_database_indexes()
        cleanup_old_cache()
        maintain_ip_locations()


def serve_app_layout():
    """Define the layout of the Dash app."""
    return dmc.MantineProvider(
        html.Div(
            [
                dmc.NotificationProvider(position="bottom-right"),
                html.Div(id="notifications-container"),
                dcc.Location(id="url", refresh=False),
                navmenu.navmenu(),
                html.Div(dash.page_container),
                create_feedback_button(),
            ]
        )
    )


@app.server.before_request
def log_request_info():
    """Log incoming requests to the telemetry database."""
    try:
        if request.path.startswith(("/static/", "/_dash-")):
            return
        client_ip = get_client_ip()
        log_request(client_ip, request.path)
    except Exception as e:
        logger.error(f"Error in request logging middleware: {str(e)}")


def check_submissions_db():
    """Verify the submissions database is accessible and properly configured."""
    try:
        session = SubmissionsSession()
        result = session.execute(text("SELECT COUNT(*) FROM submissions")).scalar()
        logger.info(f"Submissions database check passed. Current submissions: {result}")
        return True
    except Exception as e:
        logger.error(f"Submissions database check failed: {str(e)}")
        return False


def component_to_dict(component):
    """Convert a Dash component to a dictionary format."""
    if isinstance(component, (str, int, float)):
        return str(component)
    elif isinstance(component, (list, tuple)):
        return [component_to_dict(c) for c in component]
    elif hasattr(component, "children"):
        result = {
            "type": component.__class__.__name__,
            "children": component_to_dict(component.children)
            if component.children is not None
            else None,
        }
        if hasattr(component, "style") and component.style:
            result["style"] = component.style
        if hasattr(component, "className") and component.className:
            result["className"] = component.className
        return result
    return None


app.layout = serve_app_layout

with server.app_context():
    initialize_app()

if __name__ == "__main__":
    app.run_server(debug=False)
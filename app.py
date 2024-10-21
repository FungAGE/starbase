import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
import dash
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask
from flask import request
import pandas as pd

from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
import logging

from src.components import navmenu
from src.utils.telemetry import log_request

logging.basicConfig(level=logging.INFO)

_dash_renderer._set_react_version("18.2.0")

external_stylesheets = [
    dmc.styles.ALL,
    dbc.icons.BOOTSTRAP,
    dbc.themes.BOOTSTRAP,
    "/assets/lib/styles.css",
    "https://unpkg.com/tabulator-tables@6.2.5/dist/css/tabulator.min.css",
    "/assets/lib/micromodal.css",
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
    # external_scripts=external_scripts,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
    ],
)

limiter = Limiter(
    get_remote_address,
    app=server,
    default_limits=["60 per day", "20 per hour"],
)


# Log rate-limited events
@limiter.request_filter
def log_rate_limit():
    remote_addr = get_remote_address()
    logging.info(f"Rate limit hit by IP: {remote_addr}")
    # Alternatively, store the data in a database or monitoring service
    return False  # Don't exclude the request from rate limiting, just log it


@app.server.before_request
def before_request_func():
    log_request(get_remote_address(), request.path)


def serve_app_layout():
    return dmc.MantineProvider(
        html.Div(
            [
                dcc.Location(id="url", refresh=False),
                navmenu.navmenu(),
                html.Div(dash.page_container),
            ]
        )
    )


app.layout = serve_app_layout

if __name__ == "__main__":
    app.run_server(debug=True)

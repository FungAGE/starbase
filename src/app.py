import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
import dash
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask
import pandas as pd
from sqlalchemy import create_engine

from flask_caching import Cache

from src.components import navmenu
from src.components.callbacks import (
    dl_package,
    caching,
)

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
app = Dash(
    __name__,
    server=server,
    use_pages=True,
    suppress_callback_exceptions=True,
    title="starbase",
    external_stylesheets=external_stylesheets,
    # external_scripts=external_scripts,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1"},
    ],
)

CACHE_CONFIG = {
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': '/tmp/dash_cache',
    'CACHE_DEFAULT_TIMEOUT': 300  # optional, but good to specify
}

# Initialize the cache with the Flask server instance
cache = Cache(app.server, config=CACHE_CONFIG)

def serve_app_layout():
    return dmc.MantineProvider(
        html.Div(
            [
                navmenu.navmenu(),
                html.Div(dash.page_container),
                # dcc.Location(id="url", refresh=False),
                dcc.Store("joined-ships"),
                dcc.Store("curated-status"),
                dcc.Store("curated-dataset"),
                dcc.Store("paper-cache"),
                dcc.Store("pie1-cache"),
                dcc.Store("pie2-cache"),
                dcc.Store("phylogeny-cache"),
                dcc.Store("explore-table-cache"),
                dcc.Store("dl-package")                
            ]
        )
    )


app.layout = serve_app_layout

if __name__ == "__main__":
    app.run_server(debug=True)
    # Create a SQLite engine and initialize the table
    engine = create_engine('sqlite:///database_folder/starbase.sqlite')
    query = "SELECT name FROM sqlite_master WHERE type='table'"
    sql_tbls = pd.read_sql_query(query, engine)
    caching(app,engine,cache)
    dl_package(app)


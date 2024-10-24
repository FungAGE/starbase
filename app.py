import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
import dash
from dash import Dash, html, dcc, _dash_renderer
from flask import Flask
import pandas as pd

from src.components import navmenu

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
    app.run_server(debug=False)

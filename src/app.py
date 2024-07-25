import dash_bootstrap_components as dbc
from dash import Dash, html
import dash
from flask import Flask
from src.components import navmenu

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

external_stylesheets = [
    dbc.themes.FLATLY,
    dbc.icons.BOOTSTRAP,
    dbc.themes.BOOTSTRAP,
    "/assets/styles.css",
    # "https://codepen.io/chriddyp/pen/bWLwgP.css",
]

external_scripts = [
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js",
    "/lib/html2canvas.js",
    "/lib/blaster.min.js",
]

server = Flask(__name__)
app = Dash(
    __name__,
    server=server,
    use_pages=True,
    suppress_callback_exceptions=True,
    title="starbase",
    external_stylesheets=external_stylesheets,
    external_scripts=external_scripts,
)


def serve_app_layout():
    return html.Div(
        [navmenu.sidebar(), html.P(dash.page_container, style=CONTENT_STYLE)]
    )


app.layout = serve_app_layout

if __name__ == "__main__":
    app.run_server()

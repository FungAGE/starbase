import dash_bootstrap_components as dbc
from dash import Dash, html
import dash
from flask import Flask
from src.components import navmenu
from src.components.callbacks import (
    dl_package,
    update_fasta_upload,
    update_gff_upload,
    update_dataset,
)

external_stylesheets = [
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


def serve_app_layout():
    return html.Div([navmenu.navmenu(), html.Div(dash.page_container)])


app.layout = serve_app_layout

dl_package(app)
update_fasta_upload(app)
update_gff_upload(app)
update_dataset(app)

if __name__ == "__main__":
    app.run_server()

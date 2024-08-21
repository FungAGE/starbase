import dash_bootstrap_components as dbc
import dash
from dash import Dash, html, dcc
from flask import Flask
from src.components import navmenu
from src.components.callbacks import (
    dl_package,
    update_fasta_upload,
    update_gff_upload,
    load_ship_metadata,
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
    return html.Div(
        [
            navmenu.navmenu(),
            html.Div(dash.page_container),
            dcc.Location(id="url", refresh=False),
            dcc.Store(id="joined-ships"),
        ]
    )


dl_package(app)
update_fasta_upload(app)
update_gff_upload(app)
load_ship_metadata(app)
update_dataset(app)

app.layout = serve_app_layout

if __name__ == "__main__":
    app.run_server()

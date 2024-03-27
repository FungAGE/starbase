import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from dash.dependencies import Input, Output, State

from src.pages import HOME_URL

dash.register_page(
    __name__,
    title='404',
    name='404'
)

# 404 Title
def not_found():
    return html.Div([
            html.H4("Page Not Found"),
            html.P("Sorry, the page you requested could not be found."),
            dcc.Link("Back to Home", href=HOME_URL),
        ])

layout = not_found()
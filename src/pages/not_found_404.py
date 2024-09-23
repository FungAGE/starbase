import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from src.pages import HOME_URL

dash.register_page(__name__, title="404", name="404")


# 404 Title
def not_found():
    return dmc.Container(
        flex=True,
        children=[
            dmc.Grid(
                justify="content",
                children=dmc.GridCol(
                    [
                        dmc.Title("Page Not Found"),
                        dbc.Alert(
                            "Sorry, the page you requested could not be found.",
                            color="warning",
                        ),
                        dcc.Link("Back to Home", href=HOME_URL),
                    ],
                    style={"paddingTop": "20px"},
                ),
            )
        ],
    )


layout = not_found()

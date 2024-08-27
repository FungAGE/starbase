import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html

dash.register_page(__name__)

card_dict = {
    "aaron": {
        "full_name": "Aaron Vogan",
        "img": "assets/images/aaron.png",
        "role": "FungAGE group leader",
        "email": "mailto:aaron.vogan@ebc.uu.se",
    },
    "adrian": {
        "full_name": "Adrian Forsythe",
        "img": "assets/images/adrian.png",
        "role": html.P(
            [
                html.Span(
                    "starbase",
                    className="logo-text",
                ),
                " lead developer",
            ],
            className="card-text",
        ),
        "email": "mailto:adrian.e.forsythe@gmail.com",
    },
    "emile": {
        "full_name": "Emile Gluck-Thaler",
        "img": "assets/images/emile.png",
        "role": html.P(
            [
                html.Span(
                    "starfish",
                    className="logo-text",
                ),
                " lead developer, Gluck-Thaler lab group leader",
            ],
            className="card-text",
        ),
        "email": "mailto:emilegluckthaler@gmail.com",
    },
}


def make_card(name):
    card = dmc.Card(
        children=[
            dmc.CardSection(
                [
                    dmc.Image(
                        src=card_dict[name]["img"],
                        className="img-fluid rounded-start",
                        style={
                            "object-fit": "cover",
                            "width": "100%",
                        },
                    ),
                ],
            ),
            dmc.CardSection(
                [
                    dmc.Group(
                        [
                            dmc.Text(card_dict[name]["full_name"], size="lg"),
                            dmc.Text(card_dict[name]["role"], c="dimmed", size="md"),
                            dbc.Button(
                                html.I(className="bi bi-envelope"),
                                href=card_dict[name]["email"],
                                color="teal",
                                size="lg",
                                className="me-2",
                            ),
                        ],
                        justify="center",
                        mt="md",
                        mb="xs",
                    ),
                ],
            ),
        ],
        withBorder=True,
        shadow="sm",
        radius="md",
        w=250,
    )
    return card


layout = (
    dmc.Container(
        fluid=True,
        children=[
            dmc.Grid(
                justify="center",
                align="center",
                style={"paddingTop": "20px"},
                children=[
                    dmc.Text(
                        [
                            html.Span(
                                "starbase",
                                className="logo-text",
                            ),
                            " was developed by the ",
                            html.A(
                                "FungAGE lab",
                                href="https://fungage.github.io/",
                            ),
                            " in collaboration with the Gluck-Thaler lab.",
                        ],
                        size="lg",
                    ),
                ],
            ),
            dmc.Grid(
                justify="center",
                align="center",
                style={"padding": "10px"},
                children=[
                    dmc.GridCol(
                        span="content",
                        children=[
                            make_card("aaron"),
                            make_card("adrian"),
                            make_card("emile"),
                            dmc.Text(
                                [
                                    html.P(
                                        [
                                            "The source code for ",
                                            html.Span(
                                                "starbase",
                                                className="logo-text",
                                            ),
                                            " webserver will soon be available on GitHub",
                                        ],
                                        className="text-center auto-resize-600",
                                        style={"font-size": "1rem"},
                                    )
                                ],
                                size="lg",
                            ),
                            html.Div(
                                html.Img(
                                    src="assets/images/starbase-map.png",
                                    className="img-fluid auto-resize-600",
                                )
                            ),
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        html.H2(
                                            "Data Availability",
                                            className="text-center",
                                        )
                                    ),
                                    dbc.CardBody(
                                        [
                                            html.P(
                                                [
                                                    "We have been maintaining ",
                                                    html.Span(
                                                        "starbase",
                                                        className="logo-text",
                                                    ),
                                                    " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the meantime, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                ],
                                                style={"font-size": "0.875rem"},
                                            ),
                                            html.Div(
                                                [
                                                    dbc.Button(
                                                        "Download the latest version of starbase.",
                                                        id="dl-button",
                                                        color="primary",
                                                        className="mt-2",
                                                    ),
                                                    dcc.Download(id="dl-package"),
                                                ],
                                                className="text-center",
                                                style={
                                                    "font-size": "0.875rem",
                                                },
                                            ),
                                        ]
                                    ),
                                ],
                                className="auto-resize-600",
                            ),
                        ],
                    ),
                ],
            ),
        ],
    ),
)

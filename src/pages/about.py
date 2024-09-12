import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import html, dcc
from src.components.callbacks import download_ships_card, modal

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
                " lead developer",
                html.Br(),
                "Gluck-Thaler lab group leader",
            ],
        ),
        "email": "mailto:emilegluckthaler@gmail.com",
    },
}


def make_card(name):
    card = dmc.GridCol(
        # span={"xl": 4, "lg": 4, "md": 12, "sm": 12, "xs": 12},
        span="content",
        children=[
            dmc.Card(
                children=[
                    dmc.CardSection(
                        [
                            dmc.Image(
                                src=card_dict[name]["img"],
                                className="img-fluid rounded-start",
                                style={
                                    "object-fit": "cover",
                                    "width": "100%",
                                    "height": "248px",
                                },
                            ),
                        ],
                    ),
                    dmc.CardSection(
                        [
                            dmc.Stack(
                                align="stretch",
                                justify="space-around",
                                children=[
                                    dmc.Group(
                                        justify="flex-start",
                                        gap="xs",
                                        children=[
                                            dmc.Text(
                                                card_dict[name]["full_name"],
                                                size="lg",
                                                style={"padding": "10px"},
                                            ),
                                            dbc.Button(
                                                html.I(className="bi bi-envelope"),
                                                href=card_dict[name]["email"],
                                                color="teal",
                                                size="lg",
                                                className="me-2",
                                            ),
                                        ],
                                    ),
                                    dmc.Text(
                                        card_dict[name]["role"],
                                        c="dimmed",
                                        size="md",
                                        style={"padding": "10px"},
                                    ),
                                ],
                            ),
                        ],
                        style={"height": "150px"},
                    ),
                ],
                withBorder=True,
                shadow="sm",
                radius="md",
                w=250,
            )
        ],
    )
    return card


layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dmc.Grid(
            justify="center",
            align="center",
            style={"paddingTop": "20px"},
            children=[
                dmc.GridCol(
                    span=10,
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
                            className="text-center",
                        ),
                    ],
                ),
            ],
        ),
        dmc.Grid(
            justify="center",
            align="stretch",
            children=[
                make_card("aaron"),
                make_card("adrian"),
                make_card("emile"),
                dmc.GridCol(
                    span=10,
                    children=[
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
                                    className="text-center",
                                )
                            ],
                            size="lg",
                        ),
                        dmc.Center(
                            [
                                html.Div(
                                    html.Img(
                                        src="assets/images/starbase-map.png",
                                        className="auto-resize-750",
                                    )
                                ),
                            ]
                        ),
                        dmc.Center(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [
                                                dmc.Text(
                                                    "Data Availability",
                                                    size="lg",
                                                )
                                            ],
                                            className="card-header-custom",
                                        ),
                                        download_ships_card,
                                    ],
                                    className="auto-resize-750",
                                ),
                                modal,
                            ]
                        ),
                    ],
                ),
            ],
        ),
    ],
)

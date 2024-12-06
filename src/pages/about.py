import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import html, dcc

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
    return dmc.Card(
        children=[
            dmc.CardSection(
                dmc.Image(
                    src=card_dict[name]["img"],
                    fit="cover",
                    style={
                    "minHeight": "250px",
                    "height": "clamp(250px, 8vw, 350px)",
                },
                ),
            ),
            dmc.CardSection(
                dmc.Stack([
                    dmc.Group([
                        dmc.Text(
                            card_dict[name]["full_name"],
                            size="lg",
                            fw=500,
                        ),
                        dmc.Anchor(
                            dmc.ActionIcon(
                                html.I(className="bi bi-envelope"),
                                variant="light", 
                                color="indigo",
                                size="lg",
                                radius="xl",
                            ),
                            href=card_dict[name]["email"],
                        ),
                    ], pos="apart"),
                    dmc.Text(
                        card_dict[name]["role"],
                        c="dimmed",
                        size="md",
                    ),
                ], gap="xs", p="md"),
            ),
        ],
        withBorder=True,
        shadow="sm",
        radius="md",
        style={"width": 280},
    )


layout = dmc.Container(
    size="lg",
    children=[
        # Header Section
        dmc.Paper(
            children=[
                dmc.Title(["About ", html.Span("starbase", className="logo-text")], order=1, mb="md"),
                dmc.Text(
                    [
                        html.Span("starbase", className="logo-text"),
                        " was developed by the ",
                        dmc.Anchor(
                            "FungAGE lab",
                            href="https://fungage.github.io/",
                            underline=False,
                            c="indigo",
                        ),
                        " in collaboration with the Gluck-Thaler lab.",
                    ],
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        
        # Team Section
        dmc.Paper(
            children=[
                dmc.SimpleGrid(
                    cols={"md":3,"sm":1},  # Default number of columns
                    spacing="xl",
                    children=[
                        make_card("aaron"),
                        make_card("adrian"),
                        make_card("emile"),
                    ],
                ),
            ],
            p="xl",
            radius="md",
            withBorder=False,
            mb="xl",
        ),
        
        # Source Code Section
        dmc.Paper(
            children=[
                dmc.Stack([
                    dmc.Text(
                        [
                            "The source code for ",
                            html.Span("starbase", className="logo-text"),
                            " webserver will soon be available on GitHub",
                        ],
                        size="lg",
                    ),
                    dmc.Image(
                        src="assets/images/starbase-map.png",
                        fit="contain",
                        className="auto-resize-750",
                    ),
                ], gap="xl"),
            ],
            p="xl",
            radius="md",
            withBorder=True,
        ),
    ],
    py="xl",
)
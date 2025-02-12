import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash_iconify import DashIconify
import logging

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from src.components.tables import make_paper_table
from src.components.callbacks import download_ships_card
from src.database.sql_engine import check_database_connection
from src.database.sql_manager import get_database_stats

logger = logging.getLogger(__name__)

dash.register_page(__name__, title="starbase Home", name="Home", path="/")

is_connected = check_database_connection()

title = dmc.Paper([
    dmc.Group(
        [
            dmc.Image(
                src="/assets/logos/emblem.svg",
                fit="contain",
                radius="md",
                style={
                    "minWidth": "120px",
                    "minHeight": "120px",
                    "width": "clamp(120px, 8vw, 160px)",
                    "height": "clamp(120px, 8vw, 160px)",
                    "flex": "0 0 auto",  # Prevent image from growing/shrinking
                },
            ),
            dmc.Title(
                [
                    html.Span(
                        "starbase: ",
                        className="logo-text",
                        style={"color": "white"}
                    ),
                    "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                ],
                order=1,
                style={
                    "fontSize": "clamp(1.2rem, 4vw, 2.5rem)",
                    "color": "white",
                    "textAlign": "left",
                    "wordBreak": "break-word",
                    "flex": "1",  # Allow title to take remaining space
                }
            ),
        ],
        className="responsive-group",
        style={
            "width": "100%",
            "maxWidth": "1400px",
            "margin": "0 auto",
            "padding": "clamp(1rem, 2vw, 2rem)",
        }
    )], 
    shadow="sm",
    p={
        "base": "md",
        "sm": "lg",
        "md": "xl"
    },
    radius="md",
    style={
        "backgroundColor": "#2C2E33",
        "marginBottom": "2rem",
    },
    className="responsive-group"
)

working = {
    "wiki": "Catalogue/Wiki of Starship Metadata",
    "submit": "Submission of new Starship sequences",
    "blast": "BLAST/HMMER searches",
}

def create_feature_button(label, href, icon):
    return dmc.Anchor(
        dmc.Button(
            label,
            leftSection=DashIconify(icon=icon),
            variant="gradient",
            gradient={"from": "indigo", "to": "cyan"},
            size="lg",
            radius="md",
            fullWidth=True,
            disabled=not is_connected,
        ),
        href=href,
        style={"textDecoration": "none"},  # Remove underline from link
    )

working_buttons = [
    create_feature_button("BLAST Search", "/blast", "mdi:dna"),
    create_feature_button("Browse Wiki", "/wiki", "mdi:book-open-variant"),
    create_feature_button("Submit to Starbase", "/submit", "mdi:database-plus"),
]

not_working = [
    html.Div(
        ["Synteny/Genome Browser"],
    ),
    html.Div(
        [
            html.Span(
                "starfish",
                className="logo-text",
            ),
            " webserver",
        ],
    ),
]

starship_card = dmc.Paper([
    dmc.Title("What is a Starship?", order=2, mb="md"),
    dmc.Grid([
        dmc.GridCol([
            dmc.Text([
                html.Div(
                    [
                        "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes). They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in ",
                        html.Span(
                                                "Paecilomyces",
                                                style={"font-style": "italic"},
                                            ),
                                            " cheese making in ",
                                            html.Span(
                                                "Penicillium",
                                                style={
                                                    "font-style": "italic",
                                                },
                                            ),
                                            ", and enable the transfer of formaldehyde resistance in ",
                                            html.Span(
                                                "Aspergillus nidulans",
                                                style={
                                                    "font-style": "italic",
                                                },
                                            ),
                                            " and ",
                                            html.Span(
                                                "Penicillium chrysogenum.",
                                                style={
                                                    "font-style": "italic",
                                                },
                                            ),
                                            " Read more about Starships ",
                                            dmc.Anchor(
                                                "here",
                                                href="https://en.wikipedia.org/wiki/Starship_(genetics)",
                                                style={"textDecoration": "none"},
                                            ),
                                            ".",
                                        ],
                                    ),
            ], size="lg", c="dimmed"),
        ], 
        span={"sm": 12, "md": 7}
        ),
        dmc.GridCol([
            dmc.Image(
                src="assets/images/starship-model.png",
                fit="contain",
                radius="md",
                style={"maxWidth": "100%"}
            )
        ], 
        span={"sm": 12, "md": 5}
        ),
    ]),
], shadow="sm", p="xl", radius="md", withBorder=True)


working_features_card = dmc.Paper(
    [
        dmc.Title(
            [
                "What can I currently use ",
                html.Span(
                    "starbase",
                    className="logo-text",
                ),
                " for?",
            ],
            order=2,
            mb="md",
        ),
        dmc.Stack(
            working_buttons,
            gap="md",
        ),
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
    mb="xl",
)

developing_features_card = dmc.Paper(
    [
        dmc.Title(
            [
                "Functions of ",
                html.Span(
                    "starbase",
                    className="logo-text",
                ),
                " under active development:",
            ],
            order=2,
            mb="md",
        ),
        dmc.List(
            [
                dmc.ListItem(
                    dmc.Text(item, size="lg", c="dimmed")
                ) for item in not_working
            ],
            size="lg",
            spacing="sm",
        ),
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
)

def create_hero_section():
    return dmc.Container([
        dmc.Space(h=20),
        dmc.Center(title),
        dmc.Space(h=40),
        dmc.Center(starship_card),
    ], size="xl",flex=True)

def create_features_section():
    return dmc.Container([
        dmc.Grid([
            dmc.GridCol([
                working_features_card,
                developing_features_card
            ], 
            span={"sm": 12, "md": 6}
            ),
            dmc.GridCol(
                download_ships_card, 
                span={"sm": 12, "md": 6}
            )
        ], gutter="xl"),
    ], size="xl", py="xl", flex=True)

def create_publications_section():
    return dmc.Container([
        dmc.Paper([
            dmc.Title("Manuscripts Characterizing Starships", order=2, mb="xl"),
            html.Div(
                make_paper_table(),
                style={"width": "100%"}
            ),
        ], shadow="sm", p="xl", radius="md", withBorder=True, mb="xl"),
    ], size="xl", py="xl", flex=True)

def create_stats_section():
    stats = get_database_stats() if is_connected else {
        "curated_starships": "—",
        "uncurated_starships": "—",
        "species_count": "—",
        "family_count": "—"
    }
    
    return dmc.Container([
        dmc.Title(
            "Database Statistics", 
            order=2, 
            mb="xl"
        ),
        dmc.SimpleGrid([
            dmc.Paper([
                dmc.Text("Total Starships", size="lg", c="dimmed"),
                dmc.Group([
                    dmc.Stack([
                        dmc.Text("Curated", size="sm", c="dimmed"),
                        dmc.Title(
                            f"{stats['curated_starships']:,}", 
                            order=3,
                            c="green"
                        ),
                    ]),
                    dmc.Stack([
                        dmc.Text("Uncurated", size="sm", c="dimmed"),
                        dmc.Title(
                            f"{stats['uncurated_starships']:,}", 
                            order=3,
                            c="orange"
                        ),
                    ]),
                ], gap="apart", grow=True),
            ], p="xl", radius="md", shadow="sm", withBorder=True),
            dmc.Paper([
                dmc.Text("Species", size="lg", c="dimmed"),
                dmc.Title(
                    f"{stats['species_count']:,}", 
                    order=3
                ),
            ], p="xl", radius="md", shadow="sm", withBorder=True),
            dmc.Paper([
                dmc.Text("Starship Families", size="lg", c="dimmed"),
                dmc.Title(
                    f"{stats['family_count']:,}", 
                    order=3
                ),
            ], p="xl", radius="md", shadow="sm", withBorder=True),
        ], 
        cols={
            "base": 1,
            "sm": 2,
            "md": 3,
        },
        spacing="lg"
        ),
    ], size="xl", py="xl")

layout = dmc.MantineProvider([
    dcc.Location(id="url", refresh=False),
    create_hero_section(),
    dmc.Space(h=40),
    create_features_section(),
    dmc.Space(h=40),
    create_stats_section(),
    dmc.Space(h=40),
    create_publications_section(),
    # Database warning if needed
    dmc.Notification(
        title="Database Connection Failed",
        message="Many features will be disabled until connection is re-established.",
        c="red",
        style={"position": "fixed", "top": 20, "right": 20}
    ) if not is_connected else None,
], 
theme={
    "colorScheme": "light",
    "primaryColor": "indigo",
    "components": {
        "Container": {"defaultProps": {"size": "xl"}},
        "Title": {"defaultProps": {"color": "indigo"}},
    }
})

@callback(
    [Output("toggle-paper-view", "leftSection"),
     Output("desktop-table", "style"),
     Output("mobile-cards", "style")],
    Input("toggle-paper-view", "n_clicks"),
    State("toggle-paper-view", "leftSection"),
    prevent_initial_call=True
)
def toggle_paper_view(n_clicks, current_icon):
    if n_clicks is None:
        raise PreventUpdate
        
    is_table = current_icon["icon"] == "tabler:layout-list"
    
    return (
        DashIconify(icon="tabler:layout-cards" if is_table else "tabler:layout-list"),
        {"display": "none" if is_table else "block"},
        {"display": "block" if is_table else "none"}
    )

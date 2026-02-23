import dash
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash_iconify import DashIconify

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from src.components.tables import make_paper_table
from src.database.sql_manager import (
    get_database_stats,
    get_database_version,
    get_alembic_schema_version,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

dash.register_page(__name__, title="starbase Home", name="Home", path="/")

title = dmc.Paper(
    [
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
                            style={"color": "white"},
                        ),
                        "A Database and Toolkit for Exploration of ",
                        html.Span("Starship", style={"fontStyle": "italic"}),
                        " Elements in Fungi",
                    ],
                    order=1,
                    style={
                        "fontSize": "clamp(1.2rem, 4vw, 2.5rem)",
                        "color": "white",
                        "textAlign": "left",
                        "wordBreak": "break-word",
                        "flex": "1",
                    },
                ),
            ],
            className="responsive-group",
            style={
                "width": "100%",
                "maxWidth": "1400px",
                "margin": "0 auto",
                "padding": "clamp(1rem, 2vw, 2rem)",
            },
        )
    ],
    shadow="sm",
    p={"base": "md", "sm": "lg", "md": "xl"},
    radius="md",
    style={
        "backgroundColor": "#2C2E33",
        "marginBottom": "2rem",
    },
    className="responsive-group",
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
            fullWidth=False,
            style={"minWidth": "200px"},
        ),
        href=href,
        style={"textDecoration": "none"},  # Remove underline from link
    )


working_buttons = [
    create_feature_button("BLAST Search", "/blast", "mdi:dna"),
    create_feature_button("Browse Wiki and Download", "/wiki", "mdi:book-open-variant"),
    create_feature_button(
        html.Span(
            [
                "Submit to ",
                html.Span("Starbase", className="logo-text"),
            ]
        ),
        "/submit",
        "mdi:database-plus",
    ),
]

not_working = [
    html.Span(
        ["Synteny/Genome Browser"],
    ),
    html.Span(
        [
            html.Span(
                "starfish",
                className="logo-text",
            ),
            " webserver",
        ],
    ),
]

starship_card = dmc.Paper(
    [
        dmc.Title(
            html.Span(
                [
                    "What is a ",
                    html.Span("Starship", style={"fontStyle": "italic"}),
                    "?",
                ]
            ),
            order=2,
            mb="md",
        ),
        dmc.Grid(
            [
                dmc.GridCol(
                    [
                        dmc.Text(
                            [
                                html.Span("Starships", style={"fontStyle": "italic"}),
                                " are novel family of class II DNA transposons, endemic to Pezizomycotina. ",
                                html.Span("Starships", style={"fontStyle": "italic"}),
                                " can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes). They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in ",
                                html.Span(
                                    "Paecilomyces",
                                    style={"fontStyle": "italic"},
                                ),
                                ", and enable the transfer of formaldehyde resistance in ",
                                html.Span(
                                    "Aspergillus nidulans",
                                    style={
                                        "fontStyle": "italic",
                                    },
                                ),
                                " and ",
                                html.Span(
                                    "Penicillium chrysogenum.",
                                    style={
                                        "fontStyle": "italic",
                                    },
                                ),
                                " Read more about ",
                                html.Span("Starships", style={"fontStyle": "italic"}),
                                " ",
                                dmc.Anchor(
                                    "here",
                                    href="https://en.wikipedia.org/wiki/Starship_(genetics)",
                                    style={"textDecoration": "none"},
                                ),
                                ".",
                            ],
                            size="lg",
                            c="dimmed",
                        ),
                    ],
                    span={"sm": 12, "md": 7},
                ),
                dmc.GridCol(
                    [
                        dmc.Image(
                            src="assets/images/starship-model.png",
                            fit="contain",
                            radius="md",
                            style={"maxWidth": "100%"},
                        )
                    ],
                    span={"sm": 12, "md": 5},
                ),
            ]
        ),
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
)


working_features_card = dmc.Paper(
    [
        dmc.Stack(
            [
                dmc.Group(
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
                            align="center",
                            justify="center",
                            style={"width": "100%", "alignItems": "center"},
                        ),
                    ],
                    justify="space-between",
                    align="center",
                ),
                dmc.Group(
                    [
                        dmc.Title(
                            [
                                "Most functions of ",
                                html.Span(
                                    "starbase",
                                    className="logo-text",
                                ),
                                " are still under active development. We plan on implementing:",
                            ],
                            order=2,
                            mb="md",
                        ),
                        dmc.List(
                            [
                                dmc.ListItem(dmc.Text(item, c="dimmed"))
                                for item in not_working
                            ],
                            size="lg",
                            spacing="sm",
                        ),
                    ],
                    justify="space-between",
                    align="center",
                ),
            ],
            justify="space-between",
            align="center",
        )
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
    mb="xl",
    h="100%",
)

accession_card = dmc.Paper(
    [
        dmc.Title(
            html.Div(
                [
                    html.Span(
                        "starbase",
                        className="logo-text",
                    ),
                    " Accessions",
                ]
            ),
            order=2,
            mb="md",
        ),
        dmc.Text(
            [
                "To maintain data integrity, we employ an accessioning framework within ",
                html.Span(
                    "starbase", className="logo-text", style={"fontWeight": "bold"}
                ),
                ". ",
                html.Span(
                    "starbase", className="logo-text", style={"fontWeight": "bold"}
                ),
                " similar to NCBI-styled accessions. These accessions consist of a ",
                html.Span(
                    "unique six-digit numerical identifier",
                    style={"fontWeight": "bold"},
                ),
                " and provide a standardized nomenclature for keeping track of ",
                html.Span("Starships", style={"fontStyle": "italic"}),
                ". Accessions may represent 1) a ",
                html.Span("group", style={"fontWeight": "bold"}),
                " of highly similar Starship(s), that may be present in ",
                html.Span("multiple", style={"fontWeight": "bold"}),
                " different genomes or 2) a ",
                html.Span("single sequence", style={"fontWeight": "bold"}),
                " representing a Starship at a ",
                html.Span("single location", style={"fontWeight": "bold"}),
                " in a ",
                html.Span("single genome", style={"fontWeight": "bold"}),
                ".",
            ],
            c="dimmed",
        ),
        dmc.Space(h=20),
        dmc.Alert(
            "Classification/group accessions are SSA. Individual sequence accessions are SSB accessions.",
            color="blue",
            radius="md",
            mb="md",
        ),
        dmc.Space(h=20),
        dmc.Image(
            src="assets/images/accession_tag.png",
            fit="contain",
            radius="md",
            style={
                "maxWidth": "100%",
                "objectFit": "contain",
                "height": "auto",
            },
        ),
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
    h="100%",
)

# Database version card for the hero section
db_version_card = dmc.Paper(
    [
        dmc.Title("Database Version", order=2, mb="md"),
        dmc.Stack(
            [
                dmc.Group(
                    [
                        DashIconify(
                            icon="mdi:database",
                            width=24,
                            color="#868E96",
                        ),
                        dmc.Text(
                            "Current Version:",
                            size="sm",
                            c="dimmed",
                        ),
                    ],
                    gap="xs",
                    align="center",
                ),
                dmc.Title(
                    get_database_version(),
                    order=3,
                    style={"fontSize": "1.5rem", "color": "#228be6"},
                ),
                dmc.Divider(my="sm"),
                dmc.Group(
                    [
                        DashIconify(
                            icon="mdi:cog",
                            width=20,
                            color="#868E96",
                        ),
                        dmc.Text(
                            "Schema Version:",
                            size="sm",
                            c="dimmed",
                        ),
                    ],
                    gap="xs",
                    align="center",
                ),
                dmc.Text(
                    get_alembic_schema_version()[:8] + "..."
                    if len(get_alembic_schema_version()) > 8
                    else get_alembic_schema_version(),
                    size="sm",
                    c="dimmed",
                    style={"fontFamily": "monospace"},
                ),
            ],
            gap="xs",
        ),
    ],
    shadow="sm",
    p="xl",
    radius="md",
    withBorder=True,
    h="100%",
)


def create_hero_section():
    return dmc.Container(
        [
            dmc.Space(h=20),
            dmc.Center(title),
            dmc.Space(h=40),
            dmc.Grid(
                [
                    dmc.GridCol(starship_card, span={"base": 12, "md": 12}),
                ],
                gutter="xl",
                align="stretch",
            ),
        ],
        size="xl",
        flex=True,
    )


def create_features_section():
    return dmc.Container(
        [
            dmc.Grid(
                [
                    dmc.GridCol(
                        [working_features_card], span={"base": 12, "sm": 6, "lg": 4}
                    ),
                    dmc.GridCol([accession_card], span={"base": 12, "sm": 6, "lg": 4}),
                    dmc.GridCol([stats_section], span={"base": 12, "sm": 6, "lg": 4}),
                ],
                grow=True,
                gutter="xl",
                style={
                    "minHeight": "100%",
                    "alignItems": "stretch",  # Makes all grid items stretch to match heights
                    "display": "grid",
                    "gridAutoFlow": "dense",  # Fills in gaps in the grid
                },
            ),
        ],
        size="xl",
        py="xl",
        flex=True,
    )


def create_publications_section():
    return dmc.Container(
        [
            dmc.Paper(
                [
                    dmc.Title(
                        html.Div(
                            [
                                "Publications Characterizing ",
                                html.Span("Starships", style={"fontStyle": "italic"}),
                            ]
                        ),
                        order=2,
                        mb="md",
                    ),
                    html.Div(make_paper_table(), style={"width": "100%"}),
                ],
                py="xl",
                flex=True,
            )
        ],
        size="xl",
        py="xl",
        flex=True,
    )


stats = get_database_stats()

total_starships_info = dmc.Stack(
    [
        dmc.Group(
            [
                DashIconify(
                    icon="mdi:dna",
                    width=30,
                    color="#868E96",
                    style={"width": "clamp(24px, 4vw, 36px)"},
                ),
                dmc.Text(
                    [
                        "Total ",
                        html.Span(
                            "Starships",
                            style={"fontStyle": "italic"},
                        ),
                    ],
                    size="lg",
                    c="dimmed",
                    style={"fontSize": "clamp(0.9rem, 1.5vw, 1.3rem)"},
                ),
            ],
            gap="xs",
            align="center",
            justify="center",
            style={
                "flexDirection": "row",
                "flexWrap": "wrap",
                "justifyContent": "center",
                "alignItems": "center",
            },
        ),
        dmc.Title(
            f"{stats['total_starships']:,}",
            order=2,
            style={
                "fontSize": "clamp(1.8rem, 3vw, 2.5rem)",
                "textAlign": "center",
            },
        ),
    ],
    justify="center",
    gap="md",
    align="center",
)

curated_starships_info = dmc.Stack(
    [
        dmc.Group(
            [
                DashIconify(
                    icon="mdi:check-circle",
                    width=24,
                    color="#2F9E44",
                    style={"width": "clamp(20px, 3vw, 28px)"},
                ),
                dmc.Text(
                    "Curated",
                    size="sm",
                    c="dimmed",
                    style={"fontSize": "clamp(0.9rem, 1.5vw, 1.2rem)"},
                ),
            ],
            gap="xs",
            align="center",
            justify="center",
            style={
                "flexDirection": "row",
                "flexWrap": "wrap",
                "justifyContent": "center",
                "alignItems": "center",
            },
        ),
        dmc.Title(
            f"{stats['curated_starships']:,}",
            order=3,
            c="green",
            style={
                "fontSize": "clamp(1.4rem, 2vw, 1.8rem)",
                "textAlign": "center",
            },
        ),
    ],
    justify="center",
    gap="sm",
    align="center",
)

uncurated_starships_info = dmc.Stack(
    [
        dmc.Group(
            [
                DashIconify(
                    icon="mdi:clock",
                    width=24,
                    color="#E8590C",
                    style={"width": "clamp(20px, 3vw, 28px)"},
                ),
                dmc.Text(
                    "Uncurated",
                    size="sm",
                    c="dimmed",
                    style={"fontSize": "clamp(0.9rem, 1.5vw, 1.2rem)"},
                ),
            ],
            gap="xs",
            align="center",
            justify="center",
            style={
                "flexDirection": "row",
                "flexWrap": "wrap",
                "justifyContent": "center",
                "alignItems": "center",
            },
        ),
        dmc.Title(
            f"{stats['uncurated_starships']:,}",
            order=3,
            c="orange",
            style={
                "fontSize": "clamp(1.4rem, 2vw, 1.8rem)",
                "textAlign": "center",
            },
        ),
    ],
    justify="center",
    gap="sm",
    align="center",
)

species_info = dmc.Stack(
    [
        dmc.Group(
            [
                DashIconify(
                    icon="mdi:mushroom",
                    width=30,
                    color="#868E96",
                    style={"width": "clamp(24px, 4vw, 36px)"},
                ),
                dmc.Text(
                    "Species",
                    size="lg",
                    c="dimmed",
                    style={"fontSize": "clamp(0.9rem, 1.5vw, 1.3rem)"},
                ),
            ],
            gap="xs",
            align="center",
            justify="center",
            style={
                "flexDirection": "row",
                "flexWrap": "wrap",
                "justifyContent": "center",
                "alignItems": "center",
            },
        ),
        dmc.Title(
            f"{stats['species_count']:,}",
            order=2,
            style={
                "fontSize": "clamp(1.8rem, 3vw, 2.5rem)",
                "textAlign": "center",
            },
        ),
    ],
    justify="center",
    gap="md",
    align="center",
    style={"height": "100%"},
)

family_info = dmc.Stack(
    [
        dmc.Group(
            [
                DashIconify(
                    icon="mdi:family-tree",
                    width=30,
                    color="#868E96",
                    style={"width": "clamp(24px, 4vw, 36px)"},
                ),
                dmc.Text(
                    [
                        html.Span(
                            "Starship",
                            style={"fontStyle": "italic"},
                        ),
                        " Families",
                    ],
                    size="lg",
                    c="dimmed",
                    style={"fontSize": "clamp(0.9rem, 1.5vw, 1.3rem)"},
                ),
            ],
            gap="xs",
            align="center",
            justify="center",
            style={
                "flexDirection": "row",
                "flexWrap": "wrap",
                "justifyContent": "center",
                "alignItems": "center",
            },
        ),
        dmc.Title(
            f"{stats['family_count']:,}",
            order=2,
            style={
                "fontSize": "clamp(1.8rem, 3vw, 2.5rem)",
                "textAlign": "center",
            },
        ),
    ],
    justify="center",
    gap="md",
    align="center",
    style={"height": "100%"},
)


stats_section = dmc.Paper(
    [
        dmc.Title("Database Statistics", order=2, mb="xl"),
        dmc.SimpleGrid(
            [
                dmc.Stack(
                    [
                        dmc.Stack(
                            [
                                total_starships_info,
                                dmc.Group(
                                    [
                                        curated_starships_info,
                                        uncurated_starships_info,
                                    ],
                                    gap="lg",
                                    grow=True,
                                    align="center",
                                    justify="space-around",
                                ),
                                dmc.Divider(my="md", color="gray.3"),
                            ],
                            gap="md",
                        ),
                    ],
                    gap="xl",
                    justify="center",
                    style={"height": "100%"},
                ),
                dmc.Stack(
                    [
                        species_info,
                        dmc.Divider(my="md", color="gray.3"),
                        family_info,
                    ],
                    gap="xl",
                    justify="center",
                    style={"height": "100%"},
                ),
            ],
            cols={"base": 1, "sm": 2, "md": 2, "lg": 1},
            spacing="xl",
            style={"minHeight": "100%", "alignItems": "stretch"},
        ),
    ],
    p="xl",
    radius="md",
    shadow="sm",
    withBorder=True,
)


layout = dmc.MantineProvider(
    [
        dcc.Location(id="url", refresh=False),
        create_hero_section(),
        dmc.Space(h=40),
        create_features_section(),
        dmc.Space(h=40),
        create_publications_section(),
    ],
    theme={
        "colorScheme": "light",
        "primaryColor": "indigo",
        "components": {
            "Container": {"defaultProps": {"size": "xl"}},
            "Title": {"defaultProps": {"color": "indigo"}},
        },
    },
)


@callback(
    [
        Output("toggle-paper-view", "leftSection"),
        Output("desktop-table", "style"),
        Output("mobile-cards", "style"),
    ],
    Input("toggle-paper-view", "n_clicks"),
    State("toggle-paper-view", "leftSection"),
    prevent_initial_call=True,
)
def toggle_paper_view(n_clicks, current_icon):
    if n_clicks is None:
        raise PreventUpdate

    is_table = current_icon["icon"] == "tabler:layout-list"

    return (
        DashIconify(icon="tabler:layout-cards" if is_table else "tabler:layout-list"),
        {"display": "none" if is_table else "block"},
        {"display": "block" if is_table else "none"},
    )

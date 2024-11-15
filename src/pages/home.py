import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html
from dash_iconify import DashIconify
import logging

from src.components.tables import make_paper_table
from src.components.callbacks import download_ships_card
from src.components.sql_engine import sql_connected
from src.components.sql_queries import get_database_stats

logger = logging.getLogger(__name__)

dash.register_page(__name__, title="starbase Home", name="Home", path="/")

logger.info(f"sql_connected? {sql_connected}")

title = dmc.Title(
    [
        html.Span(
            "starbase: ",
            className="logo-text",
        ),
        "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
    ],
    className="text-center",
    style={"paddingTop": "20px"},
    # className="text-center text-custom text-custom-xl",
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
            disabled=not sql_connected,
        ),
        href=href,
        style={"textDecoration": "none"},  # Remove underline from link
    )

working_buttons = [
    create_feature_button("Search Database", "/search", "mdi:database-search"),
    create_feature_button("BLAST Search", "/blast", "mdi:dna"),
    create_feature_button("Browse Wiki", "/wiki", "mdi:book-open-variant"),
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
not_working_ul = html.Ul(
    [
        html.Li(
            item,
            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
        )
        for item in not_working
    ],
)

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
                                        ],
                                    ),
            ], size="lg", c="dimmed"),
        ], span=7),
        dmc.GridCol([
            dmc.Image(
                src="assets/images/starship-model.png",
                fit="contain",
                radius="md"
            )
        ], span=5),
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
                dmc.ListItem(item) for item in not_working
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
            ], span=6),
            dmc.GridCol(download_ships_card, span=6)
        ], gutter="xl"),
    ], size="xl", py="xl", flex=True)

def create_publications_section():
    return dmc.Container([
        dmc.Title("Manuscripts Characterizing Starships", order=2, mb="xl"),
        dmc.Center(make_paper_table()),
    ], size="xl", py="xl",flex=True)

def create_stats_section():
    # Get stats from database
    stats = get_database_stats() if sql_connected else {
        "total_starships": "—",
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
                dmc.Title(
                    f"{stats['total_starships']:,}", 
                    order=3
                ),
            ], p="xl", radius="md", withBorder=True),
            dmc.Paper([
                dmc.Text("Species", size="lg", c="dimmed"),
                dmc.Title(
                    f"{stats['species_count']:,}", 
                    order=3
                ),
            ], p="xl", radius="md", withBorder=True),
            dmc.Paper([
                dmc.Text("Starship Families", size="lg", c="dimmed"),
                dmc.Title(
                    f"{stats['family_count']:,}", 
                    order=3
                ),
            ], p="xl", radius="md", withBorder=True),
        ], cols=3),
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
        color="red",
        style={"position": "fixed", "top": 20, "right": 20}
    ) if not sql_connected else None,
], 
theme={
    "colorScheme": "light",
    "primaryColor": "indigo",
    "components": {
        "Container": {"defaultProps": {"size": "xl"}},
        "Title": {"defaultProps": {"color": "indigo"}},
    }
})
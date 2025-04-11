import dash_bootstrap_components as dbc
from dash import html

from src.pages import (
    HOME_URL,
    DL_URL,
    WIKI_URL,
    PGV_URL,
    BLAST_URL,
    SUBMIT_URL,
    ABOUT_URL,
)


def navmenu(buttons_disabled=False):
    navbar_content = [
        dbc.NavItem(
            dbc.NavLink(
                ["Download ", html.Span("starbase", className="logo-text")],
                href=DL_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "Wiki",
                href=WIKI_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "Starship Viewer",
                href=PGV_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "BLAST",
                href=BLAST_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "Submit",
                href=SUBMIT_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                "About",
                href=ABOUT_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
    ]

    navbar = dbc.NavbarSimple(
        children=navbar_content,
        brand=dbc.NavbarBrand(
            html.Span("starbase", className="logo-text"),
            className="ms-2",
            href=HOME_URL,
        ),
        brand_href=HOME_URL,
        brand_style={
            "align-items": "center",
            "justify-content": "center",
            "textAlign": "left",
        },
        color="primary",
        dark=True,
    )

    return html.Div(
        [navbar],
    )

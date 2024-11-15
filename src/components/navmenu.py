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
    METRICS_URL,
)


def navmenu(buttons_disabled=False):
    navbar_content = [
        dbc.NavItem(
            dbc.NavLink(
                html.P(
                    ["Download ", html.Span("starbase", className="logo-text")],
                    className="nav-item-text",
                ),
                href=DL_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                html.P(
                    [
                        "Wiki",
                    ],
                    className="nav-item-text",
                ),
                href=WIKI_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                html.P(
                    ["Starship Viewer"],
                    className="nav-item-text",
                ),
                href=PGV_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                html.P(
                    ["BLAST"],
                    className="nav-item-text",
                ),
                href=BLAST_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                html.P(
                    ["Submit"],
                    className="nav-item-text",
                ),
                href=SUBMIT_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
        dbc.NavItem(
            dbc.NavLink(
                html.P(["About"], className="nav-item-text"),
                href=ABOUT_URL,
                active="exact",
                className="nav-item-link",
                disabled=buttons_disabled,
            )
        ),
    ]



    navbar = dbc.NavbarSimple(
        children=navbar_content,
        brand=dbc.NavbarBrand(html.Span("starbase", className="logo-text"), className="ms-2"),
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

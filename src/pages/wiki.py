import dash
import dash_bootstrap_components as dbc

from src.data.joined_ships import df
from src.components.shipTable import make_table

dash.register_page(__name__)

df_sub = df[
    [
        "starshipID",
        "captain_superfamily",
        "starship_family",
        "starship_navis",
        "starship_haplotype",
        "genus",
        "species",
    ]
]


layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Stack(
            [
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            style={"padding": "20px"},
                            children=[make_table(df_sub)],
                        )
                    ],
                ),
            ],
            gap=3,
        ),
    ],
)

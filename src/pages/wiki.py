import warnings

warnings.filterwarnings("ignore")
import dash
from dash import Input, Output, State, html, dcc, callback, callback_context
from dash.dependencies import Output, Input, State
import plotly.express as px
import dash_bootstrap_components as dbc

from src.data.joined_ships import df
from src.components.shipTable import make_table
from src.utils.sunburst_plot import create_sunburst_plot
from src.utils.logomaker import make_logo

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

# Extract unique values from the column
unique_categories = df["starship_family"].dropna().unique()


def create_cards(category):
    if category == "nan":
        button = ""
        card = ""
    else:
        filtered_df = df[df["starship_family"] == category]
        n_ships = len(filtered_df["checksum"].dropna().unique())
        min_size = min(filtered_df["size"])
        max_size = max(filtered_df["size"])
        upTIRs = filtered_df["upTIR"].dropna().tolist()
        downTIRs = filtered_df["downTIR"].dropna().tolist()
        sunburst = create_sunburst_plot(
            filtered_df,
            ["genus", "species"],
            f"Genus/Species Distribution for {category}",
        )

        uplogo = make_logo(upTIRs)
        downlogo = make_logo(downTIRs)

        button = dbc.Button(
            f"{category}",
            id={"type": "collapse-button", "index": category},
            className="d-flex justify-content-start col-3 text-center",
            color="primary",
            n_clicks=0,
        )

        # Build card body content dynamically
        card_body_contents = [
            html.H5(f"Total Number of Starships in {category}: {n_ships}"),
            html.H5(f"Maximum Starship Size (bp): {max_size}"),
            html.H5(f"Minimum Starship Size (bp): {min_size}"),
            dcc.Graph(figure=sunburst),
        ]

        if uplogo:
            card_body_contents.append(
                dbc.Row(
                    [
                        html.H5(f"Sequence logo of upstream TIRs in {category}"),
                        html.Img(
                            src=f"data:image/png;base64,{uplogo}",
                            style={"width": "50%"},
                        ),
                    ]
                )
            )

        if downlogo:
            card_body_contents.append(
                dbc.Row(
                    [
                        html.H5(f"Sequence logo of downstream TIRs in {category}"),
                        html.Img(
                            src=f"data:image/png;base64,{downlogo}",
                            style={"width": "50%"},
                        ),
                    ]
                )
            )

        card = dbc.Collapse(
            dbc.Card(dbc.CardBody(card_body_contents)),
            id={"type": "collapse", "index": category},
            is_open=False,
        )

    return [button, card]


layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Row(
            justify="center",
            align="middle",
            style={"paddingTop": "20px"},
            children=[
                dbc.Col(
                    lg=6,
                    sm=12,
                    children=[
                        html.H1(
                            [
                                html.Span(
                                    "starbase",
                                    className="logo-text",
                                ),
                                " Wiki",
                            ]
                        ),
                        html.H2(
                            "Summary and characteristics of each Starship superfamily"
                        ),
                    ],
                )
            ],
        ),
        dbc.Row(
            justify="center",
            align="start",
            style={"paddingTop": "20px"},
            children=[
                dbc.Col(
                    lg=6,
                    sm=12,
                    children=[
                        dbc.Stack(
                            children=[
                                component
                                for category in unique_categories
                                for component in create_cards(category)
                            ],
                            direction="vertical",
                            gap=3,
                        )
                    ],
                ),
            ],
        ),
    ],
)


# Define callback using MATCH for dynamic callbacks
@callback(
    Output({"type": "collapse", "index": dash.dependencies.ALL}, "is_open"),
    [Input({"type": "collapse-button", "index": dash.dependencies.ALL}, "n_clicks")],
    [State({"type": "collapse", "index": dash.dependencies.ALL}, "is_open")],
)
def toggle_collapse(n_clicks, is_open):
    # Prevents callback from firing if the app is just loading
    ctx = callback_context
    if not ctx.triggered:
        return [False] * len(is_open)

    triggered_button = ctx.triggered[0]["prop_id"].split(".")[0]
    triggered_index = eval(triggered_button)["index"]  # Extract index from button ID

    return [
        not open if index == triggered_index else open
        for index, open in zip(unique_categories, is_open)
    ]

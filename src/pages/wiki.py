import warnings

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc

from src.data.joined_ships import df
from src.utils.plot_utils import create_sunburst_plot
from src.utils.plot_utils import make_logo

dash.register_page(__name__)


def create_accordion_item(category):
    if category == "nan":
        return None
    else:
        filtered_df = df[df["starship_family"] == category]
        n_ships = len(filtered_df["checksum"].dropna().unique())
        min_size = min(filtered_df["size"].dropna())
        max_size = max(filtered_df["size"].dropna())
        upTIRs = filtered_df["upTIR"].dropna().tolist()
        downTIRs = filtered_df["downTIR"].dropna().tolist()
        sunburst = create_sunburst_plot(
            df=filtered_df,
            groups=["genus", "species"],
            title=f"Genus/Species Distribution for {category}",
        )

        uplogo = make_logo(upTIRs)
        downlogo = make_logo(downTIRs)

        accordion_content = [
            html.H5(f"Total Number of Starships in {category}: {n_ships}"),
            html.H5(f"Maximum Starship Size (bp): {max_size}"),
            html.H5(f"Minimum Starship Size (bp): {min_size}"),
            dcc.Graph(figure=sunburst),
        ]

        if uplogo:
            accordion_content.append(
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
            accordion_content.append(
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

        return dbc.AccordionItem(
            title=category,
            children=[dbc.CardBody(accordion_content)],
            item_id=category,
        )


def create_accordion(df):
    unique_categories = df["starship_family"].dropna().unique().tolist()
    assert isinstance(unique_categories, list), "unique_categories must be a list"

    accordion_items = [
        create_accordion_item(category)
        for category in unique_categories
        if category != "nan"
    ]

    return dbc.Accordion(
        children=accordion_items,
        id="category-accordion",
        always_open=False,
    )


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
                        html.H2("Summary and characteristics of each Starship family"),
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
                            children=create_accordion(df),
                            direction="vertical",
                            gap=3,
                        )
                    ],
                ),
            ],
        ),
    ],
)

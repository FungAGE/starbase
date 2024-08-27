import warnings

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State
import dash_bootstrap_components as dbc

from src.utils.plot_utils import create_sunburst_plot
from src.utils.plot_utils import make_logo

import pandas as pd

dash.register_page(__name__)


def create_accordion_item(df, papers, category):
    if category == "nan":
        return None
    else:
        filtered_meta_df = df[df["familyName"] == category]
        n_ships = len(filtered_meta_df["checksum"].dropna().unique())
        min_size = min(filtered_meta_df["size"].dropna())
        max_size = max(filtered_meta_df["size"].dropna())
        upDRs = filtered_meta_df["upDR"].dropna().tolist()
        downDRs = filtered_meta_df["downDR"].dropna().tolist()
        filtered_papers_df = papers[papers["familyName"] == category]
        type_element_reference = filtered_papers_df["type_element_reference"]
        sunburst = create_sunburst_plot(df=filtered_meta_df, type="tax")

        uplogo = make_logo(upDRs)
        downlogo = make_logo(downDRs)

        accordion_content = [
            html.H5(f"Total Number of Starships in {category}: {n_ships}"),
            html.H5(f"Maximum Starship Size (bp): {max_size}"),
            html.H5(f"Minimum Starship Size (bp): {min_size}"),
        ]

        if type_element_reference is not None:

            link = filtered_papers_df["Url"]
            paper_link = html.Link(
                type_element_reference,
                href=link,
            )

            accordion_content.append(
                html.H5(["Reference for defining type element in family: ", paper_link])
            )

        accordion_content.append(dcc.Graph(figure=sunburst))

        if uplogo:
            uplogo_img = dbc.Col(
                lg=6,
                sm=12,
                children=[
                    html.H5(f"Sequence logo of upstream DRs in {category}"),
                    html.Img(
                        src=f"data:image/png;base64,{uplogo}",
                        style={"width": "100%"},
                    ),
                ],
            )
        if downlogo:
            downlogo_img = dbc.Col(
                lg=6,
                sm=12,
                children=[
                    html.H5(f"Sequence logo of downstream DRs in {category}"),
                    html.Img(
                        src=f"data:image/png;base64,{downlogo}",
                        style={"width": "100%"},
                    ),
                ],
            )
        if uplogo and downlogo:
            accordion_content.append(dbc.Row([uplogo_img, downlogo_img]))
        elif uplogo and not downlogo:
            accordion_content.append(dbc.Row([uplogo_img]))
        elif downlogo and not uplogo:
            accordion_content.append(dbc.Row([downlogo_img]))

        return dbc.AccordionItem(
            title=category,
            children=[dbc.CardBody(accordion_content)],
            item_id=category,
        )


def create_accordion(df, papers):
    unique_categories = df["familyName"].dropna().unique().tolist()
    assert isinstance(unique_categories, list), "unique_categories must be a list"

    accordion_items = [
        create_accordion_item(df, papers, category)
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
        dcc.Location(id="url", refresh=False),
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
                        dcc.Loading(
                            id="wiki-loading",
                            type="default",
                            children=[
                                dbc.Stack(
                                    children=html.Div(id="accordion"),
                                    direction="vertical",
                                    gap=3,
                                )
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@callback(
    [Output("accordion", "children"), Output("accordion", "active_item")],
    [
        Input("url", "search"),
        Input("joined-ships", "data"),
        Input("paper-cache", "data"),
    ],
    [
        State("accordion", "active_item"),
    ],
)
def load_data(search, cached_meta, cached_papers, active_item):
    accordion = None
    query_param = active_item

    if cached_meta and cached_papers:
        initial_df = pd.DataFrame(cached_meta)
        paper_df = pd.DataFrame(cached_papers)
        accordion = create_accordion(initial_df, paper_df)
        if search:
            query_param = search.split("=")[-1]
    return accordion, query_param

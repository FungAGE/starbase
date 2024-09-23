import warnings

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html, dash_table
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import pandas as pd
import os

from src.components.sqlite import engine
from src.utils.plot_utils import create_sunburst_plot
from src.utils.plot_utils import make_logo


dash.register_page(__name__)


def create_accordion_item(df, papers, category):
    if category == "nan":
        return None
    else:
        filtered_meta_df = df[df["familyName"] == category]
        n_ships = len(filtered_meta_df["accession_tag"].dropna().unique())
        min_size = min(filtered_meta_df["size"].dropna())
        max_size = max(filtered_meta_df["size"].dropna())
        upDRs = filtered_meta_df["upDR"].dropna().tolist()
        downDRs = filtered_meta_df["downDR"].dropna().tolist()
        filtered_papers_df = papers[papers["familyName"] == category]
        type_element_reference = (
            filtered_papers_df["type_element_reference"].dropna().unique().astype(str)
        )
        sunburst = create_sunburst_plot(df=filtered_meta_df, type="tax")

        uplogo_img_path = f"assets/images/DR/{category}-upDR.png"
        if not os.path.exists(uplogo_img_path):
            uplogo_img_path = make_logo(upDRs, uplogo_img_path)
        uplogo_img = dbc.Col(
            lg=6,
            sm=12,
            children=[
                dmc.Center(html.H5(f"Sequence logo of upstream DRs in {category}")),
                dmc.Center(
                    html.Img(
                        src=uplogo_img_path,
                        style={"width": "100%"},
                    ),
                ),
            ],
        )

        downlogo_img_path = f"assets/images/DR/{category}-downDR.png"
        if not os.path.exists(downlogo_img_path):
            downlogo_img_path = make_logo(downDRs, downlogo_img_path)
        downlogo_img = dbc.Col(
            lg=6,
            sm=12,
            children=[
                html.H5(f"Sequence logo of downstream DRs in {category}"),
                html.Img(
                    src=downlogo_img_path,
                    style={"width": "100%"},
                ),
            ],
        )

        accordion_content = [
            dash_table.DataTable(
                columns=[
                    {"name": "Metric", "id": "metric"},
                    {"name": "Value", "id": "value"},
                ],
                data=[
                    {
                        "metric": "Total Number of Starships in {}".format(category),
                        "value": f"{n_ships:,.0f}",
                    },
                    {
                        "metric": "Maximum Starship Size (bp)",
                        "value": f"{max_size:,.0f}",
                    },
                    {
                        "metric": "Minimum Starship Size (bp)",
                        "value": f"{min_size:,.0f}",
                    },
                ],
                style_table={"overflowX": "auto", "maxWidth": "500px"},
                style_cell={"textAlign": "left"},
                style_header={"fontWeight": "bold"},
                # You can customize more styles as needed
            )
        ]
        if type_element_reference.size > 0:
            # Get unique URLs and convert to strings
            link = filtered_papers_df["Url"].dropna().unique().astype(str)

            # Create a list of links
            paper_links = [
                html.A(ref, href=url, target="_blank")  # Use html.A for links
                for ref, url in zip(type_element_reference, link)
            ]

            # Combine links into a single H5 element
            accordion_content.append(
                html.H5(
                    ["Reference for defining type element in family: ", *paper_links]
                )
            )

        accordion_content.append(dcc.Graph(figure=sunburst))

        if uplogo_img and downlogo_img:
            accordion_content.append(dbc.Row([uplogo_img, downlogo_img]))
        elif uplogo_img and not downlogo_img:
            accordion_content.append(dbc.Row([uplogo_img]))
        elif downlogo_img and not uplogo_img:
            accordion_content.append(dbc.Row([downlogo_img]))

        return dbc.AccordionItem(
            title=category,
            children=[
                dbc.CardBody(dbc.Stack(accordion_content, direction="vertical", gap=5))
            ],
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
        always_open=True,
        active_item="Voyager",
    )


def load_data():
    # Query your database to fetch the relevant data
    meta_query = """
    SELECT j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, a.accession_tag, t."order", t.family
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN accessions a ON j.ship_id = a.id
    """
    meta_df = pd.read_sql_query(meta_query, engine)

    paper_query = """
    SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, p.shortCitation, f.familyName, f.type_element_reference
    FROM papers p
    JOIN family_names f ON p.shortCitation = f.type_element_reference
    """
    paper_df = pd.read_sql_query(paper_query, engine)

    accordion = create_accordion(meta_df, paper_df)

    return accordion


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
                            "Wiki - Summary and characteristics for each Starship family",
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
                        dcc.Loading(
                            id="wiki-loading",
                            type="circle",
                            children=[load_data()],
                        ),
                    ],
                ),
            ],
        ),
    ],
)

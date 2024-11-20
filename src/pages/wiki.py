import warnings

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html, dash_table, callback, no_update
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import pandas as pd
import pickle
import os
import logging

from sqlalchemy.exc import SQLAlchemyError

from src.components.cache import cache
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import (
    fetch_meta_data,
    cache_sunburst_plot,
    fetch_paper_data,
)
from src.components.tables import make_ship_table, make_wiki_table
from src.utils.plot_utils import make_logo
from src.utils.seq_utils import clean_contigIDs
from src.components.callbacks import create_accession_modal

dash.register_page(__name__)

logger = logging.getLogger(__name__)


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

        if len(downDRs) > 10:
            uplogo_img_path = f"assets/images/DR/{category}-upDR.png"

            if not os.path.exists(uplogo_img_path):
                uplogo_img_path = make_logo(upDRs, uplogo_img_path, type="up")
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
        else:
            uplogo_img = None

        if len(downDRs) > 10:
            downlogo_img_path = f"assets/images/DR/{category}-downDR.png"
            if not os.path.exists(downlogo_img_path):
                downlogo_img_path = make_logo(downDRs, downlogo_img_path, type="down")
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
        else:
            downlogo_img = None

        accordion_content = [
            make_wiki_table(category, n_ships, max_size, min_size)

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

modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle(id="wiki-modal-title")),
        dbc.ModalBody(id="wiki-modal-content"),
    ],
    id="wiki-modal",
    is_open=False,
)

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="meta-data"),
        dcc.Store(id="paper-data"),
        dcc.Store(id="active-item-cache"),
        
        # Header Section
        dmc.Space(h=20),
        dmc.Paper(
            children=[
                dmc.Title(
                    "Starship Family Wiki",
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    "Explore characteristics and distributions of different Starship families",
                    c="dimmed",
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
        ),
        
        # Main Content Grid
        dmc.Grid(
            children=[
                # Left Column - Family Selector
                dmc.GridCol(
                    span={"lg": 4, "md": 12},
                    children=[
                        dmc.Paper(
                            children=[
                                dmc.Title("Starship Families", order=2, mb="md"),
                                dcc.Loading(
                                    id="wiki-loading",
                                    type="circle",
                                    children=html.Div(
                                        id="accordion",
                                        children=dbc.Accordion(id="category-accordion"),
                                    ),
                                ),
                            ],
                            p="md",
                            radius="md",
                            withBorder=True,
                            style={"height": "calc(100vh - 200px)", "overflowY": "auto"},
                        ),
                    ],
                ),
                
                # Right Column - Family Details
                dmc.GridCol(
                    span={"lg": 8, "md": 12},
                    children=[                       
                        # Content for selected family
                        dmc.Paper(
                            children=[
                                dcc.Loading(
                                    id="sidebar-loading",
                                    type="circle",
                                    children=dmc.Stack([
                                        # Sunburst Plot
                                        dmc.Paper(
                                            children=[
                                                html.Div(id="sidebar-title"),
                                                html.Div(
                                                    id="sidebar",
                                                    style={
                                                        "width": "100%",
                                                        "minHeight": "400px",
                                                    },
                                                ),
                                            ],
                                            p="md",
                                            radius="md",
                                            withBorder=True,
                                        ),
                                    ], gap=3),
                                ),
                            ],
                            p="md",
                            radius="md",
                            withBorder=True,
                        ),
                    ],
                ),
            ],
            gutter="xl",
        ),
        modal,  # Keep the modal at the root level
    ],
)


@cache.memoize()
@callback(Output("meta-data", "data"), Input("url", "href"))
def load_meta_data(url):
    if url:
        try:
            logger.debug("Loading metadata from database.")
            meta_df = load_from_cache("meta_data")
            if meta_df is None:
                meta_df = fetch_meta_data()
            logger.info(f"Metadata query returned {len(meta_df)} rows.")

            # Clean 'contigID' column if present
            if "contigID" in meta_df.columns:
                meta_df["contigID"] = meta_df["contigID"].apply(clean_contigIDs)
                logger.debug("Cleaned contigID column in the metadata.")

            return meta_df.to_dict("records")
        except SQLAlchemyError as e:
            logger.error(
                "Error occurred while executing the metadata query.", exc_info=True
            )
            return []

    logger.warning("No URL provided for metadata loading.")
    raise PreventUpdate  # Skip the update if no URL is provided


# Callback to load paper data
@cache.memoize()
@callback(Output("paper-data", "data"), Input("url", "href"))
def load_paper_data(url):
    if url:
        try:
            logger.debug("Loading paper data from database.")
            paper_df = load_from_cache("paper_data")
            if paper_df is None:
                paper_df = fetch_paper_data()
            logger.info(f"Paper data query returned {len(paper_df)} rows.")
            return paper_df.to_dict("records")
        except SQLAlchemyError as e:
            logger.error(
                "Error occurred while executing the paper data query.", exc_info=True
            )
            return []

    logger.warning("No URL provided for paper data loading.")
    raise PreventUpdate


# Callback to create accordion layout
@callback(
    Output("accordion", "children"),
    [Input("meta-data", "data"), Input("paper-data", "data")],
)
def create_accordion(cached_meta, cached_papers):
    if cached_meta is None or cached_papers is None:
        logger.error("One or more inputs are None: meta-data or paper-data.")
        raise PreventUpdate

    df = pd.DataFrame(cached_meta)
    papers = pd.DataFrame(cached_papers)

    assert isinstance(
        df, pd.DataFrame
    ), f"Expected df to be a DataFrame, but got {type(df)}."
    assert isinstance(
        papers, pd.DataFrame
    ), f"Expected papers to be a DataFrame, but got {type(papers)}."

    logger.debug("Creating accordion based on metadata and paper data.")

    try:
        unique_categories = df["familyName"].dropna().unique().tolist()
        logger.info(
            f"Found {len(unique_categories)} unique categories for the accordion."
        )
        assert isinstance(unique_categories, list), "unique_categories must be a list"

        accordion_items = [
            create_accordion_item(df, papers, category)
            for category in unique_categories
            if category != "nan"
        ]
        logger.debug("Accordion items created successfully.")

        return dbc.Accordion(
            children=accordion_items,
            id="category-accordion",
            always_open=False,
            active_item="Prometheus",
        )
    except Exception as e:
        logger.error("Error occurred while creating the accordion.", exc_info=True)
        raise


# Callback to create sidebar content
@callback(
    [
        Output("sidebar", "children"),
        Output("sidebar-title", "children"),
        Output("active-item-cache", "value"),
    ],
    Input("category-accordion", "active_item"),
    State("meta-data", "data"),
)
def create_sidebar(active_item, cached_meta):
    if active_item is None or cached_meta is None:
        logger.warning("No active item or cached meta data provided.")
        raise PreventUpdate  # Prevent the callback from running if the inputs are invalid

    logger.debug(f"Creating sidebar for active item: {active_item}")
    df = pd.DataFrame(cached_meta)

    try:
        title = dmc.Title(f"Taxonomy and Genomes for Starships in {active_item}", order=2, mb="md")
        filtered_meta_df = df[df["familyName"] == active_item].sort_values(
            by="accession_tag", ascending=False
        )
        logger.info(
            f"Filtered metadata for {active_item} contains {len(filtered_meta_df)} rows."
        )

        sunburst_figure = load_from_cache(f"sunburst_{active_item}")
        if sunburst_figure is None:
            sunburst_figure = cache_sunburst_plot(
                family=active_item, df=filtered_meta_df
            )

        fig = dcc.Graph(
            figure=sunburst_figure, style={"width": "100%", "height": "100%"}
        )

        table_columns = [
            {
                "name": "Accession",
                "id": "accession_tag",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
            # {
            #     "name": "Species",
            #     "id": "species",
            #     "deletable": False,
            #     "selectable": False,
            # },
            # {
            #     "name": "Contig ID",
            #     "id": "contigID",
            #     "deletable": False,
            #     "selectable": False,
            # },
            # {
            #     "name": "Start Position in Genome",
            #     "id": "elementBegin",
            #     "deletable": False,
            #     "selectable": False,
            # },
            # {
            #     "name": "End Position in Genome",
            #     "id": "elementEnd",
            #     "deletable": False,
            #     "selectable": False,
            # },
            {
                "name": "Element Length",
                "id": "size",
                "deletable": False,
                "selectable": False,
            },
        ]
        table = make_ship_table(
            filtered_meta_df, id="wiki-table", columns=table_columns, pg_sz=15
        )
        logger.debug("Table for sidebar created successfully.")

        output = dbc.Stack(children=[fig, table, modal], direction="vertical", gap=5)
        return output, title, active_item
    except Exception as e:
        logger.error(
            f"Error occurred while creating the sidebar for {active_item}.",
            exc_info=True,
        )
        raise


@callback(
    Output("wiki-modal", "is_open"),
    Output("wiki-modal-content", "children"),
    Output("wiki-modal-title", "children"),
    Output("wiki-table", "active_cell"),
    Input("wiki-table", "active_cell"),
    State("wiki-modal", "is_open"),
    State("wiki-table", "data"),
)
def toggle_modal(active_cell, is_open, table_data):
    if active_cell:
        row = active_cell["row"]
        row_data = table_data[row]
        modal_content, modal_title = create_accession_modal(row_data["accession_tag"])

        # Return the modal open, content, title, and reset active_cell
        return True, modal_content, modal_title, None

    # Keep the modal closed if there's no active cell
    return is_open, no_update, no_update, no_update

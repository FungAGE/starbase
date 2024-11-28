import warnings

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html, dash_table, callback, no_update
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

from dash.long_callback import DiskcacheLongCallbackManager
from functools import lru_cache

import pandas as pd
import pickle
import os
import logging

from sqlalchemy.exc import SQLAlchemyError

from src.components.cache import cache
from src.components.sql_manager import load_from_cache, fetch_meta_data
from src.components.sql_manager import (
    fetch_meta_data,
    cache_sunburst_plot,
    fetch_paper_data,
)
from src.components.tables import make_ship_table, make_wiki_table
from src.utils.plot_utils import make_logo
from src.utils.seq_utils import clean_contigIDs
from src.components.callbacks import create_accession_modal, create_modal_callback

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
                    dmc.Center(html.H5("Upstream DRs")),
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
                    dmc.Center(html.H5("Downstream DRs")),
                    dmc.Center(
                        html.Img(
                            src=downlogo_img_path,
                            style={"width": "100%"},
                        ),
                    ),
                ],
            )
        else:
            downlogo_img = None

        accordion_content = [
            make_wiki_table(n_ships, max_size, min_size)

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

modal = dmc.Modal(
    id="wiki-modal",
    opened=False,
    centered=True,
    overlayProps={"blur": 3},
    size="lg",
    children=[
        dmc.Title(id="wiki-modal-title", order=3),
        dmc.Space(h="md"),
        html.Div(id="wiki-modal-content"),
    ],
)

def load_initial_data():
    """Load initial data for the page"""
    meta_data = load_from_cache("meta_data")
    if meta_data is None:
        # Fallback to fetching directly if not in cache
        meta_data = fetch_meta_data()
        dcc.Store(id="filtered-meta-data"),
    # Convert DataFrame to dictionary before storing
    if isinstance(meta_data, pd.DataFrame):
        return meta_data.to_dict('records')
    return meta_data

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="meta-data", data=load_initial_data()),
        dcc.Store(id="filtered-meta-data"),
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
                    "Search and explore the characteristics of different Starship families",
                    c="dimmed",
                    size="lg",
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
        ),
        
        dmc.Grid(
            children=[
                # Left Column - Search Section
                dmc.GridCol(
                    span={"lg": 4, "md": 12},
                    children=[

                        dmc.Paper(
                            children=[
                                dmc.Title("Search Starships", order=3, mb="md"),
                                dmc.Grid(
                                    children=[
                                        # Taxonomy Search
                                        dmc.GridCol(
                                            span={"lg": 3, "md": 6, "sm": 12},
                                            children=[
                                                dmc.MultiSelect(
                                                    id="taxonomy-search",
                                                    label="Taxonomy",
                                                    placeholder="Search by taxonomy...",
                                                    searchable=True,
                                                    clearable=True,
                                                    nothingFoundMessage="No options found",
                                                    data=[],  # Will be populated by callback
                                                    value=[],  # Initialize with empty list
                                                ),
                                            ],
                                        ),
                                        # Family Search
                                        dmc.GridCol(
                                            span={"lg": 3, "md": 6, "sm": 12},
                                            children=[
                                                dmc.MultiSelect(
                                                    id="family-search",
                                                    label="Starship Family",
                                                    placeholder="Search by family...",
                                                    searchable=True,
                                                    clearable=True,
                                                    nothingFoundMessage="No options found",
                                                    data=[],  # Will be populated by callback
                                                    value=[],  # Initialize with empty list
                                                ),
                                            ],
                                        ),
                                        # Navis Search
                                        dmc.GridCol(
                                            span={"lg": 3, "md": 6, "sm": 12},
                                            children=[
                                                dmc.MultiSelect(
                                                    id="navis-search",
                                                    label="Navis",
                                                    placeholder="Search by navis...",
                                                    searchable=True,
                                                    clearable=True,
                                                    nothingFoundMessage="No options found",
                                                    data=[],  # Will be populated by callback
                                                    value=[],  # Initialize with empty list
                                                ),
                                            ],
                                        ),
                                        # Haplotype Search
                                        dmc.GridCol(
                                            span={"lg": 3, "md": 6, "sm": 12},
                                            children=[
                                                dmc.MultiSelect(
                                                    id="haplotype-search",
                                                    label="Haplotype",
                                                    placeholder="Search by haplotype...",
                                                    searchable=True,
                                                    clearable=True,
                                                    nothingFoundMessage="No options found",
                                                    data=[],  # Will be populated by callback
                                                    value=[],  # Initialize with empty list
                                                ),
                                            ],
                                        ),
                                    ],
                                    gutter="xl",
                                ),
                                # Search Actions
                                dmc.Group(
                                    align="right",
                                    mt="md",
                                    children=[
                                        dmc.Button(
                                            "Reset",
                                            id="reset-search",
                                            variant="outline",
                                            leftSection=DashIconify(icon="tabler:refresh"),
                                        ),
                                        dmc.Button(
                                            "Search",
                                            id="apply-search",
                                            variant="filled",
                                            leftSection=DashIconify(icon="tabler:search"),
                                        ),
                                    ],
                                ),
                            ],
                            p="xl",
                            radius="md",
                            withBorder=True,
                            mb="xl",
                        ),
                    ]),
                # Right Column - Search Results
                dmc.GridCol(
                    span={"lg": 8, "md": 12},
                    children=html.Div(id="search-results")
                )
            ],
        ),
        dmc.Space(h=20),
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
                                        html.Div(id="sidebar-title"),
                                        modal,
                                        html.Div(
                                            id="sidebar",
                                            style={
                                                "width": "100%",
                                                "minHeight": "400px",
                                                "maxHeight": "calc(100vh - 300px)",
                                                "overflow": "auto",
                                            },
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


@callback(
    Output("search-results", "children"),
    Input("filtered-meta-data", "data"),
    State("meta-data", "data"),
)
def create_search_results(filtered_meta, cached_meta):
    # Use filtered data if available, otherwise use original data
    data_to_use = filtered_meta if filtered_meta is not None else cached_meta
    
    if data_to_use is None:
        return dmc.Text("Start a search to see results", size="lg", c="dimmed")

    try:
        df = pd.DataFrame(data_to_use)       
        if df.empty:
            return dmc.Alert(
                "No results match your search criteria.",
                color="blue",
                variant="filled"
            )
        # Calculate the count for each accession_tag
        genome_counts = df.groupby('accession_tag').size().reset_index(name='n_genomes')
        
        # Remove duplicates and merge with counts
        filtered_meta_df = df.drop_duplicates(subset=['accession_tag']).merge(
            genome_counts, 
            on='accession_tag', 
            how='left'
        )
        
        if filtered_meta_df.empty:
            return dmc.Text(
                "No results found",
                size="lg",
                c="dimmed"
            )

        table_columns = [
            {
                "name": "Accession",
                "id": "accession_tag",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
            {
                "name": "Starship Family",
                "id": "familyName",
                "deletable": False,
                "selectable": False,
            },
            {
                "name": "Starship Navis",
                "id": "starship_navis",
                "deletable": False,
                "selectable": False,
            },
            {
                "name": "Starship Haplotype",
                "id": "starship_haplotype",
                "deletable": False,
                "selectable": False,
            },
            {
                "name": "Number of genomes",
                "id": "n_genomes",
                "deletable": False,
                "selectable": False,
            },
            {
                "name": "Species",
                "id": "species",
                "deletable": False,
                "selectable": False,
            },
            {
                "name": "Element Length (bp)",
                "id": "size",
                "deletable": False,
                "selectable": False,
            },
        ]
        
        table = make_ship_table(
            filtered_meta_df, 
            id="wiki-table", 
            columns=table_columns, 
            select_rows=False,
            pg_sz=15
        )
        
        title = dmc.Title("Search Results", order=2, mb="md")
        return dmc.Paper(children=[
            dbc.Stack([title, table], gap=3)
        ], p="xl", radius="md", withBorder=True)

    except Exception as e:
        logger.error(f"Error in create_search_results: {str(e)}", exc_info=True)
        return dmc.Alert(
            f"An error occurred while loading the results: {str(e)}",
            color="red",
            variant="filled"
        )


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
        logger.warning("No active item or meta data provided.")
        raise PreventUpdate
    
    try:
        title = dmc.Title(f"Taxonomy Distribution for {active_item}", order=2, mb="md")
        
        # Load or create sunburst plot
        sunburst_figure = load_from_cache(f"sunburst_{active_item}")
        if sunburst_figure is None:
            df = pd.DataFrame(cached_meta)
            filtered_df = df[df["familyName"] == active_item]
            sunburst_figure = cache_sunburst_plot(
                family=active_item, df=filtered_df
            )
        
        # Make the plot responsive
        sunburst_figure.update_layout(
            autosize=True,
            margin=dict(l=0, r=0, t=30, b=0),  # Reduce margins
            height=None,  # Allow height to be determined by container
        )

        fig = dcc.Graph(
            figure=sunburst_figure,
            style={
                "width": "100%",
                "height": "100%"
            },
            config={
                'responsive': True,
                'displayModeBar': False,  # Hide the mode bar for cleaner mobile view
                'scrollZoom': False  # Disable scroll zoom on mobile
            }
        )

        return fig, title, active_item
    except Exception as e:
        logger.error(
            f"Error occurred while creating the sidebar for {active_item}.",
            exc_info=True,
        )
        raise


toggle_modal = create_modal_callback(
    "wiki-table",
    "wiki-modal",
    "wiki-modal-content",
    "wiki-modal-title"
)

# Add this cache decorator for common filter combinations
@lru_cache(maxsize=128)
def get_filtered_options(taxonomy_tuple, family_tuple, navis_tuple, haplotype_tuple, data_hash):
    """Cache-friendly version of option filtering"""
    df = pd.DataFrame(pickle.loads(data_hash))
    
    if taxonomy_tuple:
        df = df[df["genus"].isin(taxonomy_tuple)]
    if family_tuple:
        df = df[df["familyName"].isin(family_tuple)]
    if navis_tuple:
        df = df[df["starship_navis"].isin(navis_tuple)]
    if haplotype_tuple:
        df = df[df["starship_haplotype"].isin(haplotype_tuple)]
    
    return {
        "taxonomy": sorted(df["genus"].dropna().unique()),
        "family": sorted(df["familyName"].dropna().unique()),
        "navis": sorted(df["starship_navis"].dropna().unique()),
        "haplotype": sorted(df["starship_haplotype"].dropna().unique())
    }

@callback(
    [
        Output("taxonomy-search", "data"),
        Output("family-search", "data"),
        Output("navis-search", "data"),
        Output("haplotype-search", "data"),
    ],
    [
        Input("taxonomy-search", "value"),
        Input("family-search", "value"),
        Input("navis-search", "value"),
        Input("haplotype-search", "value"),
        Input("meta-data", "data"),
    ]
)
def update_search_options(taxonomy_val, family_val, navis_val, haplotype_val, meta_data):
    if not meta_data:
        empty_data = []
        return empty_data, empty_data, empty_data, empty_data
    
    try:
        df = pd.DataFrame(meta_data)
        
        # Apply filters based on current selections
        if taxonomy_val:
            df = df[df["genus"].isin(taxonomy_val)]
        if family_val:
            df = df[df["familyName"].isin(family_val)]
        if navis_val:
            df = df[df["starship_navis"].isin(navis_val)]
        if haplotype_val:
            df = df[df["starship_haplotype"].isin(haplotype_val)]
        
        # Get available options based on filtered data
        taxonomy_options = sorted(df["genus"].dropna().unique())
        family_options = sorted(df["familyName"].dropna().unique())
        navis_options = sorted(df["starship_navis"].dropna().unique())
        haplotype_options = sorted(df["starship_haplotype"].dropna().unique())
        
        # Format options for Mantine MultiSelect
        taxonomy_data = [{"value": x, "label": x} for x in taxonomy_options]
        family_data = [{"value": x, "label": x} for x in family_options]
        navis_data = [{"value": x, "label": x} for x in navis_options]
        haplotype_data = [{"value": x, "label": x} for x in haplotype_options]
        
        return taxonomy_data, family_data, navis_data, haplotype_data
    except Exception as e:
        logger.error(f"Error in update_search_options: {str(e)}", exc_info=True)
        empty_data = []
        return empty_data, empty_data, empty_data, empty_data

@callback(
    [
        Output("filtered-meta-data", "data"),
        Output("taxonomy-search", "value"),
        Output("family-search", "value"),
        Output("navis-search", "value"),
        Output("haplotype-search", "value"),
    ],
    [
        Input("apply-search", "n_clicks"),
        Input("reset-search", "n_clicks"),
    ],
    [
        State("taxonomy-search", "value"),
        State("family-search", "value"),
        State("navis-search", "value"),
        State("haplotype-search", "value"),
        State("meta-data", "data"),
    ],
    prevent_initial_call=True
)
def handle_search(search_clicks, reset_clicks, taxonomy, family, navis, haplotype, original_data):
    if not original_data:
        raise PreventUpdate
        
    triggered_id = dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    
    # If reset button clicked, clear all filters
    if triggered_id == "reset-search":
        return None, [], [], [], []  # Return None for filtered data to use original
    
    # If no search clicked, prevent update
    if not search_clicks:
        raise PreventUpdate
    
    df = pd.DataFrame(original_data)
    
    # Apply filters if they exist
    if taxonomy and len(taxonomy) > 0:
        df = df[df["genus"].isin(taxonomy)]
    if family and len(family) > 0:
        df = df[df["familyName"].isin(family)]
    if navis and len(navis) > 0:
        df = df[df["starship_navis"].isin(navis)]
    if haplotype and len(haplotype) > 0:
        df = df[df["starship_haplotype"].isin(haplotype)]
    
    # If no filters selected, return None to use original data
    if not any([taxonomy, family, navis, haplotype]):
        return None, [], [], [], []
    
    # Return filtered data and keep current filter values
    return df.to_dict("records"), taxonomy or [], family or [], navis or [], haplotype or []
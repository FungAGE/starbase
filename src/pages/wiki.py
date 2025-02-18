import warnings
import json

warnings.filterwarnings("ignore")
import dash
from dash import dcc, html, callback, clientside_callback, ClientsideFunction
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import pandas as pd
import os
import logging

from sqlalchemy.exc import SQLAlchemyError

from src.config.cache import cache
from src.database.sql_manager import fetch_meta_data
from src.database.sql_manager import (
    fetch_meta_data,
    fetch_paper_data,
)
from src.components.tables import make_ship_table, make_wiki_table
from src.utils.plot_utils import make_logo, create_sunburst_plot
from src.utils.seq_utils import clean_contigIDs
from src.components.callbacks import create_modal_callback

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
    try:
        meta_data = cache.get("meta_data")
        if meta_data is None:
            logger.debug("Cache miss for meta_data, fetching from database")
            meta_data = fetch_meta_data()
            if meta_data is not None:
                cache.set("meta_data", meta_data)
        
        if isinstance(meta_data, pd.DataFrame):
            return meta_data.to_dict('records')
        return meta_data
    except Exception as e:
        logger.error(f"Error loading initial data: {str(e)}")
        return None

layout = dmc.Container(
    fluid=True,
    children=[
        modal,
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="meta-data", data=load_initial_data()),
        dcc.Store(id="filtered-meta-data"),
        dcc.Store(id="paper-data"),
        
        # Header Section
        dmc.Space(h="md"),        
        dmc.Paper(
            children=[
                dmc.Title(
                    "Starship Wiki",
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    "Search and explore the characteristics of different Starship families",
                    c="dimmed",
                    size="lg",
                ),
                dmc.Space(h="md"),
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
        dmc.Grid(
            children=[
                # Left Column - Search Results
                dmc.GridCol(
                    span={"lg": 6, "md": 12},
                    children=[
                        html.Div(id="search-results"),
                    ]
                ),
                # Right Column - Search Section
                dmc.GridCol(
                    span={"lg": 6, "md": 12},
                    children=[
                        dmc.Paper(
                            children=[
                                dmc.Title("Taxonomic Distribution", order=2, mb="md"),
                                dcc.Loading(
                                    id="search-sunburst-loading",
                                    type="circle",
                                    children=html.Div(
                                        id="search-sunburst-plot",
                                        style={
                                            "width": "100%",
                                            "height": "100%",
                                            "overflow": "hidden",  # Prevent overflow
                                            "display": "flex",     # Use flexbox
                                            "alignItems": "center",
                                            "justifyContent": "center",
                                        },
                                    ),
                                ),
                            ],
                            withBorder=True,
                            shadow="sm",
                            radius="md",
                            p="md"
                        ),
                        dmc.Space(h="sm"),
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
                            style={
                                "minHeight": "200px",  # Minimum height when collapsed
                                "maxHeight": "calc(100vh - 200px)",  # Maximum height
                                "height": "auto",  # Allow height to adjust to content
                                "overflowY": "auto"
                            },
                        ),

                ]),
            ],
        ),
    ],
)


@cache.memoize()
@callback(Output("meta-data", "data"), Input("url", "href"))
def load_meta_data(url):
    if not url:
        raise PreventUpdate
    
    try:
        meta_data = cache.get("meta_data")
        if meta_data is None:
            meta_data = fetch_meta_data()
            if meta_data is not None:
                cache.set("meta_data", meta_data)
                
        if meta_data is None:
            logger.error("Failed to fetch metadata")
            return []
            
        if "contigID" in meta_data.columns:
            meta_data["contigID"] = meta_data["contigID"].apply(clean_contigIDs)
            
        return meta_data.to_dict("records")
    except Exception as e:
        logger.error(f"Error loading meta data: {str(e)}")
        return []


# Callback to load paper data
@cache.memoize()
@callback(Output("paper-data", "data"), Input("url", "href"))
def load_paper_data(url):
    if url:
        try:
            logger.debug("Loading paper data from database.")
            paper_df = cache.get("paper_data")
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
            active_item=[],
        )
    except Exception as e:
        logger.error("Error occurred while creating the accordion.", exc_info=True)
        raise


@callback(
    Output("search-results", "children"),
    Input("filtered-meta-data", "data"),
    Input("meta-data", "data"),  # Change State to Input to handle initial load
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

        # Updated column definitions for AG Grid
        table_columns = [
            {
                "field": "accession_tag",
                "headerName": "Accession",
                "flex": 1,
                "cellStyle": {"cursor": "pointer", "color": "#1976d2"}
            },
            {
                "field": "familyName",
                "headerName": "Starship Family",
                "flex": 1
            },
            {
                "field": "n_genomes",
                "headerName": "Number of Genomes",
                "flex": 1,
                "type": "numericColumn",
                "valueFormatter": "value.toLocaleString()",
                "sortable": True
            },
            {
                "field": "species",
                "headerName": "Species",
                "flex": 1
            },
            {
                "field": "size",
                "headerName": "Element Length (bp)",
                "flex": 1,
                "type": "numericColumn",
                "valueFormatter": "value.toLocaleString()",
                "sortable": True
            }
        ]
        
        # Convert numeric columns to appropriate type
        filtered_meta_df['n_genomes'] = pd.to_numeric(filtered_meta_df['n_genomes'], errors='coerce')
        filtered_meta_df['size'] = pd.to_numeric(filtered_meta_df['size'], errors='coerce')
        
        # Fill NA values
        filtered_meta_df = filtered_meta_df.fillna('')
        
        table = make_ship_table(
            filtered_meta_df, 
            id="wiki-table", 
            columns=table_columns, 
            select_rows=False,
            pg_sz=25
        )
        
        title = dmc.Title("Search Results", order=2, mb="md")
        return dmc.Paper(
            children=[
                dbc.Stack([title, table], gap=3)
            ],
            p="xl",
            radius="md",
            withBorder=True,
            style={
                "minHeight": "calc(100vh - 120px)",
                "height": "auto",
                "maxHeight": "none",
                "overflowY": "visible",
                "marginBottom": "2rem"
            }
        )

    except Exception as e:
        logger.error(f"Error in create_search_results: {str(e)}", exc_info=True)
        return dmc.Alert(
            f"An error occurred while loading the results: {str(e)}",
            color="red",
            variant="filled"
        )


toggle_modal = create_modal_callback(
    "wiki-table",
    "wiki-modal",
    "wiki-modal-content",
    "wiki-modal-title"
)

# Add this cache decorator for common filter combinations
@cache.memoize()
def get_filtered_options(taxonomy=None, family=None, navis=None, haplotype=None):
    """Cache-friendly version of option filtering"""
    try:
        meta_data = cache.get("meta_data")
        if meta_data is None:
            meta_data = fetch_meta_data()
            
        df = pd.DataFrame(meta_data)
        
        if taxonomy:
            df = df[df["genus"].isin(taxonomy)]
        if family:
            df = df[df["familyName"].isin(family)]
        if navis:
            df = df[df["starship_navis"].isin(navis)]
        if haplotype:
            df = df[df["starship_haplotype"].isin(haplotype)]
        
        return {
            "taxonomy": sorted(df["genus"].dropna().unique()),
            "family": sorted(df["familyName"].dropna().unique()),
            "navis": sorted(df["starship_navis"].dropna().unique()),
            "haplotype": sorted(df["starship_haplotype"].dropna().unique())
        }
    except Exception as e:
        logger.error(f"Error in get_filtered_options: {str(e)}")
        return {
            "taxonomy": [], "family": [], "navis": [], "haplotype": []
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
    ],
    prevent_initial_call=True
)
def update_search_options(taxonomy_val, family_val, navis_val, haplotype_val, meta_data):
    if not meta_data:
        empty_data = []
        return empty_data, empty_data, empty_data, empty_data
    
    try:
        # Get filtered options from cache if possible
        options = get_filtered_options(taxonomy_val, family_val, navis_val, haplotype_val)
        
        # Format options for Mantine MultiSelect
        taxonomy_data = [{"value": x, "label": x} for x in options["taxonomy"]]
        family_data = [{"value": x, "label": x} for x in options["family"]]
        navis_data = [{"value": x, "label": x} for x in options["navis"]]
        haplotype_data = [{"value": x, "label": x} for x in options["haplotype"]]
        
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
        
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate
        
    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]
    
    # If reset button clicked, clear all filters
    if triggered_id == "reset-search":
        return None, [], [], [], []
    
    # If no search clicked or no filters selected, prevent update
    if not search_clicks or not any([taxonomy, family, navis, haplotype]):
        return None, [], [], [], []
    
    try:
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
        
        # Return filtered data and keep current filter values
        return df.to_dict("records"), taxonomy or [], family or [], navis or [], haplotype or []
        
    except Exception as e:
        logger.error(f"Error in handle_search: {str(e)}", exc_info=True)
        return None, [], [], [], []

@callback(
    Output("search-sunburst-plot", "children"),
    [Input("filtered-meta-data", "data"),
     Input("meta-data", "data")],
)
def update_search_sunburst(filtered_meta, meta_data):
    # For initial load, use meta_data
    if filtered_meta is None and meta_data is not None:
        data_to_use = meta_data
    # For filtered results, use filtered_meta
    elif filtered_meta is not None:
        data_to_use = filtered_meta
    else:
        return dmc.Text("No data available", size="lg", c="dimmed")

    try:
        df = pd.DataFrame(data_to_use)
        if df.empty:
            return dmc.Text("No results to display", size="lg", c="dimmed")

        # Create sunburst plot
        sunburst_figure = create_sunburst_plot(
            df=df,
            type="tax",
            title_switch=False
        )
        
        if sunburst_figure is None:
            return dmc.Text("No data available", size="lg", c="dimmed")
        
        # Enhanced layout settings for consistent sizing
        sunburst_figure.update_layout(
            autosize=True,
            margin=dict(l=0, r=0, t=0, b=0, pad=0),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            uniformtext=dict(
                minsize=12,
                mode='hide'
            ),
            uirevision='constant',
            # Remove fixed height to allow container to control size
            width=None,
            height=None,
        )

        # Update trace settings
        sunburst_figure.update_traces(
            textinfo='label+value',
            insidetextorientation='radial',
            hoverinfo='label+value+percent parent',
            hovertemplate=(
                "<b>%{label}</b><br>" +
                "Count: %{value}<br>" +
                "Percentage: %{percentParent:.1%}<br>" +
                "<extra></extra>"
            ),
            textfont=dict(size=12),
        )

        return dcc.Graph(
                    figure=sunburst_figure,
                    style={
                        "width": "100%",
                        "height": "600px",  # Fixed height in container
                    },
                    config={
                        'responsive': True,
                        'displayModeBar': False,
                        'scrollZoom': False,
                        'staticPlot': False,
                        'doubleClick': False,
                        'showTips': False,
                        'showAxisDragHandles': False,
                        'displaylogo': False,
                        'showLink': False,
                        'animate': False,
                        'editable': False,
                        'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
                        'toImageButtonOptions': {'height': None, 'width': None},
                    },
                    clear_on_unhover=True,
                    className='plot-container'
                )
    except Exception as e:
        logger.error(f"Error in update_search_sunburst: {str(e)}")
        return dmc.Alert(
            f"An error occurred while creating the plot: {str(e)}",
            color="red",
            variant="filled"
        )
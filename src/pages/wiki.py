from datetime import datetime

import dash
from dash import dcc, html, callback, clientside_callback
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

# import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify

import pandas as pd
import os

from sqlalchemy.exc import SQLAlchemyError

from src.config.cache import cache
from src.database.sql_manager import (
    fetch_meta_data,
    fetch_paper_data,
    fetch_ships,
    dereplicate_sequences,
)
from src.components.tables import (
    make_dl_table,
    make_wiki_table,
)
from src.utils.plot_utils import make_logo, create_sunburst_plot
from src.utils.seq_utils import clean_contigIDs, create_ncbi_style_header
from src.components.callbacks import (
    curated_switch,
    dereplicated_switch,
)
from src.components.callbacks import handle_callback_error
from src.config.logging import get_logger

logger = get_logger(__name__)

dash.register_page(__name__)


def create_accordion_item(df, papers, category):
    if category == "nan":
        return None
    else:
        filtered_meta_df = df[df["familyName"] == category]
        n_ships = len(filtered_meta_df)

        # Element length consolidation: prioritize elementLength, fallback to calculated or sequence_length
        # All length fields are now properly populated in the database
        element_lengths = filtered_meta_df["elementLength"].dropna()
        if len(element_lengths) == 0:
            # Calculate lengths from begin/end positions
            begin_vals = filtered_meta_df["elementBegin"].dropna()
            end_vals = filtered_meta_df["elementEnd"].dropna()
            if len(begin_vals) > 0 and len(end_vals) > 0:
                # Use the minimum length of begin/end pairs to avoid mismatched lengths
                min_len = min(len(begin_vals), len(end_vals))
                element_lengths = abs(
                    begin_vals.iloc[:min_len] - end_vals.iloc[:min_len]
                )
            else:
                element_lengths = filtered_meta_df["sequence_length"].dropna()

        if len(element_lengths) == 0:
            min_size = 0
            max_size = 0
        else:
            min_size = int(element_lengths.min())
            max_size = int(element_lengths.max())
        upDRs = filtered_meta_df["upDR"].dropna().tolist()
        downDRs = filtered_meta_df["downDR"].dropna().tolist()
        filtered_papers_df = papers[papers["familyName"] == category]
        type_element_reference = (
            filtered_papers_df["type_element_reference"].dropna().unique().astype(str)
        )

        if len(upDRs) > 10:
            uplogo_img_path = f"assets/images/DR/{category}-upDR.png"

            if not os.path.exists(uplogo_img_path):
                created_path = make_logo(upDRs, uplogo_img_path, type="up")
                if created_path:
                    uplogo_img_path = created_path

            if os.path.exists(uplogo_img_path):
                uplogo_img = dbc.Col(
                    lg=6,
                    sm=12,
                    children=[
                        dmc.Center(html.H5("Upstream DRs")),
                        dmc.Center(
                            html.Img(
                                src=f"/{uplogo_img_path}",
                                style={"width": "100%"},
                            ),
                        ),
                    ],
                )
            else:
                uplogo_img = None
        else:
            uplogo_img = None

        if len(downDRs) > 10:
            downlogo_img_path = f"assets/images/DR/{category}-downDR.png"
            if not os.path.exists(downlogo_img_path):
                created_path = make_logo(downDRs, downlogo_img_path, type="down")
                if created_path:
                    downlogo_img_path = created_path

            if os.path.exists(downlogo_img_path):
                downlogo_img = dbc.Col(
                    lg=6,
                    sm=12,
                    children=[
                        dmc.Center(html.H5("Downstream DRs")),
                        dmc.Center(
                            html.Img(
                                src=f"/{downlogo_img_path}",
                                style={"width": "100%"},
                            ),
                        ),
                    ],
                )
            else:
                downlogo_img = None
        else:
            downlogo_img = None

        accordion_content = [make_wiki_table(n_ships, max_size, min_size)]
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
        meta_data = fetch_meta_data()
        if meta_data is not None and isinstance(meta_data, pd.DataFrame):
            return meta_data.to_dict("records")
        return meta_data
    except Exception as e:
        logger.error(f"Error loading initial data: {str(e)}")
        return None


table_columns = [
    {
        "name": "Ship Accession (SSB)",
        "id": "ship_accession_tag",
        "deletable": False,
        "selectable": True,
        "presentation": "markdown",
        "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
    },
    {
        "name": "Group Accession (SSA)",
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
        "presentation": "markdown",
    },
    {
        "name": "Order",
        "id": "order",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Family",
        "id": "family",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Species",
        "id": "name",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
    {
        "name": "Reference",
        "id": "shortCitation",
        "deletable": False,
        "selectable": False,
        "presentation": "markdown",
    },
]


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
                    [
                        html.Span(
                            "starbase ",
                            className="logo-text",
                        ),
                        " Wiki",
                    ],
                    order=1,
                    mb="md",
                ),
                dmc.Text(
                    [
                        "Search and explore the characteristics of different ",
                        html.Span("Starship", style={"fontStyle": "italic"}),
                        " families",
                    ],
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
                                dmc.Autocomplete(
                                    id="taxa-search",
                                    label="Taxonomy",
                                    placeholder="Search across family, genus, etc.",
                                    data=[],  # Will be populated by callback
                                    limit=20,  # Limit dropdown items for performance
                                    style={"width": "100%"},
                                )
                            ],
                        ),
                        # Family Search
                        dmc.GridCol(
                            span={"lg": 3, "md": 6, "sm": 12},
                            children=[
                                dmc.Autocomplete(
                                    id="family-search",
                                    label="Starship Family",
                                    placeholder="Search by family...",
                                    data=[],  # Will be populated by callback
                                    limit=20,  # Limit dropdown items for performance
                                    style={"width": "100%"},
                                ),
                            ],
                        ),
                    ],
                    gutter="xl",
                ),
                dmc.Group(
                    align="right",
                    mt="md",
                    children=[
                        dmc.Group(
                            children=[
                                curated_switch(
                                    text=html.Div(
                                        [
                                            "Only show curated ",
                                            html.Span(
                                                "Starships",
                                                style={"fontStyle": "italic"},
                                            ),
                                        ]
                                    ),
                                    size="md",
                                ),
                                dereplicated_switch(
                                    text=html.Div(
                                        [
                                            "Only show dereplicated ",
                                            html.Span(
                                                "Starships",
                                                style={"fontStyle": "italic"},
                                            ),
                                        ]
                                    ),
                                    size="md",
                                ),
                            ],
                            mt="md",
                        ),
                    ],
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
                    ],
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
                                            "display": "flex",  # Use flexbox
                                            "alignItems": "center",
                                            "justifyContent": "center",
                                        },
                                    ),
                                ),
                            ],
                            p="xl",
                            radius="md",
                            withBorder=True,
                            mb="xl",
                        ),
                        dmc.Space(h="sm"),
                        dmc.Paper(
                            children=[
                                dmc.Title(
                                    html.Div(
                                        [
                                            html.Span(
                                                "Starship",
                                                style={"fontStyle": "italic"},
                                            ),
                                            " Families",
                                        ]
                                    ),
                                    order=2,
                                    mb="md",
                                ),
                                dcc.Loading(
                                    id="wiki-loading",
                                    type="circle",
                                    children=html.Div(
                                        id="accordion",
                                        children=dbc.Accordion(id="category-accordion"),
                                    ),
                                ),
                            ],
                            p="xl",
                            radius="md",
                            withBorder=True,
                            mb="xl",
                            style={
                                "minHeight": "200px",  # Minimum height when collapsed
                                "maxHeight": "calc(100vh - 200px)",  # Maximum height
                                "height": "auto",  # Allow height to adjust to content
                                "overflowY": "auto",
                            },
                        ),
                    ],
                ),
            ],
        ),
    ],
)


@callback(Output("meta-data", "data"), Input("url", "href"))
@handle_callback_error
def load_meta_data(url):
    if not url:
        raise PreventUpdate

    try:
        meta_data = fetch_meta_data()

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
@handle_callback_error
@callback(Output("paper-data", "data"), Input("url", "href"))
def load_paper_data(url):
    if url:
        try:
            logger.debug("Loading paper data from database.")
            paper_df = cache.get("paper_data")
            if paper_df is None:
                paper_df = fetch_paper_data()
            logger.debug(f"Paper data query returned {len(paper_df)} rows.")
            return paper_df.to_dict("records")
        except SQLAlchemyError as e:
            logger.error(
                f"Error occurred while executing the paper data query: {str(e)}"
            )
            return []

    logger.warning("No URL provided for paper data loading.")
    raise PreventUpdate


# Callback to create accordion layout
@callback(
    Output("accordion", "children"),
    [Input("meta-data", "data"), Input("paper-data", "data")],
)
@handle_callback_error
def create_accordion(cached_meta, cached_papers):
    if cached_meta is None or cached_papers is None:
        logger.error("One or more inputs are None: meta-data or paper-data.")
        raise PreventUpdate

    df = pd.DataFrame(cached_meta)
    papers = pd.DataFrame(cached_papers)

    assert isinstance(df, pd.DataFrame), (
        f"Expected df to be a DataFrame, but got {type(df)}."
    )
    assert isinstance(papers, pd.DataFrame), (
        f"Expected papers to be a DataFrame, but got {type(papers)}."
    )

    logger.debug("Creating accordion based on metadata and paper data.")

    try:
        unique_categories = df["familyName"].dropna().unique().tolist()
        logger.debug(
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
        logger.error(f"Error occurred while creating the accordion: {str(e)}")
        raise


@callback(
    Output("search-results", "children"),
    [
        Input("filtered-meta-data", "data"),
        Input("meta-data", "data"),
        Input("curated-input", "checked"),
        Input("dereplicated-input", "checked"),
    ],
)
@handle_callback_error
def create_search_results(filtered_meta, cached_meta, curated, dereplicate):
    # Use filtered data if available, otherwise use original data
    data_to_use = filtered_meta if filtered_meta is not None else cached_meta

    if data_to_use is None:
        return dmc.Text("Start a search to see results", size="lg", c="dimmed")

    try:
        df = pd.DataFrame(data_to_use)

        # Check if DataFrame is empty before applying filters
        if df.empty:
            return dmc.Alert(
                "No results match your search criteria.", color="blue", variant="filled"
            )

        # Apply curated/dereplicated filters if switches are enabled
        if curated and "curated_status" in df.columns:
            df = df[df["curated_status"] == "curated"]

        # Check again after applying filters
        if df.empty:
            return dmc.Alert(
                "No results match your search criteria.", color="blue", variant="filled"
            )

        if dereplicate:
            filtered_meta_df = dereplicate_sequences(df)
        else:
            filtered_meta_df = df.copy()

        if filtered_meta_df.empty:
            return dmc.Text("No results found", size="lg", c="dimmed")

        # drop unnecessary columns
        filtered_meta_df = filtered_meta_df.drop(columns=["elementLength"])

        # Fill NA values
        filtered_meta_df = filtered_meta_df.fillna("")

        table = make_dl_table(
            filtered_meta_df,
            id="dl-table",
            table_columns=table_columns,
        )

        # Create the enhanced component similar to wiki_table_with_download
        enhanced_results = dmc.Paper(
            children=[
                # Main Content
                dmc.Title("Info Table and Downloads", order=2, mb="md"),
                dmc.Stack(
                    [
                        dmc.Text(
                            [
                                "Select individual ",
                                html.Span("Starships", style={"fontStyle": "italic"}),
                                " or download the complete dataset",
                            ],
                            size="lg",
                            c="dimmed",
                        ),
                        # Download Options
                        dmc.Center(
                            dmc.Group(
                                gap="xl",
                                children=[
                                    dmc.Button(
                                        html.Div(
                                            [
                                                "Download All ",
                                                html.Span(
                                                    "Starships",
                                                    style={"fontStyle": "italic"},
                                                ),
                                            ]
                                        ),
                                        id="download-all-btn",
                                        variant="gradient",
                                        gradient={
                                            "from": "indigo",
                                            "to": "cyan",
                                        },
                                        leftSection=html.I(
                                            className="bi bi-cloud-download"
                                        ),
                                        size="md",
                                        loaderProps={
                                            "variant": "dots",
                                            "color": "white",
                                        },
                                    ),
                                    dmc.Button(
                                        html.Div(
                                            [
                                                "Download Selected ",
                                                html.Span(
                                                    "Starships",
                                                    style={"fontStyle": "italic"},
                                                ),
                                            ]
                                        ),
                                        id="download-selected-btn",
                                        variant="gradient",
                                        gradient={"from": "teal", "to": "lime"},
                                        leftSection=html.I(className="bi bi-download"),
                                        size="md",
                                        loaderProps={
                                            "variant": "dots",
                                            "color": "white",
                                        },
                                    ),
                                ],
                            ),
                        ),
                    ],
                    gap="xl",
                ),
                dmc.Space(h="md"),
                # Notification Area
                html.Div(id="dl-notify"),
                dcc.Download(id="dl-package"),
                # Table Section
                html.Div(
                    children=[
                        # Table Header with Stats
                        dmc.Group(
                            gap="apart",
                            mt="md",
                            children=[
                                dmc.Text(
                                    id="table-stats",
                                    size="sm",
                                    c="dimmed",
                                ),
                                dmc.Text(
                                    [
                                        "Click rows to select ",
                                        html.Span(
                                            "Starships", style={"fontStyle": "italic"}
                                        ),
                                    ],
                                    size="sm",
                                    c="dimmed",
                                    style={"fontStyle": "italic"},
                                ),
                            ],
                        ),
                        dmc.Modal(
                            id="accession-modal",
                            opened=False,
                            centered=True,
                            overlayProps={"blur": 3},
                            size="lg",
                            children=[
                                dmc.Title(id="modal-title", order=3),
                                dmc.Space(h="md"),
                                html.Div(id="modal-content"),
                            ],
                        ),
                        html.Div(
                            id="dl-table-container",
                            children=[table],
                        ),
                        # Dummy div for clientside callback
                        html.Div(id="dummy-output", style={"display": "none"}),
                    ],
                ),
            ],
            p="xl",
            radius="md",
            withBorder=True,
            mb="xl",
            style={
                "minHeight": "calc(100vh - 120px)",
                "height": "auto",
                "maxHeight": "none",
                "overflowY": "visible",
                "marginBottom": "2rem",
            },
        )

        return enhanced_results

    except Exception as e:
        logger.error(f"Error in create_search_results: {str(e)}")
        return dmc.Alert(
            f"An error occurred while loading the results: {str(e)}",
            color="red",
            variant="filled",
        )


# Callback to populate autocomplete components with initial data
@callback(
    [
        Output("taxa-search", "data"),
        Output("family-search", "data"),
    ],
    Input("meta-data", "data"),
    prevent_initial_call=False,  # Allow initial call to populate on page load
)
@handle_callback_error
def populate_search_components(meta_data):
    if not meta_data:
        return [], []

    try:
        df = pd.DataFrame(meta_data)

        # Get taxonomy search data
        search_columns = [
            "name",
            "subkingdom",
            "phylum",
            "subphylum",
            "class",
            "subclass",
            "order",
            "suborder",
            "family",
            "genus",
        ]

        # Collect all unique values across specified columns
        all_taxa_values = set()

        for col in search_columns:
            if col in df.columns:
                # Get non-null values and add to set
                values = df[col].dropna().astype(str).unique()
                all_taxa_values.update(values)

        # Convert to sorted list and format for Autocomplete
        # sort by length instead of alphabetically
        taxa_search_data = [
            {"value": val, "label": val}
            for val in sorted(all_taxa_values, key=lambda s: len(s))
        ]

        # Get family search data
        family_values = []
        if "familyName" in df.columns:
            family_values = df["familyName"].dropna().astype(str).unique()
            family_search_data = [
                {"value": val, "label": val}
                for val in sorted(family_values, key=lambda s: len(s))
            ]
        else:
            family_search_data = []

        return taxa_search_data, family_search_data

    except Exception as e:
        logger.error(f"Error in populate_search_components: {str(e)}")
        return [], []


@callback(
    Output("filtered-meta-data", "data"),
    [
        Input("apply-search", "n_clicks"),
        Input("reset-search", "n_clicks"),
    ],
    [
        State("taxa-search", "value"),
        State("family-search", "value"),
        State("meta-data", "data"),
    ],
    prevent_initial_call=True,
)
@handle_callback_error
def handle_taxa_and_family_search(
    search_clicks, reset_clicks, taxa_search_value, family_search_value, original_data
):
    if not original_data:
        raise PreventUpdate

    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    triggered_id = ctx.triggered[0]["prop_id"].split(".")[0]

    # Reset functionality
    if triggered_id == "reset-search":
        return None

    # If no search clicked, prevent update
    if not search_clicks:
        return None

    try:
        df = pd.DataFrame(original_data)
        filtered_df = df.copy()

        # Apply taxonomy search if value is provided
        if taxa_search_value and taxa_search_value.strip():
            # Define columns to search across
            taxa_search_columns = [
                "name",
                "subkingdom",
                "phylum",
                "subphylum",
                "class",
                "subclass",
                "order",
                "suborder",
                "family",
                "genus",
            ]

            # Create a mask for rows that contain the search value in any of the specified columns
            mask = pd.Series([False] * len(filtered_df))

            for col in taxa_search_columns:
                if col in filtered_df.columns:
                    # Case-insensitive partial matching
                    mask |= (
                        filtered_df[col]
                        .astype(str)
                        .str.contains(taxa_search_value, case=False, na=False)
                    )

            filtered_df = filtered_df[mask]

        # Apply family search if value is provided
        if family_search_value and family_search_value.strip():
            if "familyName" in filtered_df.columns:
                # Case-insensitive exact matching for family
                filtered_df = filtered_df[
                    filtered_df["familyName"].astype(str).str.lower()
                    == family_search_value.lower()
                ]

        # Return empty list if no results found, otherwise return the filtered data
        if filtered_df.empty:
            return []
        else:
            return filtered_df.to_dict("records")

    except Exception as e:
        logger.error(f"Error in handle_taxa_and_family_search: {str(e)}")
        return []


# TODO: only change sunburst plot if the "Search" button is clicked? instead of if the switches are changed?
@callback(
    Output("search-sunburst-plot", "children"),
    [
        Input("filtered-meta-data", "data"),
        Input("meta-data", "data"),
        Input("curated-input", "checked"),
        Input("dereplicated-input", "checked"),
    ],
    prevent_initial_call=False,  # Allow initial call
)
@handle_callback_error
def update_search_sunburst(filtered_meta, meta_data, curated, dereplicate):
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

        if curated:
            df = df[df["curated_status"] == "curated"]

        # Deduplicate data to match table processing
        if dereplicate:
            df = dereplicate_sequences(df)

        # Create sunburst plot
        sunburst_figure = create_sunburst_plot(
            df=df,
            type="tax",
            title_switch=False,
        )

        if sunburst_figure is None:
            return dmc.Text("No data available", size="lg", c="dimmed")

        # Enhanced layout settings for consistent sizing
        sunburst_figure.update_layout(
            autosize=True,
            margin=dict(l=0, r=0, t=0, b=0, pad=0),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            uniformtext=dict(minsize=12, mode="hide"),
            uirevision="constant",
            # Remove fixed height to allow container to control size
            width=None,
            height=None,
        )

        # Update trace settings
        sunburst_figure.update_traces(
            textinfo="label+value",
            insidetextorientation="radial",
            hoverinfo="label+value+percent parent",
            hovertemplate=(
                "<b>%{label}</b><br>"
                + "Count: %{value}<br>"
                + "Percentage: %{percentParent:.1%}<br>"
                + "<extra></extra>"
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
                "responsive": True,
                "displayModeBar": False,
                "scrollZoom": False,
                "staticPlot": False,
                "doubleClick": False,
                "showTips": False,
                "showAxisDragHandles": False,
                "displaylogo": False,
                "showLink": False,
                "editable": False,
                "modeBarButtonsToRemove": ["select2d", "lasso2d"],
                "toImageButtonOptions": {"height": None, "width": None},
            },
            clear_on_unhover=True,
            className="plot-container",
        )
    except Exception as e:
        logger.error(f"Error in update_search_sunburst: {str(e)}")
        return dmc.Alert(
            f"An error occurred while creating the plot: {str(e)}",
            color="red",
            variant="filled",
        )


@callback(
    Output("table-stats", "children"),
    Input("filtered-meta-data", "data"),
    Input("meta-data", "data"),
    Input("curated-input", "checked"),
    Input("dereplicated-input", "checked"),
)
@handle_callback_error
def update_table_stats(filtered_meta, cached_meta, curated, dereplicate):
    """Update the table statistics based on current data and filters"""
    data_to_use = filtered_meta if filtered_meta is not None else cached_meta

    if data_to_use is None:
        return "No data available"

    try:
        df = pd.DataFrame(data_to_use)

        # Apply curated/dereplicated filters if switches are enabled
        if curated:
            df = df[df["curated_status"] == "curated"]

        if dereplicate:
            df = dereplicate_sequences(df)

        if df.empty:
            return "No records found"

        # Count unique accession ships
        unique_ships = len(df)
        return f"Showing {unique_ships} records"

    except Exception as e:
        logger.error(f"Error updating table stats: {str(e)}")
        return "Error loading data"


def generate_download_helper(rows, curated, dereplicate):
    """Helper function containing the common download logic"""
    import re

    try:
        if not rows:
            raise ValueError("No rows selected for download")

        accessions = [
            re.sub(pattern=r"\..*", repl="", string=row["accession_tag"])
            for row in rows
        ]
        dl_df = fetch_ships(
            accession_tags=accessions,
            curated=curated,
            dereplicate=dereplicate,
            with_sequence=True,
        )

        if dl_df is None or dl_df.empty:
            raise ValueError("No sequences found for download")

        # Count occurrences of each accession tag
        accession_counts = dl_df["accession_tag"].value_counts()

        fasta_content = []
        for _, row in dl_df.drop_duplicates(
            subset=["accession_tag", "sequence"]
        ).iterrows():
            count = accession_counts[row["accession_tag"]]
            header = create_ncbi_style_header(row, count)
            # Skip if header creation failed (returns None)
            if header is None:
                logger.warning(
                    f"Skipping sequence with accession_tag={row.get('accession_tag', 'unknown')} due to header creation failure"
                )
                continue
            fasta_content.append(f"{header}\n{row['sequence']}")

        fasta_str = "\n".join(fasta_content)
        logger.debug(
            f"FASTA content created successfully for {len(fasta_content)} sequences."
        )

        # Return both the FASTA content and download info
        return dict(
            content=fasta_str,
            filename=f"starships_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
            type="text/plain",
        ), len(fasta_content)

    except ValueError as e:
        logger.warning(f"Download helper warning: {str(e)}")
        return None, 0
    except Exception as e:
        logger.error(f"Failed to execute database query. Details: {e}")
        return None, 0


@callback(
    [
        Output("dl-package", "data", allow_duplicate=True),
        Output("dl-notify", "children", allow_duplicate=True),
    ],
    [Input("download-all-btn", "n_clicks")],
    [
        State("dl-table", "rowData"),
        State("curated-input", "checked"),
        State("dereplicated-input", "checked"),
    ],
    prevent_initial_call=True,
    running=[
        (Output("download-all-btn", "loading"), True, False),
    ],
)
@handle_callback_error
def generate_download_all(dl_all_clicks, table_data, curated, dereplicate):
    if not dl_all_clicks:
        raise dash.exceptions.PreventUpdate

    download_data, num_sequences = generate_download_helper(
        table_data, curated, dereplicate
    )

    if not download_data:
        return None, dmc.Alert(
            children="No sequences found for download", color="red", title="Error"
        )

    if num_sequences > 0:
        if num_sequences == 1:
            download_message = f"Downloading {num_sequences} sequence"
        else:
            download_message = f"Downloading {num_sequences} sequences"
    else:
        download_message = "No sequences found for download"

    return download_data, dmc.Alert(
        children=download_message,
        color="green",
        title="Success",
    )


@callback(
    [Output("dl-package", "data"), Output("dl-notify", "children")],
    [Input("download-selected-btn", "n_clicks")],
    [
        State("dl-table", "rowData"),
        State("dl-table", "selectedRows"),
        State("curated-input", "checked"),
        State("dereplicated-input", "checked"),
    ],
    prevent_initial_call=True,
    running=[
        (Output("download-selected-btn", "loading"), True, False),
    ],
)
@handle_callback_error
def generate_download_selected(
    dl_select_clicks, table_data, selected_rows, curated, dereplicate
):
    if not dl_select_clicks or not selected_rows:
        raise dash.exceptions.PreventUpdate

    # Filter table_data to get only the selected rows
    selected_data = [
        row
        for row in table_data
        if row.get("accession_tag")
        in [selected.get("accession_tag") for selected in selected_rows]
    ]

    download_data, num_sequences = generate_download_helper(
        selected_data, curated, dereplicate
    )

    if not download_data:
        return None, dmc.Alert(
            children="No sequences found for selected rows", color="red", title="Error"
        )

    return download_data, dmc.Alert(
        children=f"Downloading {num_sequences} sequences",
        color="green",
        title="Success",
    )


@callback(
    Output("download-selected-btn", "disabled"),
    [Input("dl-table", "selectedRows")],  # Only selectedRows, no rowData
)
@handle_callback_error
def update_download_selected_button(selected_rows):
    # Disable the button if no rows are selected
    return not selected_rows or len(selected_rows) == 0


# Add clientside callback to handle accession modal clicks
clientside_callback(
    """
    function(cellClicked, activeCell, tableData, pageCurrent, pageSize) {
        if (!cellClicked && !activeCell) {
            return window.dash_clientside.no_update;
        }

        let accession = null;

        // Handle AG Grid cell clicks
        if (cellClicked && (cellClicked.colId === 'ship_accession_tag' || cellClicked.colId === 'ship_accession_display')) {
            accession = cellClicked.value;
        }
        // Handle DataTable active cell
        else if (activeCell && (activeCell.column_id === 'ship_accession_tag' || activeCell.column_id === 'ship_accession_display')) {
            const actualRowIdx = (pageCurrent || 0) * pageSize + activeCell.row;
            if (tableData && actualRowIdx < tableData.length) {
                accession = tableData[actualRowIdx][activeCell.column_id];
            }
        }

        if (accession) {
            // Clean and standardize the accession tag
            accession = accession.toString().trim().split('/').pop().trim();

            // Show the universal modal
            showAccessionModal(accession);
        }

        return window.dash_clientside.no_update;
    }
    """,
    Output(
        "dummy-output", "children"
    ),  # Dummy output since we don't need to update any Dash components
    [
        Input("dl-table", "cellClicked"),
        Input("dl-table", "active_cell"),
    ],
    [
        State("dl-table", "derived_virtual_data"),
        State("dl-table", "page_current"),
        State("dl-table", "page_size"),
    ],
    prevent_initial_call=True,
)

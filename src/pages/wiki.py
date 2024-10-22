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
import os

from src.components.sqlite import engine
from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot, make_logo

dash.register_page(__name__)


layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="meta-data"),
        dcc.Store(id="paper-data"),
        dcc.Store(id="active-item-cache"),
        dmc.Grid(
            justify="center",
            align="begin",
            style={"paddingTop": "20px"},
            children=[
                dmc.GridCol(
                    span={"lg": 6, "sm": 12},
                    style={
                        "overflow": "hidden",
                        "padding": "10px",
                    },
                    children=[
                        html.H1(
                            "Summary and characteristics for each Starship family",
                        ),
                        dcc.Loading(
                            id="wiki-loading",
                            type="circle",
                            children=html.Div(
                                id="accordion",
                                children=dbc.Accordion(
                                    id="category-accordion",
                                    children=[
                                        dbc.AccordionItem("Phoenix", item_id="Phoenix"),
                                        dbc.AccordionItem(
                                            "Hephaestus", item_id="Hephaestus"
                                        ),
                                        dbc.AccordionItem("Tardis", item_id="Tardis"),
                                        dbc.AccordionItem(
                                            "Serenity", item_id="Serenity"
                                        ),
                                        dbc.AccordionItem(
                                            "Prometheus", item_id="Prometheus"
                                        ),
                                        dbc.AccordionItem(
                                            "Enterprise", item_id="Enterprise"
                                        ),
                                        dbc.AccordionItem(
                                            "Galactica", item_id="Galactica"
                                        ),
                                        dbc.AccordionItem("Moya", item_id="Moya"),
                                        dbc.AccordionItem("Arwing", item_id="Arwing"),
                                        dbc.AccordionItem("Voyager", item_id="Voyager"),
                                        dbc.AccordionItem(
                                            "Family-11", item_id="Family-11"
                                        ),
                                    ],
                                ),
                            ),
                        ),
                        dcc.Loading(
                            id="loading-content",
                            children=html.Div(id="accordion-content"),
                        ),
                    ],
                ),
                dmc.GridCol(
                    span={"lg": 6, "sm": 12},
                    style={
                        "overflow": "hidden",
                        "padding": "10px",
                    },
                    children=[
                        html.Div(id="sidebar-title"),
                        dcc.Loading(
                            id="sidebar-loading",
                            type="circle",
                            children=[
                                dmc.Center(
                                    html.Div(
                                        id="sidebar",
                                        style={
                                            "width": "100%",
                                            "overflowY": "auto",
                                            "overflowX": "auto",
                                        },
                                    )
                                )
                            ],
                        ),
                    ],
                ),
            ],
        ),
    ],
)


def clean_contigIDs(string):
    # Removing omes from contigIDs
    if string is None:
        return None
    else:
        parts = string.split("_", 1)
        if len(parts) > 1:
            prefix = parts[0]
            suffix = parts[1]
            if 7 <= len(prefix) <= 9:
                return suffix
            else:
                return string


def get_logo_or_generate(category, upDRs, downDRs):
    uplogo_img_path = f"assets/images/DR/{category}-upDR.png"
    downlogo_img_path = f"assets/images/DR/{category}-downDR.png"

    if not os.path.exists(uplogo_img_path):
        uplogo_img_path = make_logo(upDRs, uplogo_img_path)
    if not os.path.exists(downlogo_img_path):
        downlogo_img_path = make_logo(downDRs, downlogo_img_path)

    return uplogo_img_path, downlogo_img_path


def create_accordion_item(df, papers, category):
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

    # Generate content
    accordion_content = [
        dash_table.DataTable(
            columns=[
                {"name": "Metric", "id": "metric"},
                {"name": "Value", "id": "value"},
            ],
            data=[
                {
                    "metric": f"Total Starships in {category}",
                    "value": f"{n_ships:,.0f}",
                },
                {"metric": "Max Starship Size (bp)", "value": f"{max_size:,.0f}"},
                {"metric": "Min Starship Size (bp)", "value": f"{min_size:,.0f}"},
            ],
            style_table={"overflowX": "auto"},
            style_cell={"textAlign": "left"},
            style_header={"fontWeight": "bold"},
        )
    ]

    # Add paper references, if available
    if type_element_reference.size > 0:
        paper_links = [
            html.A(ref, href=url, target="_blank")
            for ref, url in zip(
                type_element_reference, filtered_papers_df["Url"].dropna().unique()
            )
        ]
        accordion_content.append(
            html.H5(["Reference for defining type element: ", *paper_links])
        )

    # Sequence logos
    if upDRs and downDRs:
        uplogo_img_path, downlogo_img_path = get_logo_or_generate(
            category, upDRs, downDRs
        )
        accordion_content.append(
            dbc.Row(
                [
                    dbc.Col(html.Img(src=uplogo_img_path)),
                    dbc.Col(html.Img(src=downlogo_img_path)),
                ]
            )
        )

    return accordion_content


@callback(Output("meta-data", "data"), Input("url", "href"))
def load_meta_data(url):
    if url:
        meta_query = """
        SELECT j.ship_family_id, j.curated_status, j.taxid, j.ship_id, j.genome_id, j.ome, j.size, j.upDR, j.downDR, f.familyName, j.contigID, j.elementBegin, j.elementEnd, j."size", f.type_element_reference, a.accession_tag, t."order", t.family, t.genus, t.species, g.version, g.genomeSource, g.citation
        FROM joined_ships j
        LEFT JOIN taxonomy t ON j.taxid = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        WHERE j.orphan IS NULL
        """
        meta_df = pd.read_sql_query(meta_query, engine)

        # clean ome from contigID
        if "contigID" in meta_df.columns:
            meta_df["contigID"] = meta_df["contigID"].apply(clean_contigIDs)

        return meta_df.to_dict("records")


@callback(Output("paper-data", "data"), Input("url", "href"))
def load_paper_data(url):
    if url:
        paper_query = """
        SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, p.shortCitation, f.familyName, f.type_element_reference
        FROM papers p
        LEFT JOIN family_names f ON p.shortCitation = f.type_element_reference
        """
        paper_df = pd.read_sql_query(paper_query, engine)
        return paper_df.to_dict("records")


# Callback to load content only when an accordion item is expanded
@callback(
    Output("accordion-content", "children"),
    [Input("accordion", "active_item")],
    [State("meta-data", "data"), State("paper-data", "data")],
)
def load_accordion_content(active_item, cached_meta, cached_papers):
    if not active_item:
        return "Please expand an accordion item to see the content."

    # Convert the cached data back to DataFrame
    df = pd.DataFrame(cached_meta)
    papers = pd.DataFrame(cached_papers)

    # Generate the content for the active accordion item (lazy loading)
    if active_item == "cat_a":
        category = "Category A"
    elif active_item == "cat_b":
        category = "Category B"
    elif active_item == "cat_c":
        category = "Category C"
    else:
        return "No data available."

    # Call your create_accordion_item function for the active category
    accordion_item = create_accordion_item(df, papers, category)

    # Return the generated accordion item content
    return accordion_item


@callback(
    Output("accordion", "children"),
    [Input("meta-data", "data"), Input("paper-data", "data")],
)
def create_accordion(cached_meta, cached_papers):
    df = pd.DataFrame(cached_meta)
    papers = pd.DataFrame(cached_papers)
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
        active_item="Voyager",
    )


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
        raise PreventUpdate  # Prevent the callback from running if the inputs are invalid

    df = pd.DataFrame(cached_meta)
    title = html.H1(f"Taxonomy and Genomes for Starships in {active_item}")

    filtered_meta_df = df[df["familyName"] == active_item].sort_values(
        by="accession_tag", ascending=False
    )
    sunburst = create_sunburst_plot(df=filtered_meta_df, type="tax", title_switch=False)
    fig = dcc.Graph(figure=sunburst, style={"width": "100%", "height": "100%"})

    # TODO: add links to genome browser, and extra classification info for selected genome

    table_columns = [
        {
            "name": "Accession",
            "id": "accession_tag",
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
            "name": "Contig ID",
            "id": "contigID",
            "deletable": False,
            "selectable": False,
        },
        {
            "name": "Start Position in Genome",
            "id": "elementBegin",
            "deletable": False,
            "selectable": False,
        },
        {
            "name": "End Position in Genome",
            "id": "elementEnd",
            "deletable": False,
            "selectable": False,
        },
        {
            "name": "Element Length",
            "id": "size",
            "deletable": False,
            "selectable": False,
        },
        # {
        #     "name": "Source",
        #     "id": "genomeSource",
        #     "deletable": False,
        #     "selectable": False,
        # },
        # {
        #     "name": "Citation/Release Date",
        #     "id": "citation",
        #     "deletable": False,
        #     "selectable": False,
        # },
    ]
    table = make_ship_table(
        filtered_meta_df, id="wiki-table", columns=table_columns, pg_sz=15
    )

    output = dbc.Stack(children=[fig, table], direction="vertical", gap=5)

    return output, title, active_item
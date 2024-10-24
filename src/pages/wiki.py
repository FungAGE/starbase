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

from src.components.mariadb import engine
from src.components.tables import make_ship_table
from src.utils.plot_utils import create_sunburst_plot, make_logo

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
                                children=dbc.Accordion(id="category-accordion"),
                            ),
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


@callback(Output("meta-data", "data"), Input("url", "href"))
def load_meta_data(url):
    if url:
        meta_query = """
        SELECT 
            j.ship_family_id, 
            j.curated_status, 
            j.taxid, 
            j.ship_id, 
            j.genome_id, 
            j.ome, 
            j.size, 
            j.upDR, 
            j.downDR, 
            f.familyName, 
            j.contigID, 
            j.elementBegin, 
            j.elementEnd, 
            j.size, 
            f.type_element_reference, 
            a.accession_tag, 
            t.`order`,
            t.family, 
            t.genus, 
            t.species, 
            g.version, 
            g.genomeSource, 
            g.citation
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

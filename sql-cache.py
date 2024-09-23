import dash
from dash import dcc, html, callback, no_update, dash_table
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
from flask_caching import Cache
import pandas as pd
from sqlalchemy import create_engine

from src.utils.plot_utils import create_sunburst_plot
from src.utils.tree import plot_tree

CACHE_CONFIG = {
    "CACHE_TYPE": "filesystem",
    "CACHE_DIR": "/tmp/dash_cache",
    "CACHE_DEFAULT_TIMEOUT": 300,  # optional, but good to specify
}

external_stylesheets = [
    # dbc.themes.BOOTSTRAP,
    "https://codepen.io/chriddyp/pen/bWLwgP.css",
    "https://codepen.io/chriddyp/pen/brPBPO.css",
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Initialize the cache with the Flask server instance
cache = Cache(app.server, config=CACHE_CONFIG)

app.layout = html.Div(
    [
        dbc.Button(
            "Reset",
            id="reset-button",
            n_clicks=0,
            color="success",
            class_name="mx-auto",
        ),
        dcc.Dropdown(
            id="dropdown",
            options=[
                {"label": "curated", "value": "curated"},
                {"label": "all", "value": "all"},
            ],
            value="all",
        ),
        html.Div(
            [
                dcc.Loading(
                    id="loading-1",
                    type="circle",  # or "dot", "default"
                    children=[
                        html.Div(dcc.Graph(id="pie-1"), className="six columns"),
                    ],
                ),
                dcc.Loading(
                    id="loading-2",
                    type="circle",  # or "dot", "default"
                    children=[
                        html.Div(dcc.Graph(id="pie-2"), className="six columns"),
                    ],
                ),
            ],
            className="row",
        ),
        html.Div(
            [
                html.Div(
                    dash_table.DataTable(
                        id="explore-table",
                        data=[],
                        columns=[
                            {
                                "name": i,
                                "id": i,
                                "deletable": False,
                                "selectable": True,
                            }
                            for i in [
                                "accession_tag",
                                "familyName",
                                "order",
                                "family",
                                "genus",
                                "species",
                            ]
                        ],
                        selected_columns=[],
                        selected_rows=[],
                        editable=False,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        row_selectable="multi",
                        row_deletable=False,
                        page_action="native",
                        page_current=0,
                        page_size=10,
                    ),
                    className="six columns",
                ),
                html.Div(dcc.Graph(id="tree"), className="six columns"),
            ],
            className="row",
        ),
        dcc.Location(id="url", refresh=False),
        dcc.Store("curated-status"),
    ]
)


@cache.memoize()
def global_store(value):
    print(f"Computing {value} ships")
    query = """
    SELECT j.*, t."order", t.family, f.longFamilyID, f.familyName, a.accession_tag
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN accessions a ON j.ship_id = a.id
    """
    df = pd.read_sql_query(query, engine)
    if value == "curated":
        df_filtered = df[df["curated_status"] == "curated"]
    else:
        df_filtered = df
    unique_df = df_filtered.drop_duplicates(subset=["accession_tag"])
    return unique_df


@app.callback(Output("curated-status", "data"), Input("dropdown", "value"))
def compute_value(value):
    global_store(value)
    return value


@cache.memoize()
def update_graph_1(data):
    return create_sunburst_plot(df=data, type="ship")


@cache.memoize()
def update_graph_2(data):
    return create_sunburst_plot(df=data, type="tax")


@cache.memoize()
def update_ship_table_data(data):
    return data.to_dict("records")


@cache.memoize()
def update_phylogeny(highlight_families="all"):
    tree = plot_tree(highlight_families)
    return tree


@app.callback(
    [
        Output("pie-1", "figure"),
        Output("pie-2", "figure"),
        Output("tree", "figure"),
        Output("explore-table", "data"),
    ],
    [
        Input("url", "href"),
        Input("reset-button", "n_clicks"),
        Input("pie-1", "clickData"),
        Input("pie-2", "clickData"),
        Input("tree", "clickData"),
        Input("explore-table", "derived_virtual_data"),
        Input("explore-table", "derived_virtual_selected_rows"),
        Input("curated-status", "data"),
    ],
)
def update_ui(
    url,
    n_clicks,
    clickData1,
    clickData2,
    clickData_phylo,
    table_data,
    table_rows,
    value,
):
    df = pd.DataFrame(global_store(value))

    # Initialize outputs to 'no_update' to avoid unnecessary updates
    pie1_output = no_update
    pie2_output = no_update
    tree_output = no_update
    table_output = no_update

    # Check if reset or URL navigation triggers update
    if url or n_clicks:
        pie1_output = update_graph_1(df).to_dict()
        pie2_output = update_graph_2(df).to_dict()
        tree_output = update_phylogeny().to_dict()
        table_output = df.to_dict("records")

    # Handle clicks on Pie 1
    if clickData1:
        path1 = [point["label"] for point in clickData1["points"]]
        plot1_df = df[df["familyName"].isin(path1)]
        print(plot1_df.head())

        pie2_output = update_graph_2(plot1_df)
        tree_output = update_phylogeny(plot1_df["familyName"].unique())
        table_output = update_ship_table_data(plot1_df)

    # Handle clicks on Pie 2
    if clickData2:
        path2 = [point["label"] for point in clickData2["points"]]
        plot2_df = df[df["order"].isin(path2)]
        print(plot2_df.head())

        pie1_output = update_graph_1(plot2_df)
        tree_output = update_phylogeny(plot2_df["familyName"].unique())
        table_output = update_ship_table_data(plot2_df)

    # Handle clicks on the Phylogeny tree
    if clickData_phylo:
        selected_clades = [clickData_phylo["points"][0]["text"]]
        print(selected_clades)
        plot_df = df[df["familyName"].isin(selected_clades)]
        print(plot_df.head())

        pie1_output = update_graph_1(plot_df)
        pie2_output = update_graph_2(plot_df)
        table_output = update_ship_table_data(plot_df)

    # Handle row selections from the table
    if table_rows:
        table_df = pd.DataFrame(table_data).iloc[table_rows]
        print(table_df.head())
        selected_clades = table_df["familyName"].tolist()
        selected_ships = table_df["accession_tag"]
        plot_df = df[df["accession_tag"].isin(selected_ships)]

        pie1_output = update_graph_1(plot_df)
        pie2_output = update_graph_2(plot_df)
        tree_output = update_phylogeny(selected_clades)

    # Return the updated or unchanged figures and table
    return pie1_output, pie2_output, tree_output, table_output


if __name__ == "__main__":
    engine = create_engine("sqlite:///database_folder/starbase.sqlite")
    query = "SELECT name FROM sqlite_master WHERE type='table'"
    sql_tbls = pd.read_sql_query(query, engine)
    app.run_server(debug=True, threaded=True, host="localhost", port=8005)

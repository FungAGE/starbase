import warnings

warnings.filterwarnings("ignore")

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag

from src.config.cache import cache
from src.database.sql_manager import fetch_paper_data


def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


# Function to convert URL string to HTML link
def url_to_link(url, label):
    return f"[{label}]({url})"


def make_ship_table(df, id, columns=None, select_rows=False, pg_sz=None):
    print("Data shape:", df.shape)
    print("Columns:", df.columns.tolist())
    print("First row:", df.iloc[0].to_dict() if not df.empty else "Empty")
    print("Column definitions:", columns)
    
    if pg_sz is None:
        pg_sz = 10
        
    defaultColDef = {
        "resizable": True,
        "sortable": True,
        "filter": True,
        "minWidth": 100,
        "flex": 1
    }
    
    # If no columns provided, create them from DataFrame columns
    if columns is None:
        grid_columns = []
        if select_rows:
            grid_columns.append({
                "headerCheckboxSelection": True,
                "checkboxSelection": True,
                "width": 50,
                "pinned": "left",
                "lockPosition": True,
                "suppressMenu": True,
                "headerName": "",
                "flex": 0
            })
        grid_columns.extend([
            {
                "field": col,  # Use column name directly from DataFrame
                "headerName": col.replace("_", " ").title(),  # Create readable header
                "flex": 1
            }
            for col in df.columns
        ])
    else:
        # Ensure columns are in AG Grid format
        grid_columns = []
        if select_rows:
            grid_columns.append({
                "headerCheckboxSelection": True,
                "checkboxSelection": True,
                "width": 50,
                "pinned": "left",
                "lockPosition": True,
                "suppressMenu": True,
                "headerName": "",
                "flex": 0
            })
        grid_columns.extend([
            {
                "field": col.get("field"),  # Use the field directly
                "headerName": col.get("headerName"),
                "flex": 1
            }
            for col in columns
        ])
    
    print("Final grid columns:", grid_columns)
    print("Row data sample:", df.head(1).to_dict("records") if not df.empty else "Empty")
    
    grid = dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=df.to_dict("records") if df is not None else [],
        defaultColDef=defaultColDef,
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": pg_sz,
            "rowSelection": "multiple" if select_rows else None,
        },
        className="ag-theme-alpine",
        style={"height": "60vh", "width": "100%"}
    )
    return html.Div(grid)

def make_paper_table():
    df = cache.get("paper_data")
    if df is None:
        df = fetch_paper_data()
    if df is not None:
        df_summary = (
            df.groupby("Title")
            .agg({
                "familyName": lambda x: ", ".join(sorted(filter(None, x.unique()))),
                "Author": "first",
                "PublicationYear": "first",
                "DOI": "first",
                "Url": "first",
            })
            .reset_index()
        )

        sub_df = df_summary.sort_values(by="PublicationYear", ascending=False)

        # sub_df["Title"] = sub_df["Title"].apply(lambda x: truncate_string(x, length=40))
        # sub_df["Author"] = sub_df["Author"].apply(lambda x: truncate_string(x, length=40))
        sub_df["DOI"] = sub_df["DOI"].apply(lambda x: url_to_link(x, label=x))
        sub_df["Url"] = sub_df["Url"].apply(lambda x: url_to_link(x, label="full text"))

        # Modify columns for mobile responsiveness
        sub_df_columns = [
            {
                "name": "Title",
                "id": "Title",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
            {
                "name": "Publication Year",
                "id": "PublicationYear",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
            {
                "name": "Authors",
                "id": "Author",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
            {
                "name": "DOI",
                "id": "DOI",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
        ]

        # Create card view for mobile
        cards = dmc.Stack(
            children=[
                dmc.Paper(
                    children=[
                        dmc.Title(row["Title"], order=4, mb="sm"),
                        dmc.Text(f"Authors: {row['Author']}", size="sm", mb="xs"),
                        dmc.Text(f"Year: {row['PublicationYear']}", size="sm", mb="xs"),
                        dmc.Group([
                            dmc.Anchor("DOI", href=row["DOI"].split("](")[1][:-1], size="sm"),
                            dmc.Anchor("Full Text", href=row["Url"].split("](")[1][:-1], size="sm"),
                        ]),
                    ],
                    p="md",
                    radius="md",
                    withBorder=True,
                )
                for _, row in sub_df.iterrows()
            ],
            gap="md",
            style={"display": "none"},
            className="mobile-cards"
        )

        # Create table view (your existing table code)
        table = html.Div(
            dash_table.DataTable(
                data=sub_df.to_dict("records"),
                columns=sub_df_columns,
                id="papers-table",
                markdown_options={"html": True},
                style_table={
                    "overflowX": "auto",
                    "overflowY": "auto",
                    "maxHeight": "60vh",
                    "minWidth": "300px",  # Ensure minimum width on mobile
                },
                style_data={
                    "height": "auto",
                    "lineHeight": "20px",
                    "padding": "10px",
                    "whiteSpace": "normal",  # Allow text wrapping
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                },
                style_cell={
                    "fontFamily": "Arial, sans-serif",
                    "textAlign": "left",
                    "minWidth": "100px",
                    "maxWidth": {  # Responsive column widths
                        "Title": "300px",
                        "PublicationYear": "100px",
                        "Author": "200px",
                        "DOI": "150px",
                    },
                    "width": {  # Default widths for different columns
                        "Title": "300px",
                        "PublicationYear": "100px",
                        "Author": "200px",
                        "DOI": "150px",
                    },
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                },
                style_cell_conditional=[  # Hide certain columns on small screens
                    {
                        "if": {"column_id": "Author"},
                        "@media screen and (max-width: 768px)": {"display": "none"},
                    },
                    {
                        "if": {"column_id": "PublicationYear"},
                        "@media screen and (max-width: 480px)": {"display": "none"},
                    },
                ],
                style_header={
                    "backgroundColor": "#f8f9fa",
                    "fontWeight": "bold",
                    "borderBottom": "2px solid #dee2e6",
                    "textAlign": "left",
                    "padding": "12px",
                },
                style_filter={
                    "backgroundColor": "#f8f9fa",
                    "padding": "8px",
                },
                style_data_conditional=[
                    {
                        "if": {"row_index": "odd"},
                        "backgroundColor": "#f8f9fa",
                    },
                    {
                        "if": {"state": "selected"},
                        "backgroundColor": "#e3f2fd",
                        "border": "1px solid #2196f3",
                    },
                ],
            ),
            className="desktop-table"
        )

        # Add CSS to handle view switching
        return html.Div([
            # View toggle button
            dmc.Group(
                [
                    dmc.Button(
                        "Toggle View",
                        id="toggle-paper-view",
                        variant="outline",
                        size="sm",
                        leftSection=DashIconify(icon="tabler:layout-list"),
                        className="mobile-only"
                    )
                ],
                pos="right",
            mb="md"
            ),
            cards,
            table
        ])

def make_ship_blast_table(ship_blast_results, id, df_columns):
    return html.Div(dash_table.DataTable(
        columns=df_columns,
        data=ship_blast_results.to_dict("records"),
        id=id,
        editable=False,
        markdown_options=None,
        sort_action="native",
        sort_by=[{"column_id": "evalue", "direction": "asc"}],
        sort_mode="single",
        row_selectable="single",
        selected_rows=[0],
        row_deletable=False,
        page_action="native",
        page_current=0,
        page_size=10,
        derived_virtual_selected_rows=[0],
        derived_virtual_indices=[],
        derived_virtual_data=ship_blast_results.to_dict("records"),
        style_table={
            "overflowX": "auto",
            "overflowY": "auto",
            "maxHeight": "60vh",
        },
        style_data={
            "height": "auto",
            "lineHeight": "20px",
            "padding": "10px",
        },
        style_cell={
            "fontFamily": "Arial, sans-serif",
            "textAlign": "left",
            "minWidth": "100px",
            "maxWidth": "300px",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
        },
        style_header={
            "backgroundColor": "#f8f9fa",
            "fontWeight": "bold",
            "borderBottom": "2px solid #dee2e6",
            "textAlign": "left",
            "padding": "12px",
        },
        style_filter={
            "backgroundColor": "#f8f9fa",
            "padding": "8px",
        },
        style_data_conditional=[
            {
                "if": {"row_index": "odd"},
                "backgroundColor": "#f8f9fa",
            },
            {
                "if": {"state": "selected"},
                "backgroundColor": "#e3f2fd",
                "border": "1px solid #2196f3",
            },
            {
                "if": {"column_id": "accession_tag"},
                "color": "blue",
                "textDecoration": "underline",
                "cursor": "pointer",
            }
        ],
    ))

def make_dl_table(df, id, table_columns):
    if isinstance(df, pd.DataFrame):
        data = df.to_dict("records")
    elif isinstance(df, list):
        data = df
    else:
        data = []
        
    return html.Div(dash_table.DataTable(
        id=id,
        columns=table_columns,
        data=data,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="multi",
        page_action="native",
        page_current=0,
        page_size=20,
        cell_selectable=True,
        style_data={
            "height": "auto",
            "lineHeight": "20px",
            "padding": "10px",
        },
        style_cell={
            "fontFamily": "Arial, sans-serif",
            "textAlign": "left",
            "minWidth": "100px",
            "maxWidth": "300px",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
        },
        style_header={
            "backgroundColor": "#f8f9fa",
            "fontWeight": "bold",
            "borderBottom": "2px solid #dee2e6",
            "textAlign": "left",
            "padding": "12px",
        },
        style_filter={
            "backgroundColor": "#f8f9fa",
            "padding": "8px",
        },
        style_data_conditional=[
            {
                "if": {"row_index": "odd"},
                "backgroundColor": "#f8f9fa",
            },
            {
                "if": {"state": "selected"},
                "backgroundColor": "#e3f2fd",
                "border": "1px solid #2196f3",
            },
            {
                "if": {"column_id": "accession_tag"},
                "color": "blue",
                "textDecoration": "underline",
                "cursor": "pointer",
            }
        ],
    )
    )

def make_wiki_table(n_ships, max_size, min_size):
    """Create a summary table for a Starship family."""
    data = [
        {
            "metric": "Total Number of Starships",
            "value": f"{n_ships:,.0f}",
        },
        {
            "metric": "Maximum Size (bp)",
            "value": f"{max_size:,.0f}",
        },
        {
            "metric": "Minimum Size (bp)",
            "value": f"{min_size:,.0f}",
        },
    ]

    return html.Div(
        dash_table.DataTable(
            data=data,
            columns=[
                {"name": col, "id": col} for col in ["metric", "value"]
            ],
            style_table={"overflowX": "auto", "maxWidth": "500px"},
            style_cell={"textAlign": "left"},
            style_header={"fontWeight": "bold"},
        )
    )
import warnings

warnings.filterwarnings("ignore")

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd

from src.components.sql_manager import load_from_cache
from src.components.sql_manager import fetch_paper_data


def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


# Function to convert URL string to HTML link
def url_to_link(url, label):
    return f"[{label}]({url})"


def make_ship_table(df, id, columns=None, select_rows=False, pg_sz=None):
    if pg_sz is None:
        pg_sz = 10

    if columns:
        table_columns = columns
    else:
        table_columns = [
            {
                "name": i,
                "id": i,
                "deletable": False,
                "selectable": select_rows,
            }
            for i in df.columns
        ]

    if df is not None:
        table_df = df.to_dict("records")
    else:
        table_df = []
    table = html.Div(dash_table.DataTable(
        id=id,
        columns=table_columns,
        data=df.to_dict("records"),
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable=select_rows,
        page_action="native",
        page_current=0,
        page_size=20,
        cell_selectable=True,  # Added to enable cell clicking
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
    return table


def make_paper_table():
    df = load_from_cache("paper_data")
    if df is None:
        df = fetch_paper_data()
    if df is not None:
        df_summary = (
            df.groupby("Title")
            .agg(
                {
                    "familyName": lambda x: ", ".join(sorted(filter(None, x.unique()))),
                    "Author": "first",
                    "PublicationYear": "first",
                    "DOI": "first",
                    "Url": "first",
                }
            )
            .reset_index()
        )

        sub_df = df_summary.sort_values(by="PublicationYear", ascending=False)

        # sub_df["Title"] = sub_df["Title"].apply(lambda x: truncate_string(x, length=40))
        # sub_df["Author"] = sub_df["Author"].apply(lambda x: truncate_string(x, length=40))
        sub_df["DOI"] = sub_df["DOI"].apply(lambda x: url_to_link(x, label=x))
        sub_df["Url"] = sub_df["Url"].apply(lambda x: url_to_link(x, label="full text"))

        # rename columns
        sub_df_columns = [
            {
                "name": "Starship Families Described",
                "id": "familyName",
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
                "name": "Title",
                "id": "Title",
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
            {
                "name": "Url",
                "id": "Url",
                "deletable": False,
                "selectable": False,
                "presentation": "markdown",
            },
        ]

        paper_table = html.Div(dash_table.DataTable(
            data=sub_df.to_dict("records"),
            sort_action="none",
            columns=sub_df_columns,
            id="papers-table",
            markdown_options={"html": True},
            style_table={
                "overflowX": "auto",
                "overflowY": "auto",
                "maxHeight": "60vh",
            },
            style_data={
                "height": "auto",
                "lineHeight": "20px",
                "padding": "10px",
                "whiteSpace": "normal",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
            },
            style_cell={
                "fontFamily": "Arial, sans-serif",
                "textAlign": "left",
                "minWidth": "120px",
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
            ],
        ))
        return paper_table

def make_ship_blast_table(ship_blast_results, id,df_columns):
    return html.Div(dash_table.DataTable(
        columns=df_columns,
        data=ship_blast_results.to_dict("records"),
        id=id,
        editable=False,
        sort_action="native",
        sort_by=[{"column_id": "evalue", "direction": "asc"}],
        sort_mode="single",
        row_selectable="single",
        selected_rows=[0],
        row_deletable=False,
        selected_columns=[],
        page_action="native",
        page_current=0,
        page_size=10,
        export_format="tsv",
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
import warnings

warnings.filterwarnings("ignore")

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd

from src.components.cache_manager import load_from_cache
from src.components.sql_queries import fetch_paper_data


def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


# Function to convert URL string to HTML link
def url_to_link(url, label):
    return f"[{label}]({url})"


def make_ship_table(df, id, columns=None, pg_sz=None):
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
                "selectable": True,
            }
            for i in df.columns
        ]

    if df is not None:
        table_df = df.to_dict("records")
    else:
        table_df = []
    table = dash_table.DataTable(
        id="pgv-table",
        columns=table_columns,
        data=df.to_dict("records"),
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="multi",
        page_action="native",
        page_current=0,
        page_size=20,
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
        ],
    )
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

        paper_table = dash_table.DataTable(
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
        )
        return paper_table

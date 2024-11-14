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
        id=id,
        columns=table_columns,
        data=table_df,
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        row_selectable="multi",
        row_deletable=False,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        cell_selectable=True,
        active_cell=None,
        page_current=0,
        page_size=pg_sz,
        style_table={
            "width": "100%",
            "height": "100%",
            "overflowX": "auto",
            "overflowY": "auto",
        },
        style_data={
            "whiteSpace": "minimal",
        },
        style_data_conditional=[
            {
                "if": {"column_id": "accession_tag"},
                "color": "blue",
                "textDecoration": "underline",
                "cursor": "pointer",
            }
        ],
        style_cell={
            "minWidth": "0px",
            "maxWidth": "100%",
            "textAlign": "left",
        },
        style_header={
            "backgroundColor": "lightgrey",
            "fontWeight": "bold",
            "textAlign": "left",
        },
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

        paper_table = dbc.Card(
            [
                dbc.CardHeader(
                    html.Div(
                        ["Manuscripts Characterizing Starships"],
                        className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                    )
                ),
                dbc.CardBody(
                    [
                        dash_table.DataTable(
                            data=sub_df.to_dict("records"),
                            sort_action="none",
                            columns=sub_df_columns,
                            id="papers-table",
                            markdown_options={"html": True},
                            style_table={
                                "overflowX": "auto",
                            },
                            style_data={
                                "height": "auto",
                                "whiteSpace": "normal",
                                "overflow": "hidden",
                                "textOverflow": "ellipsis",
                            },
                            style_cell={
                                "minWidth": "120px",
                                "maxWidth": "300px",
                                "textAlign": "left",
                                "padding": "5px",
                            },
                            style_header={
                                "backgroundColor": "lightgrey",
                                "fontWeight": "bold",
                                "textAlign": "left",
                            },
                        ),
                    ]
                ),
            ],
            # className="auto-resize-900",
        )

        return paper_table

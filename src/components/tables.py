import warnings

warnings.filterwarnings("ignore")

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd


def truncate_string(s, length=10):
    return s if len(s) <= length else s[:length] + "..."


# Function to convert URL string to HTML link
def url_to_link(url, label):
    return f'<a href="{url}" target="_blank">{label}</a>'


def make_ship_table(df, id, columns=None, pg_sz=None):
    if pg_sz is None:
        pg_sz = 10
    if columns is None:
        # columns = df.columns
        columns = [
            "accession_tag",
            "familyName",
            "starship_navis",
            "starship_haplotype",
            "genus",
            "species",
        ]

    table = html.Div(
        [
            dash_table.DataTable(
                id=id,
                columns=[
                    {
                        "name": i,
                        "id": i,
                        "deletable": False,
                        "selectable": True,
                    }
                    for i in columns
                ],
                data=df.to_dict("records"),
                editable=False,
                filter_action="native",
                sort_action="native",
                sort_mode="multi",
                row_selectable="multi",
                row_deletable=False,
                selected_columns=[],
                selected_rows=[],
                page_action="native",
                page_current=0,
                page_size=pg_sz,
                style_table={
                    "overflowX": "auto",
                },
                style_cell={
                    "minWidth": "150px",
                    "width": "150px",
                    "maxWidth": "150px",
                    "whiteSpace": "normal",
                },
            ),
        ],
        className="auto-resize-750",
    )
    return table


def make_paper_table(data):

    # df = pd.read_csv(data)
    # sub_df = df[["Title", "Author", "PublicationYear", "DOI", "Url"]].sort_values(
    #     by="PublicationYear", ascending=False
    # )

    df = pd.DataFrame(data)

    df_summary = (
        df.groupby("Title")
        .agg(
            {
                "familyName": lambda x: ", ".join(sorted(x.unique())),
                "Author": "first",
                "PublicationYear": "first",
                "DOI": "first",
                "Url": "first",
            }
        )
        .reset_index()
    )

    sub_df = df_summary.sort_values(by="PublicationYear", ascending=False)

    sub_df["Title"] = sub_df["Title"].apply(lambda x: truncate_string(x, length=20))
    sub_df["Author"] = sub_df["Author"].apply(lambda x: truncate_string(x, length=20))
    sub_df["DOI"] = sub_df["DOI"].apply(lambda x: url_to_link(x, label=x))
    sub_df["Url"] = sub_df["Url"].apply(lambda x: url_to_link(x, label="full text"))

    # rename columns
    sub_df_columns = [
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
            "name": "Publication Year",
            "id": "PublicationYear",
            "deletable": False,
            "selectable": False,
            "presentation": "markdown",
        },
        {
            "name": "Starship Families Described",
            "id": "familyName",
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
                            "overflowX": "auto",  # Enable horizontal scrolling
                        },
                        style_data={
                            "whiteSpace": "normal",  # Allow text to wrap
                            "overflow": "hidden",
                            "textOverflow": "ellipsis",
                        },
                        style_cell={
                            "minWidth": "120px",  # Minimum column width
                            "maxWidth": "100%",  # Maximum column width
                            "textAlign": "left",  # Align text to the left
                            "padding": "5px",  # Add padding for better spacing
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

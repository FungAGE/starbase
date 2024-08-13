import dash
from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd


def truncate_string(s, length=10):
    return s if len(s) <= length else s[:length] + "..."


df = pd.read_csv("src/data/papers.csv")
sub_df = df[["Title", "Author", "PublicationYear", "DOI", "Url"]].sort_values(
    by="PublicationYear", ascending=False
)


# Function to convert URL string to HTML link
def url_to_link(url):
    return f'<a href="{url}" target="_blank">{url}</a>'


sub_df["Url"] = sub_df["Url"].apply(url_to_link)
sub_df["Author"] = sub_df["Author"].apply(lambda x: truncate_string(x, length=20))

# TODO: rename columns

table = dbc.Card(
    [
        dbc.CardHeader(
            html.H3(
                "Manuscripts Characterizing Starships",
            )
        ),
        dbc.CardBody(
            [
                html.P(
                    [
                        dash_table.DataTable(
                            data=sub_df.to_dict("records"),
                            sort_action="none",
                            columns=[
                                {
                                    "name": i,
                                    "id": i,
                                    "deletable": False,
                                    "selectable": False,
                                    "presentation": "markdown",
                                }
                                for i in sub_df.columns
                            ],
                            id="papers-table",
                            markdown_options={"html": True},
                            style_data={
                                "whiteSpace": "normal",
                                "height": "auto",
                                "minWidth": "150px",
                                "width": "150px",
                                "maxWidth": "150px",
                                "overflow": "hidden",
                                "overflowX": "auto",
                                "textOverflow": "ellipsis",
                            },
                        )
                    ],
                    className="box-body",
                ),
            ]
        ),
    ]
)

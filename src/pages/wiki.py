import dash
from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd

dash.register_page(__name__)


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

layout = dbc.Container(
    fluid=True,
    children=[
        dbc.Stack(
            [
                dbc.Row(
                    justify="center",
                    align="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=12,
                            style={"padding": "20px"},
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            html.H2(
                                                ["What is a Starship?"],
                                            ),
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.P(
                                                    [
                                                        "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes). They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in ",
                                                        html.Span(
                                                            "Paecilomyces",
                                                            style={
                                                                "font-style": "italic"
                                                            },
                                                        ),
                                                        " cheese making in ",
                                                        html.Span(
                                                            "Penicillium",
                                                            style={
                                                                "font-style": "italic",
                                                            },
                                                        ),
                                                        ", and enable the transfer of formaldehyde resistance in ",
                                                        html.Span(
                                                            "Aspergillus nidulans",
                                                            style={
                                                                "font-style": "italic",
                                                            },
                                                        ),
                                                        " and ",
                                                        html.Span(
                                                            "Penicillium chrysogenum.",
                                                            style={
                                                                "font-style": "italic",
                                                            },
                                                        ),
                                                    ],
                                                    style={"fontSize": "1vw"},
                                                ),
                                                html.Br(),
                                                html.Div(
                                                    children=[
                                                        html.Img(
                                                            src="assets/images/starship-model.png",
                                                            style={
                                                                "backgroundColor": "white",
                                                                "margin": "auto",
                                                                "display": "block",
                                                                "width": "85%",
                                                            },
                                                        )
                                                    ],
                                                    className="box-body",
                                                ),
                                            ]
                                        ),
                                    ],
                                    color="primary",
                                    inverse=True,
                                    style={
                                        "width": "100%",
                                    },
                                ),
                            ],
                        )
                    ],
                ),
                dbc.Row(
                    align="center",
                    justify="center",
                    children=[
                        dbc.Col(
                            lg=8,
                            sm=10,
                            children=[
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
                                                        "textOverflow": "ellipsis",
                                                    },
                                                )
                                            ],
                                            className="box-body",
                                        ),
                                    ]
                                ),
                            ],
                        ),
                    ],
                ),
            ],
            gap=3,
        ),
    ],
)

import dash
from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd
import sqlite3

dash.register_page(__name__)

try:
    # Create SQLite database connection
    conn = sqlite3.connect("database_folder/starbase.sqlite")
    c = conn.cursor()

    # Check if the table already exists
    c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='papers'")

except sqlite3.Error as error:
    print("Error connecting to database.", error)


# Read data from SQLite table into a DataFrame
df = pd.read_csv("src/assets/papers.csv")
sub_df = df[["Title", "Author", "PublicationYear", "PublicationTitle", "DOI", "Url"]]


# Function to convert URL string to HTML link
def url_to_link(url):
    return f'<a href="{url}" target="_blank">{url}</a>'


# Apply the function to the 'Website' column
sub_df["Url"] = sub_df["Url"].apply(url_to_link)

# TODO: rename columns

layout = html.Div(
    [
        dbc.Container(
            fluid=True,
            children=[
                dbc.Row(
                    dbc.Col(
                        width=6,
                        className="align-self-center",
                        children=[
                            dbc.Card(
                                [
                                    dbc.CardHeader(
                                        html.H3(
                                            ["What is a Starship?"],
                                            style={"fontSize": "1vw"},
                                        ),
                                    ),
                                    dbc.CardBody(
                                        [
                                            html.P(
                                                [
                                                    "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes) [2]. They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in Paecilomyces, cheese making in Penicillium, and enable the reansfer of formaldehyde resistance in Aspergillus nidulans and Penicillium chrysogenum."
                                                ],
                                                style={"fontSize": "0.6vw"},
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
                    ),
                    justify="center",
                    align="center",
                ),
                html.Br(),
                dbc.Row(
                    dbc.Col(
                        width=6,
                        className="align-self-center",
                        children=[
                            dbc.CardHeader(
                                html.H3(
                                    "Research Papers Characterizing Starships",
                                    style={"fontSize": "1vw"},
                                )
                            ),
                            dbc.CardBody(
                                [
                                    html.Div(
                                        [
                                            dash_table.DataTable(
                                                data=sub_df.to_dict("records"),
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
                                                style_data={
                                                    "whiteSpace": "normal",
                                                    "height": "auto",
                                                },
                                                markdown_options={"html": True},
                                            ),
                                            # html.Div(
                                            #     id="papers-table-interactivity-container"
                                            # ),
                                        ]
                                    ),
                                ]
                            ),
                        ],
                        style={"fontSize": "0.5vw"},
                    ),
                    align="center",
                    justify="center",
                ),
            ],
        )
    ]
)

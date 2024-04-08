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
    c.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='submissions_new'"
    )

except sqlite3.Error as error:
    print("Error connecting to database.", error)


# Read data from SQLite table into a DataFrame
df = pd.read_sql_query(
    "SELECT Title, Author, PublicationYear, PublicationTitle, DOI FROM papers", conn
)

# TODO: rename columns
# TODO: include hyperlinks

layout = html.Div(
    [
        html.Div(
            style={
                "display": "flex",
                "justify-content": "center",
                "align-items": "center",
            },
            children=[
                html.Div(
                    [
                        html.Div(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(html.H3("What is a Starship?")),
                                        dbc.CardBody(
                                            [
                                                html.P(
                                                    "Starships are novel family of class II DNA transposons, endemic to Pezizomycotina. Starships can be extremely large (~20-700kb), making up to 2% of fungal genomes. These elements replicate within the host genome via tyrosine recombinases (captain genes) [2]. They can also pick up and carry relevant genetic 'cargo', including genes for metal resistance in Paecilomyces, cheese making in Penicillium, and enable the reansfer of formaldehyde resistance in Aspergillus nidulans and Penicillium chrysogenum."
                                                ),
                                                html.Br(),
                                                html.Div(
                                                    style={
                                                        "display": "flex",
                                                        "justify-content": "center",
                                                        "align-items": "center",
                                                        "backgroundColor": "white",
                                                    },
                                                    children=[
                                                        html.Img(
                                                            src="assets/images/starship-model.png",
                                                            width="85%",
                                                        )
                                                    ],
                                                    className="box-body",
                                                ),
                                            ]
                                        ),
                                    ],
                                    color="primary",
                                    inverse=True,
                                ),
                                html.Br(),
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            html.H3(
                                                "Research Papers Characterizing Starships"
                                            )
                                        ),
                                        dbc.CardBody(
                                            [
                                                html.Div(
                                                    [
                                                        dash_table.DataTable(
                                                            data=df.to_dict("records"),
                                                            columns=[
                                                                {
                                                                    "name": i,
                                                                    "id": i,
                                                                    "deletable": False,
                                                                    "selectable": False,
                                                                }
                                                                for i in df.columns
                                                            ],
                                                            id="ship_blast_table",
                                                            style_data={
                                                                "whiteSpace": "normal",
                                                                "height": "auto",
                                                            },
                                                            editable=False,
                                                            # filter_action="native",
                                                            # sort_action="native",
                                                            # sort_mode="multi",
                                                            # column_selectable="single",
                                                            # row_selectable="multi",
                                                            row_deletable=False,
                                                            selected_columns=[],
                                                            selected_rows=[],
                                                            page_action="native",
                                                            page_current=0,
                                                            page_size=10,
                                                        ),
                                                    ]
                                                )
                                            ]
                                        ),
                                    ]
                                ),
                            ]
                        ),
                    ],
                    className="box",
                    style={"width": "75%"},
                )
            ],
            className="row",
        )
    ]
)

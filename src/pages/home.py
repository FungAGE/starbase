import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate
from src.components.callbacks import download_ships_button
from src.components.tables import make_paper_table, make_ship_table
import sqlite3
import pandas as pd

dash.register_page(__name__, title="Home", name="Home", path="/")

working = {
    "wiki": "Catalogue/Wiki of Starship Metadata",
    "submit": "Submission of new Starship sequences",
    "blast": "BLAST/HMMER searches",
}
working_buttons = [
    dbc.Button(
        value,
        href=f"/{key}",
        external_link=False,
        color="primary",
        class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
    )
    for key, value in working.items()
]

not_working = [
    html.Div(
        ["Synteny/Genome Browser"],
    ),
    html.Div(
        [
            html.Span(
                "starfish",
                className="logo-text",
            ),
            " webserver",
        ],
    ),
]
not_working_ul = html.Ul(
    [
        html.Li(
            item,
            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
        )
        for item in not_working
    ],
)

columns = [
    "starshipID",
    "familyName",
    "genus",
    "species",
]

modal = dbc.Modal(
    [
        dbc.ModalHeader(dbc.ModalTitle("Download Starships")),
        dbc.ModalBody(
            html.Div(
                [
                    html.Div(id="download-table"),
                    dbc.Button("Download FASTA", id="download-btn-table"),
                    dcc.Download(id="download-fasta"),
                ]
            )
        ),
        dbc.ModalFooter(
            dbc.Button(
                "Close",
                id="close",
                className="ms-auto",
                n_clicks=0,
            )
        ),
    ],
    id="download-modal",
    is_open=False,
    size="lg",
)


layout = dmc.Container(
    fluid=True,
    children=[
        dmc.Center(
            children=[
                dcc.Location(id="url", refresh=False),
                dmc.Title(
                    [
                        html.Span(
                            "starbase: ",
                            className="logo-text",
                        ),
                        "A database and toolkit for exploring large eukaryotic transposable elements in Fungi",
                    ],
                    className="text-center",
                    style={"paddingTop": "20px"},
                    # className="text-center text-custom text-custom-xl",
                ),
            ],
        ),
        dmc.Grid(
            justify="center",
            align="center",
            style={"paddingTop": "20px"},
            gutter="xl",
            children=[
                dmc.GridCol(
                    span="content",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    [
                                        html.Div(
                                            [
                                                "What can I currently use ",
                                                html.Span(
                                                    "starbase",
                                                    className="logo-text",
                                                ),
                                                " for?",
                                            ],
                                            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                        )
                                    ],
                                    style={
                                        "justify-content": "center",
                                    },
                                    className="card-header-custom",
                                ),
                                dbc.CardBody(
                                    [
                                        dbc.Stack(
                                            working_buttons,
                                            # direction="horizontal",
                                            gap=3,
                                            className="justify-content-center",
                                        )
                                    ],
                                    className="d-flex align-items-center",
                                ),
                            ],
                            className="w-100 mb-3",
                        ),
                    ],
                ),
                dmc.GridCol(
                    span="content",
                    children=[
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    [
                                        html.Div(
                                            [
                                                "Functions of ",
                                                html.Span(
                                                    "starbase",
                                                    className="logo-text",
                                                ),
                                                " under active development:",
                                            ],
                                            className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                        )
                                    ],
                                    style={
                                        "justify-content": "center",
                                    },
                                    className="card-header-custom",
                                ),
                                dbc.CardBody(
                                    [not_working_ul],
                                    className="d-flex align-items-center",
                                ),
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    span=12,
                    children=[
                        dmc.Center(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [
                                                html.Div(
                                                    "Data Availability",
                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                )
                                            ],
                                            className="card-header-custom",
                                        ),
                                        dbc.CardBody(
                                            [
                                                dmc.Text(
                                                    [
                                                        "We have been maintaining ",
                                                        html.Span(
                                                            "starbase",
                                                            className="logo-text",
                                                        ),
                                                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~100Mb)",
                                                    ],
                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                ),
                                                download_ships_button,
                                            ],
                                        ),
                                    ],
                                    className="auto-resize-750",
                                    # style={"height": "350px"},
                                    #
                                ),
                                modal,
                            ]
                        ),
                    ],
                ),
            ],
        ),
        dmc.Grid(
            justify="center",
            align="center",
            style={"paddingTop": "20px"},
            grow=True,
            children=[
                dmc.GridCol(
                    span=12,
                    children=[
                        dmc.Center(
                            [
                                dbc.Card(
                                    [
                                        dbc.CardHeader(
                                            [
                                                html.Div(
                                                    ["What is a Starship?"],
                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                                                ),
                                            ],
                                            className="card-header-custom",
                                        ),
                                        dbc.CardBody(
                                            [
                                                dmc.Grid(
                                                    [
                                                        dmc.GridCol(
                                                            span="content",
                                                            children=[
                                                                html.Div(
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
                                                                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl align-items-center",
                                                                    style={
                                                                        "justify-content": "center"
                                                                    },
                                                                ),
                                                            ],
                                                        ),
                                                        dmc.GridCol(
                                                            span="content",
                                                            children=[
                                                                dmc.Image(
                                                                    src="assets/images/starship-model.png",
                                                                    style={
                                                                        "backgroundColor": "white",
                                                                        "maxWidth": "1000px",
                                                                    },
                                                                )
                                                            ],
                                                        ),
                                                    ]
                                                ),
                                            ],
                                        ),
                                    ],
                                    color="primary",
                                    inverse=True,
                                    className="auto-resize-750",
                                ),
                            ]
                        ),
                    ],
                ),
                dmc.GridCol(
                    span="content",
                    children=[
                        dmc.Center(
                            [
                                html.Div(id="paper-table"),
                            ]
                        ),
                    ],
                ),
            ],
        ),
    ],
)


# Callback to handle modal opening and table creation
@callback(
    [
        Output("download-modal", "is_open"),
        Output("download-table", "children"),
    ],
    [
        Input("open-modal", "n_clicks"),
        Input("joined-ships", "data"),
        Input("close", "n_clicks"),
    ],
    State("download-modal", "is_open"),
)
def create_modal_table(dl_click, cached_data, close_modal, is_open):
    # Toggle modal state when the download button is clicked
    if dl_click and not close_modal:
        modal = not is_open  # Toggle state based on current modal state
    elif close_modal:
        modal = False  # Close modal
    else:
        raise PreventUpdate  # No button was clicked

    # Create the table when modal is opening
    if modal:  # Modal is about to be opened
        initial_df = pd.DataFrame(cached_data)
        table = make_ship_table(df=initial_df, id="download-table", columns=columns)
        return modal, table

    # If modal is closing, clear the table
    return modal, None


# Callback to handle FASTA download
@callback(
    Output("download-fasta", "data"),
    Input("download-btn-table", "n_clicks"),
    [
        State("download-table", "derived_virtual_data"),
        State("download-table", "derived_virtual_selected_rows"),
    ],
    prevent_initial_call=True,
)
def download_fasta(n_clicks, rows, selected_rows):
    if n_clicks:
        # Select all rows if no rows are selected
        if not selected_rows:
            selected_rows = list(range(len(rows)))

        ship_names = [
            rows[selected_row]["starshipID"] for selected_row in selected_rows
        ]

        # Query the database with selected ship names
        try:
            conn = sqlite3.connect("database_folder/starbase.sqlite")
            query = f"SELECT genome_name, genome_sequence FROM genome_genome WHERE genome_name "
            if len(ship_names) > 1:
                placeholders = ",".join(["?"] * len(ship_names))
                query += f"IN ({placeholders})"
            else:
                query += "= ?"

            # Execute the query
            df = pd.read_sql_query(query, conn, params=ship_names)

            # Create FASTA content
            fasta_content = [
                f">{row['genome_name']}\n{row['genome_sequence']}"
                for _, row in df.iterrows()
            ]
            fasta_str = "\n".join(fasta_content)

            # Send the FASTA file for download
            return dcc.send_string(fasta_str, filename="starships.fasta")

        except sqlite3.Error as error:
            print("Failed to retrieve data from SQLite table:", error)
            return None

        finally:
            if conn:
                conn.close()


@callback(
    Output("paper-table", "children"),
    [Input("paper-cache", "data"), Input("url", "href")],
)
def load_paper_table(data, url):
    if url:
        table = make_paper_table(data)
        return table

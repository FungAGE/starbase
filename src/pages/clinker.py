import warnings

warnings.filterwarnings("ignore")

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import dcc, html, callback
from dash.dependencies import Output, Input, State

import subprocess
import tempfile
import os
import tempfile
import pandas as pd

from src.components.sqlite import engine
from src.components.tables import make_ship_table
from src.utils.gff2genbank import gff2genbank

dash.register_page(__name__)

layout = dmc.Container(
    fluid=True,
    children=[
        dcc.Location(id="url", refresh=False),
        dmc.Grid(
            justify="center",
            align="top",
            children=[
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading",
                            type="circle",
                            children=[
                                html.Div(
                                    id="clinker-table", className="center-content"
                                ),
                                dbc.Button(
                                    "Show Starship(s) in Viewer",
                                    id="update-button",
                                    n_clicks=0,
                                    className="d-grid gap-2 mx-auto",
                                    style={"fontSize": "1rem"},
                                ),
                            ],
                        ),
                    ],
                ),
                dmc.GridCol(
                    style={"padding": "40px"},
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="loading-1",
                            type="circle",
                            children=[
                                html.Div(
                                    id="clinker-figure",
                                    style={
                                        "height": "800px",
                                        "width": "100%",
                                        "overflow": "auto",
                                        "border": "1px solid #ccc",
                                    },
                                ),
                            ],
                        )
                    ],
                ),
            ],
        ),
    ],
)


def write_tmp(df, seqid, tmp, type=None):
    if type == "gff":
        df.to_csv(tmp, sep="\t", header=False, index=False)
    elif type == "fa":
        seq = df["ship_sequence"][0]
        with open(tmp, "w") as f:
            f.write(f">{seqid}\n{seq}\n")


def fetch_gff(accession):
    query = """
    SELECT g.*
    FROM gff g
    JOIN accessions a ON g.ship_id = a.id
    JOIN ships s on s.accession = a.id
    JOIN joined_ships j ON j.ship_id = a.id
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    WHERE a.accession_tag = :accession_tag AND j.orphan IS NULL
    """

    print("Executing GFF query with accession:", accession)
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})

    if df.empty:
        print(f"No GFF records found for accession_tag: {accession}")

    return df


def fetch_fa(accession):
    query = """
    SELECT s.*
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
    WHERE a.accession_tag = :accession_tag
    """

    # print("Executing FA query with accession:", accession)
    df = pd.read_sql_query(query, engine, params={"accession_tag": accession})

    if df.empty:
        print(f"No FA records found for accession_tag: {accession}")

    return df


def multi_clinker(gb_files, tmp_file):
    subprocess.run(
        [
            "clinker",
            *gb_files,
            "-p",
            tmp_file,
            "-i",
            "0.3",
            "-j",
            "2",
        ],
        check=True,
    )


@callback(
    Output("clinker-figure", "children"),
    Input("update-button", "n_clicks"),
    [
        State("clinker-table", "derived_virtual_selected_rows"),
        State("clinker-table", "derived_virtual_data"),
    ],
)
def update_clinker(n_clicks, selected_rows, table_data):
    if n_clicks > 0:
        tmp_clinker = tempfile.NamedTemporaryFile(suffix=".html", delete=False).name

        if table_data and selected_rows is not None:
            try:
                table_df = pd.DataFrame(table_data)
                # Debug: Check columns and selected_rows
                # print("Columns in table_df:", table_df.columns.tolist())
                # print("Selected rows:", selected_rows)

                if isinstance(selected_rows, list) and all(
                    isinstance(idx, int) for idx in selected_rows
                ):
                    try:
                        rows = table_df.iloc[selected_rows]
                    except IndexError as e:
                        return html.Div("Index out of bounds")

                    tmp_gbs = []
                    for index, row in rows.iterrows():
                        accession = row["accession_tag"]

                        gff_df = fetch_gff(accession)
                        tmp_gff = tempfile.NamedTemporaryFile(
                            suffix=f".gff", delete=False
                        ).name
                        write_tmp(gff_df, accession, tmp_gff, "gff")

                        fa_df = fetch_fa(accession)
                        tmp_fa = os.path.splitext(tmp_gff)[0] + ".fa"
                        write_tmp(fa_df, accession, tmp_fa, "fa")

                        tmp_gb = gff2genbank(gff_file=tmp_gff, fasta_file=tmp_fa)
                        tmp_gbs.append(str(tmp_gb))

                        output = html.P("Select up to four Starships to compare.")
                    if len(selected_rows) > 1 and len(selected_rows) <= 4:

                        multi_clinker(tmp_gbs, tmp_clinker)
                        try:
                            with open(tmp_clinker, "r") as file:
                                pgv_content = file.read()
                        except IOError:
                            output = html.P("Failed to read the temporary file.")

                        output = html.Iframe(
                            srcDoc=pgv_content,
                            style={
                                "width": "100%",
                                "height": "100%",
                                "border": "none",
                            },
                        )
                    elif len(selected_rows) == 1:
                        output = html.P("Select more Starships to compare.")
                    else:
                        output = html.P("No valid selection.")
                else:
                    output = html.H4("Invalid row selection.")
            except Exception as e:
                print("Exception:", e)
                output = html.H4("Error processing data.")
        else:
            output = html.H4("Select Starship(s) to visualize.")
    else:
        output = html.H4(
            "Select the Starship(s) in the table above and click the button to visualize."
        )
    return html.Div(
        [output],
        className="center-content text-center",
    )


@callback(
    Output("clinker-table", "children"),
    Input("url", "href"),
)
def load_ship_table(href):
    query = """
    SELECT a.accession_tag, f.familyName, t.species
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN accessions a ON j.ship_id = a.id
    JOIN ships s on s.accession = a.id
    JOIN gff g ON a.id = g.ship_id
    WHERE s.ship_sequence is NOT NULL AND g.ship_id is NOT NULL AND j.orphan IS NULL
    """
    table_df = pd.read_sql_query(query, engine)
    table_df = table_df.drop_duplicates(subset=["accession_tag"])

    if href:
        table = make_ship_table(df=table_df, id="clinker-table", columns=None, pg_sz=15)
        return table

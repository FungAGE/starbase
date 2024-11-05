from dash import html
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

import os
import logging

from src.components.cache_manager import load_from_cache

logger = logging.getLogger(__name__)

download_ships_button = dbc.Button(
    html.Div(
        [
            "Download Starships from the latest version of ",
            html.Span(
                "starbase",
                className="logo-text",
            ),
            ".",
        ],
    ),
    href="/download",
    color="primary",
    class_name="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl mx-auto",
)

download_ships_card = dbc.Card(
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
                        " data on our GitHub repo (currently private). We are currently in the process of migrating to a new back-end, which will provide more options for data export",
                    ],
                    className="text-custom text-custom-sm text-custom-md text-custom-lg text-custom-xl",
                    style={"paddingBottom": "20px"},
                ),
                dmc.Center(download_ships_button),
            ],
        ),
    ],
    className="auto-resize-750",
    # style={"height": "350px"},
    #
)

download_starbase_button = html.Div(
    [
        dbc.Button(
            html.P(
                [
                    "Download Starships from the latest version of ",
                    html.Span(
                        "starbase",
                        className="logo-text",
                    ),
                    ".",
                ]
            ),
            id="dl-button",
            color="primary",
            className="mt-2",
        ),
    ],
    className="text-center",
    style={
        "font-size": "0.875rem",
    },
)


def curated_switch(text, size="normal"):
    if size == "normal":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "justify-content": "center",
        }
    if size == "large":
        style = {
            "display": "flex",
            "alignItems": "baseline",
            "transform": "scale(1.5)",
            "justify-content": "center",
        }
    switch = dbc.Row(
        justify="center",
        align="start",
        style={"paddingTop": "20px"},
        children=[
            dbc.Col(
                lg=6,
                sm=8,
                children=[
                    dbc.Switch(
                        id="curated-input",
                        label=text,
                        value=False,
                        style=style,
                    ),
                ],
            )
        ],
    )
    return switch


def create_accession_modal(accession):
    # import tempfile
    # from src.pages.pgv import load_fa, write_tmp, load_gff, single_pgv

    initial_df = load_from_cache("meta_data")
    modal_data = initial_df[initial_df["accession_tag"] == accession]

    modal_title = html.H2(f"Ship Accession: {accession}")

    # tmp_pgv = tempfile.NamedTemporaryFile(suffix=".html", delete=True).name

    # with tempfile.TemporaryDirectory() as temp_dir:
    #     tmp_gffs = []
    #     tmp_fas = []
    #     for index, row in modal_data.iterrows():
    #         logger.info(f"Fetching FA for accession: {accession}")
    #         fa_df = load_fa(accession)
    #         tmp_fa = write_tmp(fa_df, accession, "fa", temp_dir)
    #         tmp_fas.append(str(tmp_fa))

    #         logger.info(f"Fetching GFF for accession: {accession}")
    #         gff_df = load_gff(accession)

    #         tmp_gff = write_tmp(gff_df, accession, "gff", temp_dir)
    #         tmp_gffs.append(tmp_gff)

    #         output = html.P("Select up to four Starships to compare.")
    #     single_pgv(tmp_gffs[0], tmp_pgv)
    #     try:
    #         with open(tmp_pgv, "r") as file:
    #             pgv_content = file.read()
    #     except IOError:
    #         output = html.P("Failed to read the temporary file.")

    #     output = html.Iframe(
    #         srcDoc=pgv_content,
    #         style={
    #             "width": "100%",
    #             "height": "100%",
    #             "border": "none",
    #         },
    #     )

    # modal_content = output

    modal_content = html.Div(
        [
            html.Div(
                [html.Strong("ship_id: "), html.Span(modal_data["ship_id"].iloc[0])]
            ),
            html.Div(
                [
                    html.Strong("curated_status: "),
                    html.Span(modal_data["curated_status"].iloc[0]),
                ]
            ),
            html.Div(
                [
                    html.Strong("Number of genomes with ship present:"),
                    html.Span(len(modal_data)),
                ]
            ),
            html.Hr(),
            # TODO: update genome accessions field
            html.Div(
                [
                    html.Strong("Assembly Accession(s): "),
                    html.Span(""),
                ]
            ),
            html.Div(
                [
                    html.Strong("Genome Source: "),
                    html.Span(modal_data["genomeSource"].iloc[0]),
                ]
            ),
            html.Div(
                [
                    html.Strong("Citation: "),
                    html.Span(modal_data["citation"].iloc[0]),
                ]
            ),
            html.Div(
                [html.Strong("contigID: "), html.Span(modal_data["contigID"].iloc[0])]
            ),
            html.Div(
                [
                    html.Strong("elementBegin: "),
                    html.Span(modal_data["elementBegin"].iloc[0]),
                ]
            ),
            html.Div(
                [
                    html.Strong("elementEnd: "),
                    html.Span(modal_data["elementEnd"].iloc[0]),
                ]
            ),
            html.Div([html.Strong("size: "), html.Span(modal_data["size"].iloc[0])]),
            html.Hr(),
            html.Div([html.Strong("order: "), html.Span(modal_data["order"].iloc[0])]),
            html.Div(
                [html.Strong("family: "), html.Span(modal_data["family"].iloc[0])]
            ),
            html.Div(
                [html.Strong("species: "), html.Span(modal_data["species"].iloc[0])]
            ),
            html.Div(
                [html.Strong("strain: "), html.Span(modal_data["strain"].iloc[0])]
            ),
            html.Div(
                [
                    html.Strong("taxid: "),
                    html.A(
                        (
                            int(modal_data["taxid"].iloc[0])
                            if isinstance(modal_data["taxid"].iloc[0], (float, int))
                            else modal_data["taxid"].iloc[0]
                        ),
                        href=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={int(modal_data['taxid'].iloc[0]) if isinstance(modal_data['taxid'].iloc[0], (float, int)) else modal_data['taxid'].iloc[0]}",
                    ),
                ]
            ),
        ]
    )

    return modal_content, modal_title

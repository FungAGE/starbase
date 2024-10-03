import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.DEBUG)


import subprocess
import tempfile
import os
import glob
import base64

import dash
import dash_mantine_components as dmc

from dash import dcc, html, callback
from dash.dependencies import Output, Input

from src.utils.blast_utils import check_input, guess_seq_type, write_temp_fasta
from src.utils.tree import plot_tree
from src.utils.parsing import parse_fasta
from src.components.callbacks import MOUNTED_DIRECTORY_PATH

dash.register_page(__name__)


def add_to_tree(input_seq, fasta, tree):
    if input_seq:
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        mafft_cmd = f"mafft --thread 2 --addfragments {input_seq} {fasta} > {tmp_fasta}"
        subprocess.run(
            mafft_cmd,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # temp_dir = tempfile.TemporaryDirectory().name
        # iqtree_cmd = f"iqtree -T 2 -s {tmp_fasta} -g {tree} --prefix {temp_dir}"
        # subprocess.run(iqtree_cmd, shell=True)

        # search_pattern = os.path.join(temp_dir, "*.treefile")

        # # Use glob to find all files matching the pattern
        # matching_files = glob.glob(search_pattern)

        # # Return the first matching file or None if no match is found
        # if matching_files:
        #     return matching_files[0]  # or return the full list if you need all matches
        # else:
        #     return None
        tmp_tree = tempfile.NamedTemporaryFile(suffix=".treefile", delete=False).name
        fasttree_cmd = f"fasttree {tmp_fasta} > {tmp_tree}"
        subprocess.run(
            fasttree_cmd,
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return tmp_tree


layout = dmc.Container(
    fluid=True,
    children=[
        dmc.Grid(
            justify="start",
            align="start",
            grow=True,
            style={"paddingTop": "20px"},
            gutter="xl",
            children=[
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Upload(
                            id="fasta-upload",
                            children=html.Div(
                                "Drag and drop or click to select a FASTA file.",
                                id="fasta-sequence-upload",
                                style={"fontSize": "1rem"},
                            ),
                            className="upload-box",
                            multiple=False,
                            accept=".fa, .fas, .fasta, .fna",
                        ),
                    ],
                ),
                dmc.GridCol(
                    span={
                        "sm": 12,
                        "lg": 8,
                    },
                    children=[
                        dcc.Loading(
                            id="tree-loading",
                            type="default",
                            children=html.Div(
                                id="phylogeny",
                            ),
                        ),
                    ],
                ),
            ],
        )
    ],
)


@callback(
    Output("fasta-sequence-upload", "children"),
    [
        Input("fasta-upload", "contents"),
        Input("fasta-upload", "filename"),
    ],
)
def update_fasta_details(seq_content, seq_filename):
    if seq_content is None:
        return [
            html.Div(
                html.P(
                    ["Select a FASTA file to upload"],
                )
            )
        ]
    else:
        try:
            # "," is the delimeter for splitting content_type from content_string
            content_type, content_string = seq_content.split(",")
            query_string = base64.b64decode(content_string).decode("utf-8")
            children = parse_fasta(query_string, seq_filename)
            return children

        except Exception as e:
            logging.error(e)
            return html.Div(["There was an error processing this file."])


@callback(
    Output("phylogeny", "children"),
    Input("fasta-upload", "contents"),
)
def update_ui(fasta_upload):
    input_type, input_header, input_seq = check_input(None, fasta_upload)
    query_type = guess_seq_type(input_seq)
    tmp_query_fasta = write_temp_fasta(input_header, input_seq)
    if query_type == "prot":
        new_tree = add_to_tree(
            tmp_query_fasta,
            f"{MOUNTED_DIRECTORY_PATH}/Starships/captain/tyr/faa/alignments/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit",
            "src/data/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile",
        )
        if new_tree is not None:
            output = dcc.Graph(
                id="phylogeny",
                className="div-card",
                figure=plot_tree(new_tree, highlight_families=input_header),
            )
            return output
    else:
        return None

import warnings

warnings.filterwarnings("ignore")

import logging

logging.basicConfig(level=logging.ERROR)

import os
import glob
import subprocess
import tempfile
import json
import base64

import dash
import dash_mantine_components as dmc

from dash import dcc, html, callback, callback_context
from dash.dependencies import Output, Input

from Bio import SeqIO

from src.config.settings import PHYLOGENY_PATHS
from src.utils.blast_utils import check_input, guess_seq_type, write_temp_fasta
from src.utils.tree import plot_tree
from src.utils.seq_utils import parse_fasta, load_fasta_to_dict

dash.register_page(__name__)


def run_mafft(query, ref_msa):
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    MODEL = "PROTGTR+G+F"
    fasta_dict = load_fasta_to_dict(query)
    tmp_headers = list(fasta_dict.keys())

    mafft_cmd = f"mafft --thread 2 --auto --addfragments {query} --keeplength {ref_msa} | seqkit grep -n -f {tmp_headers} > {tmp_fasta}"
    subprocess.run(
        mafft_cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return tmp_headers, tmp_fasta


def add_to_tree(query_msa, tree, ref_msa, model, tmp_dir):
    if query_msa:
        # Construct the epa-ng command
        epa_cmd = f"epa-ng -T 2 --redo --ref-msa {ref_msa} --tree {tree} --query {query_msa} --model {model} --out-dir {tmp_dir}"

        # Run the command
        result = subprocess.run(
            epa_cmd,
            shell=True,
            stdout=subprocess.DEVNULL,  # Uncomment to suppress stdout
            stderr=subprocess.DEVNULL,  # Uncomment to suppress stderr
        )

        # Check for errors in command execution
        if result.returncode != 0:
            raise RuntimeError(
                f"epa-ng command failed with return code {result.returncode}"
            )

        tmp_tree = os.path.join(tmp_dir, "epa_result.jplace")
        return tmp_tree


def gappa(tmp_tree):
    out_dir = os.path.dirname(tmp_tree)
    out_tree = os.path.join(out_dir, "epa_result.newick")
    gappa_cmd = (
        f"gappa examine graft --out-dir {out_dir} --threads 2 --jplace-path {tmp_tree}"
    )
    subprocess.run(
        gappa_cmd,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    logging.info(f"Tree output: {out_tree}")

    return out_tree


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
                    span="content",
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
    output = None
    headers_list = None
    if fasta_upload:
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_type, inputs = check_input(None, fasta_upload)
            query_type = guess_seq_type(inputs)
            tmp_query_fasta = write_temp_fasta(inputs)
            fasta_dict = load_fasta_to_dict(tmp_query_fasta)
            headers_list = list(fasta_dict.keys())
            
            if query_type == "prot":
                tmp_headers, query_msa = run_mafft(
                    query=tmp_query_fasta,
                    ref_msa=PHYLOGENY_PATHS["msa"],
                )
                new_tree = add_to_tree(
                    query_msa=query_msa,
                    tree=PHYLOGENY_PATHS["tree"],
                    ref_msa=PHYLOGENY_PATHS["msa"],
                    model="PROTGTR+G+F",
                    tmp_dir=tmp_dir,
                )
                out_tree = gappa(new_tree)
                with open(out_tree, "r") as file:
                    file_contents = file.read()
                logging.info(file_contents)

                if out_tree is not None:
                    output = dcc.Graph(
                        id="phylogeny",
                        className="div-card",
                        figure=plot_tree(
                            out_tree, highlight_families="all", tips=headers_list
                        ),
                    )
        return output

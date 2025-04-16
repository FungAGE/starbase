# dash page for classification workflow

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_core_components as dcc

from dash import html, callback
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

import os
import logging
import tempfile
import pandas as pd

from src.utils.seq_utils import guess_seq_type
from src.config.settings import BLAST_DB_PATHS
from src.utils.classification_utils import create_classification_callback
from src.tasks import check_contained_matches_task, check_similar_matches_task, check_exact_matches_task, run_family_classification_task, run_navis_classification_task, run_haplotype_classification_task
from src.components.callbacks import create_file_upload
from src.components.error_boundary import handle_callback_error
from src.utils.seq_utils import parse_fasta_from_file, write_fasta
from src.database.sql_manager import fetch_captains, fetch_ships


dash.register_page(__name__)

logger = logging.getLogger(__name__)

def make_progress_bar(message, id, color):
    return dmc.Group([
        dmc.Text(message, size="sm", w=100),
        dbc.Progress(
            id=id,
            value=0,
            color=color,
            animated=False,
            striped=False,
            style={"width": "100%", "marginBottom": "10px"}
        ),
        html.Div(id=id)
    ])


layout = dmc.Container(
    size="md",
    children=[
        html.H1("Classification Workflow"),

        # Stores
        # file upload
        dcc.Store(id="classification-upload"),

        # text for stage
        dcc.Store(id="classification-stage"),

        # matches store
        dcc.Store(id="check-exact-matches"),
        dcc.Store(id="check-exact-matches-progress"),
        dcc.Store(id="check-exact-matches-color"),
        dcc.Store(id="check-contained-matches"),
        dcc.Store(id="check-contained-matches-progress"),
        dcc.Store(id="check-contained-matches-color"),
        dcc.Store(id="check-similar-matches"),
        dcc.Store(id="check-similar-matches-progress"),
        dcc.Store(id="check-similar-matches"),
        dcc.Store(id="check-similar-matches"),

        # denovo annotation store
        dcc.Store(id="classification-denovo-annotation-progress"),
        dcc.Store(id="classification-denovo-annotation-color"),
        dcc.Store(id="classification-denovo-annotation-name"),
        dcc.Store(id="classification-denovo-annotation-output"),
        dcc.Store(id="classification-denovo-annotation-data"),

        # annotate store
        dcc.Store(id="classification-annotate-progress"),
        dcc.Store(id="classification-annotate-color"),
        dcc.Store(id="classification-codon-fasta"),
        dcc.Store(id="classification-fasta"),
        dcc.Store(id="classification-gff"),

        # progress for family classification
        dcc.Store(id="classification-family-progress"),
        dcc.Store(id="classification-family-color"),
        dcc.Store(id="classification-family-name"),
        dcc.Store(id="classification-family-output"),
        dcc.Store(id="classification-family-data"),

        # progress for navis classification
        dcc.Store(id="classification-navis-progress"),
        dcc.Store(id="classification-navis-color"),
        dcc.Store(id="classification-navis-name"),
        dcc.Store(id="classification-navis-output"),
        dcc.Store(id="classification-navis-data"),

        # progress for haplotype classification
        dcc.Store(id="classification-haplotype-progress"),
        dcc.Store(id="classification-haplotype-color"),
        dcc.Store(id="classification-haplotype-name"),
        dcc.Store(id="classification-haplotype-output"),
        dcc.Store(id="classification-haplotype-data"),

        dmc.Stack([
            dmc.Paper(
                children=create_file_upload(
                    upload_id="classification-fasta-upload",
                    output_id="classification-fasta-sequence-upload",
                    accept_types=[".fa", ".fas", ".fasta", ".fna"],
                    placeholder_text="Drag and drop or click to select a FASTA file"
                ),
                withBorder=False,
                radius="md",
                style={"cursor": "pointer"}
            ),
            dmc.Stack([
                dmc.Group([
                    dmc.Text("Current Stage:", size="lg", fw=500),
                    dmc.Text(id="classification-stage", size="lg", c="blue")
                ]),
                dmc.Stack([
                    make_progress_bar("Checks for exact matches", "check-exact-matches-progress", "green"),
                    make_progress_bar("Checks for contained matches", "check-contained-matches-progress", "blue"),
                    make_progress_bar("Checks for similar matches", "check-similar-matches-progress", "violet"),
                    make_progress_bar("Running denovo annotation", "run-denovo-annotation-progress", "orange"),
                    make_progress_bar("Running family classification", "run-family-classification-progress", "green"),
                    make_progress_bar("Running navis classification", "run-navis-classification-progress", "blue"),
                    make_progress_bar("Running haplotype classification", "run-haplotype-classification-progress", "violet"),
                ]),
            ]),
        ], gap="md"),
        html.Div(id='classification-output', className="mt-4"),
    ]
)

# Callbacks
# handle upload and start family classification progress
@callback(
    [
        Output("classification-upload", "data"),
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-family-progress", "value", allow_duplicate=True),
        Output("classification-family-progress", "animated", allow_duplicate=True),
        Output("classification-family-progress", "striped", allow_duplicate=True),
        Output("classification-family-color", "color", allow_duplicate=True),
    ],
    Input("classification-fasta-upload", "contents"),
    prevent_initial_call=True
)
@handle_callback_error
def setup_classification(seq_content):
    if seq_content is None:
        raise PreventUpdate
        
    # Parse and validate FASTA
    _, seq, fasta_error = parse_fasta_from_file(seq_content)
    if fasta_error or not seq:
        raise ValueError(f"FASTA parsing error: {fasta_error or 'No sequence found'}")
        
    seq_type = guess_seq_type(seq)
    if not seq_type:
        raise ValueError("Could not determine sequence type")

    # Save sequence to temporary file
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_fasta({"query_sequence": seq}, tmp_fasta)
    
    return (
        {
            "seq_type": seq_type, 
            "fasta": tmp_fasta, 
            "fetch_ship_params": {
                "curated": False,
                "with_sequence": True,
                "dereplicate": True
            },
            "fetch_captain_params": {
                "curated": True,
                "with_sequence": True
            }

        },
        "Checking for exact matches",
        10,
        True,
        True,
        "green",
    )

def make_alert_using_match(match_class=None, classification=None, type=None):    
    if match_class:
        if classification == "Family":
            color = "green"
        elif classification == "Navis":
            color = "blue"
        elif classification == "Haplotype":
            color = "violet"

        return (
            dmc.Alert(
                title=f"{classification} found using {type} match",
                children=f"{classification}: {match_class}",
                color=color,
                variant="light",
                className="mb-3"
            )
        )
    else:
        return None

create_classification_callback(
    task_name="Exact",
    task_function=check_exact_matches_task,
    input_store="classification-upload",
    output_store="classification-exact-matches",
    next_stage="Checking for Contained Matches",
    active_progress="family",
    custom_task_args=lambda upload_data, input_data: {
        "sequence": upload_data["fasta"],
        "ships_dict": fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
    }
)

create_classification_callback(
    task_name="Contained",
    task_function=check_contained_matches_task,
    input_store="classification-exact-matches",
    output_store="classification-contained-matches",
    next_stage="Checking for Similar Matches",
    active_progress="family",
    custom_task_args=lambda upload_data, input_data: {
        "sequence": upload_data["fasta"],
        "ships_dict": fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
    }
)

create_classification_callback(
    task_name="Similar",
    task_function=check_similar_matches_task,
    input_store="classification-contained-matches",
    output_store="classification-similar-matches",
    next_stage="Running Family Classification",
    active_progress="family",
    custom_task_args=lambda upload_data, input_data: {
        "sequence": upload_data["fasta"],
        "ships_dict": fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
    }
)

# Usage example for family classification:
create_classification_callback(
    task_name="Family",
    task_function=run_family_classification_task,
    input_store="classification-similar-matches",
    output_store="classification-family-data",
    next_stage="Running Navis Classification",
    active_progress="family",
    progress_color="green",
    next_progress="navis",
    next_progress_color="blue",
    custom_task_args=lambda upload_data, input_data: {
        "fasta": upload_data["fasta"],
        "seq_type": upload_data.get("seq_type", "nucleotide"),
        "db_list": BLAST_DB_PATHS
    }
)

# Usage example for navis classification:
create_classification_callback(
    task_name="Navis",
    task_function=run_navis_classification_task,
    input_store="classification-family-data",
    output_store="classification-navis-data",
    next_stage="Running Haplotype Classification",
    active_progress="navis",
    progress_color="blue",
    next_progress="haplotype",
    next_progress_color="violet",
    custom_task_args=lambda upload_data, input_data: {
        "protein": input_data["protein"],
        "existing_captains": fetch_captains(**upload_data["fetch_captain_params"]),
        "threads": 1
    }
)

# Usage example for haplotype classification:
create_classification_callback(
    task_name="Haplotype",
    task_function=run_haplotype_classification_task,
    input_store="classification-navis-data",
    output_store="classification-haplotype-data",
    next_stage="Classification Complete",
    active_progress="haplotype",
    progress_color="violet",
    custom_task_args=lambda upload_data, input_data: {
        "fasta": upload_data["fasta"],
        "existing_ships": fetch_ships(**upload_data["fetch_ship_params"])[
            fetch_ships(**upload_data["fetch_ship_params"])["curated"] == True
        ],
        "navis": input_data["navis"]
    }
)
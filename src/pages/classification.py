# dash page for classification workflow

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_core_components as dcc

from dash import html, callback
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

import logging
import tempfile
import pandas as pd

from src.utils.seq_utils import guess_seq_type
from src.config.settings import BLAST_DB_PATHS
from src.utils.classification_utils import check_exact_match, check_contained_match, check_similar_match, classify_family, classify_navis, classify_haplotype
from src.components.callbacks import create_file_upload
from src.components.error_boundary import handle_callback_error
from src.utils.seq_utils import parse_fasta_from_file, write_fasta
from src.database.sql_manager import fetch_captains, fetch_ships


dash.register_page(__name__)

logger = logging.getLogger(__name__)

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
        dcc.Store(id="classification-exact-matches"),
        dcc.Store(id="classification-contained-matches"),
        dcc.Store(id="classification-similar-matches"),

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
                    dmc.Group([
                        dmc.Text("Family:", size="sm", w=100),
                        dbc.Progress(
                            id="classification-family-progress",
                            value=0,
                            color="green",
                            animated=False,
                            striped=False,
                            style={"width": "100%", "marginBottom": "10px"}
                        ),
                        html.Div(id="classification-family-output")
                    ]),
                    dmc.Group([
                        dmc.Text("Navis:", size="sm", w=100),
                        dbc.Progress(
                            id="classification-navis-progress",
                            value=0,
                            color="blue",
                            animated=False,
                            striped=False,
                            style={"width": "100%", "marginBottom": "10px"}
                        ),
                        html.Div(id="classification-navis-output")
                    ]),
                    dmc.Group([
                        dmc.Text("Haplotype:", size="sm", w=100),
                        dbc.Progress(
                            id="classification-haplotype-progress",
                            value=0,
                            color="violet",
                            animated=False,
                            striped=False,
                            style={"width": "100%", "marginBottom": "10px"}
                        ),
                        html.Div(id="classification-haplotype-output")
                    ]),
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
            "protein": tmp_fasta,
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

# First, check for exact matches
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
        Output("classification-navis-output", "children", allow_duplicate=True),
        Output("classification-haplotype-output", "children", allow_duplicate=True),
        Output("classification-exact-matches", "data"),
    ],
    [Input("classification-upload", "data")],
    prevent_initial_call=True
)
@handle_callback_error
def check_exact_matches(data):
    print("Starting exact match check...")
    if data is None:
        print("No data provided for exact match check")
        raise PreventUpdate
    
    # Fetch ships data using query parameters
    fetch_ship_params = data["fetch_ship_params"]
    existing_ships = fetch_ships(**fetch_ship_params)
    print(f"Checking exact matches against {len(existing_ships)} existing ships")
    
    exact_match = check_exact_match(data["protein"], existing_ships)
    matches_found = False

    if exact_match:
        print(f"Found exact match: {exact_match}")
        match_row = existing_ships[existing_ships["accession_tag"] == exact_match]
        match_family = match_row["familyName"].values[0]
        match_navis = match_row["starship_navis"].values[0]
        match_haplotype = match_row["starship_haplotype"].values[0]
        matches_found = True
        print(f"Exact match details - Family: {match_family}, Navis: {match_navis}, Haplotype: {match_haplotype}")
    else:
        print("No exact matches found")
        match_family = match_navis = match_haplotype = None

    return (
        "Checking for contained matches",
        make_alert_using_match(match_family, "Family", "exact"),
        make_alert_using_match(match_navis, "Navis", "exact"),
        make_alert_using_match(match_haplotype, "Haplotype", "exact"),
        {"found": matches_found} if matches_found else None
    )

# then, check for contained matches
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
        Output("classification-navis-output", "children", allow_duplicate=True),
        Output("classification-haplotype-output", "children", allow_duplicate=True),
        Output("classification-contained-matches", "data"),
    ],
    [Input("classification-exact-matches", "data")],
    [State("classification-upload", "data")],
    prevent_initial_call=True
)
@handle_callback_error
def check_contained_matches(exact_matches, data):
    print("Starting contained match check...")
    if exact_matches or data is None:
        logger.info("Skipping contained match check - exact match was found or no data")
        raise PreventUpdate
    
    # Fetch ships data using query parameters
    fetch_ship_params = data["fetch_ship_params"]
    existing_ships = fetch_ships(**fetch_ship_params)
    print(f"Checking contained matches against {len(existing_ships)} existing ships")
    
    contained_match = check_contained_match(data["protein"], existing_ships)
    matches_found = contained_match is not None

    if matches_found:
        print(f"Found contained match: {contained_match}")
    else:
        print("No contained matches found")

    return (
        "Checking for similar matches",  # stage
        make_alert_using_match(contained_match, "Family", "contained"),  # family output
        make_alert_using_match(contained_match, "Navis", "contained"),  # navis output
        make_alert_using_match(contained_match, "Haplotype", "contained"),  # haplotype output
        {"found": matches_found} if matches_found else None  # contained matches data
    )

# then, check for similar matches
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
        Output("classification-navis-output", "children", allow_duplicate=True),
        Output("classification-haplotype-output", "children", allow_duplicate=True),
        Output("classification-similar-matches", "data"),
    ],
    [Input("classification-contained-matches", "data")],
    [State("classification-upload", "data")],
    prevent_initial_call=True
)
@handle_callback_error
def check_similar_matches(contained_matches, data):
    logging.info("Starting similar match check...")
    if contained_matches or data is None:
        logging.info("Skipping similar match check - contained match was found or no data")
        raise PreventUpdate
    
    # Fetch ships data using query parameters
    fetch_ship_params = data["fetch_ship_params"]
    existing_ships = fetch_ships(**fetch_ship_params)
    logging.info(f"Checking similar matches against {len(existing_ships)} existing ships")
    
    # Add threshold parameter (you may want to adjust this value)
    similar_match = check_similar_match(data["protein"], existing_ships, threshold=0.9)
    matches_found = similar_match is not None

    if matches_found:
        logging.info(f"Found similar match: {similar_match}")
    else:
        logging.warning("No similar matches found")

    # Make sure to return a tuple matching the number of outputs
    return (
        "Running Family Classification",  # stage
        make_alert_using_match(similar_match, "Family", "similar"),  # family output
        make_alert_using_match(similar_match, "Navis", "similar"),  # navis output
        make_alert_using_match(similar_match, "Haplotype", "similar"),  # haplotype output
        {"found": matches_found} if matches_found else None  # similar matches data
    )

# Run Family Classification
# outputs will be return when work is complete
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-family-progress", "value", allow_duplicate=True),
        Output("classification-family-progress", "animated", allow_duplicate=True),
        Output("classification-family-progress", "striped", allow_duplicate=True),
        Output("classification-family-color", "color", allow_duplicate=True),
        Output("classification-family-name", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
        Output("classification-family-data", "data"),
        Output("classification-navis-progress", "value", allow_duplicate=True),
        Output("classification-navis-progress", "animated", allow_duplicate=True),
        Output("classification-navis-progress", "striped", allow_duplicate=True),
        Output("classification-navis-color", "color", allow_duplicate=True),
    ],
    [
        Input("classification-similar-matches", "data"),
    ],
    [
        State("classification-upload", "data"),
        State("classification-exact-matches", "data"),
        State("classification-contained-matches", "data"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def run_family_classification(similar_matches, data, exact_matches, contained_matches):
    # Check for required data
    if data is None:
        raise PreventUpdate
        
    # Check if any previous matches were found
    if (exact_matches and exact_matches.get("found")) or \
       (contained_matches and contained_matches.get("found")) or \
       (similar_matches and similar_matches.get("found")):
        raise PreventUpdate

    try:
        # Run the classification
        family_dict, tmp_protein = classify_family(
            fasta=data["fasta"],
            seq_type=data["seq_type"],
            db_list=BLAST_DB_PATHS,
            threads=1
        )

        if not family_dict:
            raise ValueError("Family classification failed")

        # Success case
        return (
            "Running Navis Classification",  # stage
            100,  # family progress complete
            False,  # family not animated
            True,  # family stays striped
            "green",  # family success color
            family_dict['family'],  # family name
            dmc.Alert(  # output
                title="Family Classification",
                children=[
                    f"Family: {family_dict['family']}",
                    html.Br(),
                    f"Alignment Length: {family_dict['aln_length']}",
                    html.Br(),
                    f"E-value: {family_dict['evalue']}"
                ],
                color="green",
                variant="light",
                className="mb-3"
            ),
            {"protein": tmp_protein},  # data for next step
            10,  # navis progress starts
            True,  # navis animated
            True,  # navis striped
            "blue",  # navis color
        )

    except Exception as e:
        # Failure case
        return (
            "Family Classification Failed",  # stage
            100,  # family progress complete
            False,  # family not animated
            False,  # family not striped
            "red",  # family failure color
            "Error",  # family name
            dmc.Alert(  # error output
                title="Family Classification Error",
                children=str(e),
                color="red",
                variant="light",
                className="mb-3"
            ),
            None,  # no data for next step
            0,  # navis progress reset
            False,  # navis not animated
            False,  # navis not striped
            "gray",  # navis disabled color
        )

  
# Run Navis Classification
# outputs will be return when work is complete
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-navis-progress", "value", allow_duplicate=True),
        Output("classification-navis-progress", "animated", allow_duplicate=True),
        Output("classification-navis-progress", "striped", allow_duplicate=True),
        Output("classification-navis-color", "color", allow_duplicate=True),
        Output("classification-navis-name", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
        Output("classification-navis-data", "data"),
        Output("classification-haplotype-progress", "value", allow_duplicate=True),
        Output("classification-haplotype-progress", "animated", allow_duplicate=True),
        Output("classification-haplotype-progress", "striped", allow_duplicate=True),
        Output("classification-haplotype-color", "color", allow_duplicate=True),
    ],
    [
        Input("classification-family-name", "children"),
    ],
    [
        State("classification-upload", "data"),
        State("classification-family-data", "data"),
        State("classification-exact-matches", "data"),
        State("classification-contained-matches", "data"),
        State("classification-similar-matches", "data"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def run_navis_classification(family_name, data, family_data, exact_matches, contained_matches, similar_matches):
    # Check for required data
    if data is None or family_data is None:
        raise PreventUpdate
        
    # Check if any previous matches were found
    if (exact_matches and exact_matches.get("found")) or \
       (contained_matches and contained_matches.get("found")) or \
       (similar_matches and similar_matches.get("found")):
        raise PreventUpdate
        
    # Check if family classification failed
    if family_name == "Error":
        raise PreventUpdate

    try:
        fetch_captain_params = data["fetch_captain_params"]
        existing_captains = fetch_captains(**fetch_captain_params)
        navis_name = classify_navis(
            fasta=data["protein"],
            existing_captains=existing_captains,
            threads=1
        )

        if not navis_name:
            raise ValueError("Navis classification failed")

        # Success case
        return (
            "Running Haplotype Classification",  # stage
            100,  # navis progress complete
            False,  # navis not animated
            True,  # navis stays striped
            "green",  # navis success color
            navis_name,  # navis name
            dmc.Alert(  # output
                title="Navis Classification",
                children=f"Navis: {navis_name}",
                color="blue",
                variant="light",
                className="mb-3"
            ),
            {"protein": data["protein"], "navis": navis_name},  # data for next step
            10,  # haplotype progress starts
            True,  # haplotype animated
            True,  # haplotype striped
            "violet",  # haplotype color
        )

    except Exception as e:
        # Failure case
        return (
            "Navis Classification Failed",  # stage
            100,  # navis progress complete
            False,  # navis not animated
            False,  # navis not striped
            "red",  # navis failure color
            "Error",  # navis name
            dmc.Alert(  # error output
                title="Navis Classification Error",
                children=str(e),
                color="red",
                variant="light",
                className="mb-3"
            ),
            None,  # no data for next step
            0,  # haplotype progress reset
            False,  # haplotype not animated
            False,  # haplotype not striped
            "gray",  # haplotype disabled color
        )

# Run Haplotype Classification
# outputs will be return when work is complete
@callback(
    [
        Output("classification-stage", "children", allow_duplicate=True),
        Output("classification-haplotype-progress", "value", allow_duplicate=True),
        Output("classification-haplotype-progress", "animated", allow_duplicate=True),
        Output("classification-haplotype-progress", "striped", allow_duplicate=True),
        Output("classification-haplotype-color", "color", allow_duplicate=True),
        Output("classification-haplotype-name", "children", allow_duplicate=True),
        Output("classification-family-output", "children", allow_duplicate=True),
    ],
    [
        Input("classification-navis-data", "data"),
    ],
    [
        State("classification-upload", "data"),
        State("classification-exact-matches", "data"),
        State("classification-contained-matches", "data"),
        State("classification-similar-matches", "data"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def run_haplotype_classification(navis_data, data, exact_matches, contained_matches, similar_matches):
    # Check for required data
    if data is None or navis_data is None:
        raise PreventUpdate
        
    # Check if any previous matches were found
    if (exact_matches and exact_matches.get("found")) or \
       (contained_matches and contained_matches.get("found")) or \
       (similar_matches and similar_matches.get("found")):
        raise PreventUpdate

    try:
        fetch_ship_params = data["fetch_ship_params"]
        existing_ships = fetch_ships(**fetch_ship_params)
        curated_ships = existing_ships[existing_ships["curated"] == True]

        haplotype_name = classify_haplotype(
            fasta=data["protein"],
            existing_ships=curated_ships,
            navis=data["navis"]
        )

        if not haplotype_name:
            raise ValueError("Haplotype classification failed")

        # Success case
        return (
            "Classification Complete",  # stage
            100,  # haplotype progress complete
            False,  # haplotype not animated
            True,  # haplotype stays striped
            "green",  # haplotype success color
            haplotype_name,  # haplotype name
            dmc.Alert(  # output
                title="Haplotype Classification",
                children=f"Haplotype: {haplotype_name}",
                color="violet",
                variant="light",
                className="mb-3"
            ),
        )

    except Exception as e:
        # Failure case
        return (
            "Haplotype Classification Failed",  # stage
            100,  # haplotype progress complete
            False,  # haplotype not animated
            False,  # haplotype not striped
            "red",  # haplotype failure color
            "Error",  # haplotype name
            dmc.Alert(  # error output
                title="Haplotype Classification Error",
                children=str(e),
                color="red",
                variant="light",
                className="mb-3"
            ),
        )
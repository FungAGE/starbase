# dash page for classification workflow

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_core_components as dcc

from dash import html, callback
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

import base64
import logging
import tempfile
import os

from src.utils.seq_utils import guess_seq_type
from src.config.settings import BLAST_DB_PATHS
from src.utils.classification_utils import classify_sequence, classify_family, classify_navis, classify_haplotype
from src.components.callbacks import create_file_upload
from src.components.error_boundary import handle_callback_error
from src.utils.seq_utils import parse_fasta, parse_fasta_from_file, write_fasta
from src.database.sql_manager import fetch_captains, fetch_ships


dash.register_page(__name__)

logger = logging.getLogger(__name__)

layout = dmc.Container(
    size="md",
    children=[
        html.H1("Classification Workflow"),
        # Hidden stores for passing data between callbacks
        dcc.Store(id="classification-data-store"),
        dcc.Store(id="classification-temp-files-store"),
        
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
                dbc.Progress(
                    [
                        dbc.Progress(
                            id="family-progress",
                            value=0,
                            color="green",
                            bar=True,
                            striped=False,
                            animated=False,
                            label="Family",
                            style={"width": "33%"}
                        ),
                        dbc.Progress(
                            id="navis-progress",
                            value=0,
                            color="blue",
                            bar=True,
                            striped=False,
                            animated=False,
                            label="Navis",
                            style={"width": "33%"}
                        ),
                        dbc.Progress(
                            id="haplotype-progress",
                            value=0,
                            color="violet",
                            bar=True,
                            striped=False,
                            animated=False,
                            label="Haplotype",
                            style={"width": "34%"}
                        ),
                    ],
                    style={"width": "100%"},
                    className="mb-3",
                ),
            ]),
        ], gap="md"),
        html.Div(id='classification-output', className="mt-4"),
    ]
)

@callback(
    [
        Output("classification-data-store", "data"),
        Output("classification-temp-files-store", "data"),
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
        {"seq_type": seq_type, "results": [], "stage": "family"},  # data store
        {"fasta": tmp_fasta, "protein": None}  # temp files store
    )

# Step 2: Family Classification
@callback(
    [
        Output("classification-stage", "children"),
        Output("family-progress", "value"),
        Output("family-progress", "color"),
        Output("family-progress", "striped"),
        Output("family-progress", "animated"),
        Output("navis-progress", "value"),
        Output("navis-progress", "color"),
        Output("navis-progress", "striped"),
        Output("navis-progress", "animated"),
        Output("haplotype-progress", "value"),
        Output("haplotype-progress", "color"),
        Output("haplotype-progress", "striped"),
        Output("haplotype-progress", "animated"),
        Output("classification-output", "children"),
    ],
    [
        Input("classification-data-store", "data"),
        Input("classification-temp-files-store", "data"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def update_classification_status(data, temp_files):
    if not data or not temp_files:
        raise PreventUpdate
        
    stage = data.get("stage", "family")
    results = []
    
    # Initialize all progress bars as inactive
    family_progress = navis_progress = haplotype_progress = 0
    family_color = navis_color = haplotype_color = "secondary"
    family_striped = navis_striped = haplotype_striped = False
    family_animated = navis_animated = haplotype_animated = False
    
    # Start with family classification
    family_color = "green"
    family_striped = family_animated = True
    
    # Run family classification
    family_dict, tmp_protein = classify_family(
        fasta=temp_files["fasta"],
        seq_type=data["seq_type"],
        db_list=BLAST_DB_PATHS,
        threads=1
    )
    
    if not family_dict:
        # Family classification failed
        family_color = "secondary"
        return (
            "Classification Failed at Family Stage",
            0, family_color, False, False,
            0, "secondary", False, False,
            0, "secondary", False, False,
            html.Div([])
        )
        
    # Family classification succeeded
    family_progress = 100
    family_striped = family_animated = False
    try:
        evalue = float(family_dict['evalue'])
        evalue_str = f"{evalue:.2e}"
    except (ValueError, TypeError):
        evalue_str = str(family_dict['evalue'])
        
    results.append(dmc.Alert(
        title="Family Classification",
        children=[
            f"Family: {family_dict['family']}",
            html.Br(),
            f"Alignment Length: {family_dict['aln_length']}",
            html.Br(),
            f"E-value: {evalue_str}"
        ],
        color="green",
        variant="light",
        className="mb-3"
    ))
    
    if not tmp_protein:
        return (
            "Classification Stopped after Family Stage",
            100, "green", False, False,
            0, "secondary", False, False,
            0, "secondary", False, False,
            html.Div(results)
        )
        
    # Start Navis classification
    navis_color = "blue"
    navis_striped = navis_animated = True
    
    navis_name = classify_navis(
        fasta=tmp_protein,
        existing_captains=fetch_captains(curated=True, with_sequence=True),
        threads=1
    )
    
    if not navis_name:
        return (
            "Classification Stopped at Navis Stage",
            100, "green", False, False,
            0, "secondary", False, False,
            0, "secondary", False, False,
            html.Div(results)
        )
        
    # Navis classification succeeded
    navis_progress = 100
    navis_striped = navis_animated = False
    results.append(dmc.Alert(
        title="Navis Classification",
        children=f"Navis: {navis_name}",
        color="blue",
        variant="light",
        className="mb-3"
    ))
    
    # Start Haplotype classification
    haplotype_color = "violet"
    haplotype_striped = haplotype_animated = True
    
    try:
        ships_df = fetch_ships(curated=True, with_sequence=True)
        haplotype_name = classify_haplotype(
            fasta=temp_files["fasta"],
            existing_ships=ships_df,
            navis=navis_name
        )
        
        if haplotype_name:
            haplotype_progress = 100
            haplotype_striped = haplotype_animated = False
            results.append(dmc.Alert(
                title="Haplotype Classification",
                children=f"Haplotype: {haplotype_name}",
                color="violet",
                variant="light",
                className="mb-3"
            ))
    except Exception as e:
        logger.error(f"Haplotype classification error: {str(e)}")
        return (
            "Classification Stopped at Haplotype Stage",
            100, "green", False, False,
            100, "blue", False, False,
            0, "secondary", False, False,
            html.Div(results)
        )
    
    # Cleanup temporary files
    try:
        os.unlink(temp_files["fasta"])
        if temp_files.get("protein"):
            os.unlink(temp_files["protein"])
    except:
        pass
    
    if not results:
        results.append(dmc.Alert(
            title="No Classifications Found",
            children="Could not classify the sequence at any level",
            color="yellow",
            variant="light"
        ))
    
    return (
        "Classification Complete",
        family_progress, family_color, family_striped, family_animated,
        navis_progress, navis_color, navis_striped, navis_animated,
        haplotype_progress, haplotype_color, haplotype_striped, haplotype_animated,
        html.Div(results)
    )
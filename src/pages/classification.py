# dash page for classification workflow

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import dash_core_components as dcc

from dash import html, callback, Input, Output, State, MATCH
from dash.dependencies import Output, Input, State
from dash.exceptions import PreventUpdate

import os
import logging
import tempfile
import pandas as pd
from typing import Callable
import time

from src.utils.seq_utils import guess_seq_type, load_fasta_to_dict
from src.config.settings import BLAST_DB_PATHS
from src.utils.classification_utils import WORKFLOW_STAGES, create_workflow_stage_callback, create_workflow_stage_callback_old
from src.tasks import check_contained_matches_task, check_similar_matches_task, check_exact_matches_task, run_family_classification_task, run_navis_classification_task, run_haplotype_classification_task
from src.components.callbacks import create_file_upload
from src.components.error_boundary import handle_callback_error
from src.utils.seq_utils import parse_fasta_from_file, write_fasta
from src.database.sql_manager import fetch_captains, fetch_ships


dash.register_page(__name__)

logger = logging.getLogger(__name__)

def make_progress_bar(message, id_prefix, color):
    """Create a progress bar with given prefix and color"""
    return dmc.Group([
        dmc.Text(message, size="sm"),
        dbc.Progress(
            id=f"{id_prefix}-progress",
            value=0,
            color=color,
            animated=False,
            striped=False,
            style={"width": "100%", "marginBottom": "5px"}
        ),
        html.Div(id=f"{id_prefix}-progress-spacer")  # Spacer div
    ])


layout = dmc.Container(
    size="md",
    children=[
        dmc.Space(h=20),
        dmc.Title("Classification Workflow"),

        # Single data store with all workflow state
        dcc.Store(id="classification-upload"),
        dcc.Store(id="classification-stage", data="Upload a sequence"),
        dcc.Store(id="classification-workflow-state", data={"current_stage": None, "complete": False}),
        
        # Store for each workflow stage
        *[
            dcc.Store(id={"type": "classification-stage-data", "index": stage['id']}) 
            for stage in WORKFLOW_STAGES
        ],

        dmc.Stack([
            # File upload section
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
            
            # File info display and submit button
            html.Div(
                id="classification-file-info",
                style={"display": "none"}
            ),
            
            # Submit button - initially hidden
            dmc.Button(
                "Start Classification",
                id="classification-submit-button",
                color="blue",
                size="md",
                style={"display": "none"}
            ),
            
            # Progress section - initially hidden
            dmc.Stack([
                dmc.Group([
                    dmc.Text("Current Stage:", size="lg", fw=500),
                    dmc.Text(id="classification-stage-display", size="lg", c="blue")
                ]),
                # Progress bar for overall progress
                dmc.Group([
                    dmc.Text("Classification Progress:", size="sm"),
                    dbc.Progress(
                        id="classification-progress",
                        value=0,
                        color="blue",
                        animated=False,
                        striped=False,
                        style={"width": "100%", "marginBottom": "5px"}
                    ),
                ]),
                # Individual stage progress indicators
                *[make_progress_bar(
                    stage["label"], 
                    f"classification-{stage['id']}", 
                    stage["color"]
                ) for stage in WORKFLOW_STAGES],
            ], gap="md", id="classification-progress-section", style={"display": "none"}),
        ], gap="md"),
        html.Div(id='classification-output', className="mt-4"),
        
        # Interval for polling workflow state
        dcc.Interval(
            id="classification-workflow-interval",
            interval=1000,  # 1 second
            disabled=True
        ),
    ]
)

# First callback: Handle file upload
@callback(
    [
        Output("classification-upload", "data"),
        Output("classification-file-info", "children"),
        Output("classification-file-info", "style"),
        Output("classification-submit-button", "style"),
    ],
    [
        Input("classification-fasta-upload", "contents"),
        Input("classification-fasta-sequence-upload", "contents"),
    ],
    prevent_initial_call=True
)
@handle_callback_error
def handle_file_upload(seq_content, upload_contents):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    logger.info(f"File upload triggered by: {trigger_id}")
    
    # Use whichever content is provided
    content_to_use = seq_content if seq_content is not None else upload_contents
    
    if content_to_use is None:
        raise PreventUpdate
        
    # Parse and validate FASTA
    header, seq, fasta_error = parse_fasta_from_file(content_to_use)
    if fasta_error or not seq:
        raise ValueError(f"FASTA parsing error: {fasta_error or 'No sequence found'}")
        
    seq_type = guess_seq_type(seq)
    if not seq_type:
        raise ValueError("Could not determine sequence type")

    # Save sequence to temporary file
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_fasta({"query_sequence": seq}, tmp_fasta)
        
    # Prepare upload data
    upload_data = {
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
        }
        
    # Create file info display
    file_info = dmc.Alert(
        title="Sequence uploaded successfully",
        children=[
            html.P(f"Sequence header: {header}"),
            html.P(f"Sequence type: {seq_type}"),
            html.P(f"Sequence length: {len(seq)} bp"),
        ],
        color="green",
        variant="light",
    )
        
    return [
        upload_data,
        file_info,
        {"display": "block"},  # Show file info
        {"display": "block"},  # Show submit button
    ]

# Define a class for tracking workflow status
class WorkflowStatus:
    def __init__(self):
        self.start_time = None
        self.current_stage = None
        self.current_stage_idx = None
        self.stage_progress = None
        self.found_match = False
        self.match_stage = None
        self.match_result = None
        self.complete = False
        self.error = None
        self.stage_values = [0] * len(WORKFLOW_STAGES)
        self.stage_animated = [False] * len(WORKFLOW_STAGES)
        self.stage_striped = [True] * len(WORKFLOW_STAGES)
        
    def to_dict(self):
        """Convert status to dictionary for storing in dcc.Store"""
        return {
            "start_time": self.start_time,
            "current_stage": self.current_stage,
            "current_stage_idx": self.current_stage_idx,
            "stage_progress": self.stage_progress,
            "found_match": self.found_match,
            "match_stage": self.match_stage,
            "match_result": self.match_result,
            "complete": self.complete,
            "error": self.error,
            "stage_values": self.stage_values,
            "stage_animated": self.stage_animated,
            "stage_striped": self.stage_striped,
        }

# Global variable to store workflow status 
# This is the key to sharing progress between the long-running process and the UI
workflow_tracker = None

# Second callback: Initialize workflow and start the background process
@callback(
    [
        Output("classification-workflow-state", "data"),
        Output("classification-progress-section", "style"),
        Output("classification-workflow-interval", "disabled"),
    ],
    Input("classification-submit-button", "n_clicks"),
    [
        State("classification-upload", "data"),
        State("classification-workflow-interval", "disabled"),
    ],
    prevent_initial_call=True
)
def run_workflow_background(n_clicks, upload_data, interval_disabled):
    if n_clicks is None or upload_data is None:
        raise PreventUpdate
    
    global workflow_tracker
    workflow_tracker = WorkflowStatus()
    workflow_tracker.start_time = time.time()
    
    logger.info("Starting classification workflow in background")
    
    try:
        # Process each stage in sequence
        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]
            
            # Update status to show progress
            workflow_tracker.current_stage = stage_id
            workflow_tracker.current_stage_idx = i
            workflow_tracker.stage_progress = 0
            workflow_tracker.stage_values[i] = 10  # Starting progress
            workflow_tracker.stage_animated[i] = True
            
            logger.info(f"Processing stage {i+1}/{len(WORKFLOW_STAGES)}: {stage_id}")
            
            # Define task functions based on stage
            result = None
            
            if stage_id == "exact":
                result = check_exact_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Exact match result: {result}")
            
            elif stage_id == "contained":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = check_contained_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Contained match result: {result}")
            
            elif stage_id == "similar":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = check_similar_matches_task.delay(
                    fasta=upload_data["fasta"],
                    ships_dict=fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Similar match result: {result}")
            
            elif stage_id == "family":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = run_family_classification_task.delay(
                    fasta=upload_data["fasta"],
                    seq_type=upload_data["seq_type"],
                    db_list=fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Family classification result: {result}")
                
                # For family stage, if result is None, stop the pipeline with appropriate message
                if result is None:
                    logger.info("Family classification did not find a match, stopping pipeline")
                    workflow_tracker.stage_values[i] = 100  # Mark this stage as completed
                    workflow_tracker.stage_animated[i] = False
                    workflow_tracker.complete = True
                    workflow_tracker.match_stage = "family"  # Not really a match but for display purposes
                    workflow_tracker.match_result = "No family match found"
                    # Not setting error as this is an expected outcome
                    break
                
                # Store family result for use in subsequent stages
                if isinstance(result, dict) and "protein" in result:
                    upload_data["protein_file"] = result["protein"]
            
            elif stage_id == "navis":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                
                # Check if we have a protein file from the family stage
                if "protein_file" not in upload_data or not upload_data["protein_file"]:
                    logger.warning("No protein file available for navis classification")
                    workflow_tracker.stage_values[i] = 100
                    workflow_tracker.stage_animated[i] = False
                    workflow_tracker.complete = True
                    workflow_tracker.match_stage = "navis"
                    workflow_tracker.match_result = "No protein data available for navis"
                    break
                
                result = run_navis_classification_task.delay(
                    fasta=upload_data["protein_file"],
                    existing_ships=fetch_captains(**upload_data["fetch_captain_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Navis classification result: {result}")
                
                # For navis stage, if result is None, stop the pipeline with appropriate message
                if result is None:
                    logger.info("Navis classification did not find a match, stopping pipeline")
                    workflow_tracker.stage_values[i] = 100
                    workflow_tracker.stage_animated[i] = False
                    workflow_tracker.complete = True
                    workflow_tracker.match_stage = "navis"
                    workflow_tracker.match_result = "No navis match found"
                    break
            
            elif stage_id == "haplotype":
                workflow_tracker.stage_progress = 30
                workflow_tracker.stage_values[i] = 30
                result = run_haplotype_classification_task.delay(
                    fasta=upload_data["fasta"],
                    existing_ships=fetch_ships(**upload_data["fetch_ship_params"]).to_dict('records'),
                    navis=fetch_captains(**upload_data["fetch_captain_params"]).to_dict('records')
                ).get(timeout=180)
                logger.debug(f"Haplotype classification result: {result}")
                
                # For haplotype stage, if result is None, stop the pipeline with appropriate message
                if result is None:
                    logger.info("Haplotype classification did not find a match, stopping pipeline")
                    workflow_tracker.stage_values[i] = 100
                    workflow_tracker.stage_animated[i] = False
                    workflow_tracker.complete = True
                    workflow_tracker.match_stage = "haplotype"
                    workflow_tracker.match_result = "No haplotype match found"
                    break
            
            else:
                logger.warning(f"No task defined for stage {stage_id}, skipping")
                workflow_tracker.stage_values[i] = 100
                workflow_tracker.stage_animated[i] = False
                continue
            
            # Update the status to show completion
            workflow_tracker.stage_progress = 100
            workflow_tracker.stage_values[i] = 100
            workflow_tracker.stage_animated[i] = False
            
            # For the exact, contained, and similar stages, continue to next stage even if result is None
            # For other stages, if we got a result, we found a match
            if result:
                if stage_id in ["family", "navis", "haplotype"]:
                    # Match found for one of the classification stages
                    logger.info(f"Match found in stage {stage_id}: {result}")
                    workflow_tracker.found_match = True
                    workflow_tracker.match_stage = stage_id
                    workflow_tracker.match_result = result
                    break
                
        # Mark as complete if we haven't already done so
        workflow_tracker.complete = True
    
    except Exception as e:
        # Handle errors - only unexpected errors
        error_message = f"Error during classification: {str(e)}"
        logger.error(error_message)
        logger.exception("Full traceback:")  # Log full traceback for debugging
        workflow_tracker.error = error_message
        workflow_tracker.complete = True
        
        # Mark current stage as failed
        if workflow_tracker.current_stage_idx is not None:
            i = workflow_tracker.current_stage_idx
            workflow_tracker.stage_values[i] = 100
            workflow_tracker.stage_animated[i] = False
            workflow_tracker.stage_striped[i] = False
    
    # Return the final state
    return workflow_tracker.to_dict(), {"display": "block"}, False

# Third callback: Update UI periodically based on current status
@callback(
    [
        Output("classification-stage", "data", allow_duplicate=True),
        Output("classification-stage-display", "children", allow_duplicate=True),
        Output("classification-progress", "value", allow_duplicate=True),
        Output("classification-progress", "animated", allow_duplicate=True),
        Output("classification-progress", "striped"),
        *[Output(f"classification-{stage['id']}-progress", "value") for stage in WORKFLOW_STAGES],
        *[Output(f"classification-{stage['id']}-progress", "animated") for stage in WORKFLOW_STAGES],
        *[Output(f"classification-{stage['id']}-progress", "striped") for stage in WORKFLOW_STAGES],
        Output("classification-output", "children"),
    ],
    Input("classification-workflow-interval", "n_intervals"),
    prevent_initial_call=True
)
def update_ui_from_status(n_intervals):
    """Update UI based on the current workflow status."""
    global workflow_tracker
    
    if workflow_tracker is None:
        raise PreventUpdate
    
    # Get current status
    status = workflow_tracker
    
    # Calculate stage display message
    stage_message = "Starting classification..."
    result_display = None
    
    if status.complete:
        if status.error:
            # True error condition
            stage_message = "Error occurred during classification"
            result_display = dmc.Alert(
                title="Classification Error",
                children=status.error,
                color="red",
                variant="light",
                className="mb-3"
            )
        elif status.found_match:
            # Successful match
            stage_info = next((s for s in WORKFLOW_STAGES if s["id"] == status.match_stage), {})
            stage_message = f"{stage_info.get('label', 'Match')} found: {status.match_result}"
            result_display = dmc.Alert(
                title=f"{stage_info.get('label', 'Match')} Found",
                children=f"Found {status.match_stage} match: {status.match_result}",
                color=stage_info.get("color", "blue"),
                variant="light",
                className="mb-3"
            )
        elif status.match_stage:
            # No match found, but we know which stage stopped the pipeline
            stage_info = next((s for s in WORKFLOW_STAGES if s["id"] == status.match_stage), {})
            no_match_message = f"No {status.match_stage} match found"
            
            if status.match_result and "No" in status.match_result:
                no_match_message = status.match_result
                
            stage_message = f"Classification stopped: {no_match_message}"
            result_display = dmc.Alert(
                title=f"Classification Result",
                children=no_match_message,
                color="yellow",
                variant="light",
                className="mb-3"
            )
        else:
            # Generic no match case (all stages completed)
            stage_message = "Classification complete - No matches found"
            result_display = dmc.Alert(
                title="Classification Complete",
                children="No matches were found for this sequence",
                color="yellow",
                variant="light",
                className="mb-3"
            )
    elif status.current_stage:
        # Still processing a stage
        stage_info = next((s for s in WORKFLOW_STAGES if s["id"] == status.current_stage), {})
        stage_message = f"Processing: {stage_info.get('label', status.current_stage)}"
    
    # Calculate overall progress
    if status.complete:
        # If complete, show 100%
        overall_progress = 100
        overall_animated = False
    elif status.current_stage_idx is not None:
        # Calculate progress based on current stage
        overall_progress = ((status.current_stage_idx) / len(WORKFLOW_STAGES)) * 100
        if status.stage_progress:
            # Add partial progress from current stage
            stage_contribution = (1 / len(WORKFLOW_STAGES)) * (status.stage_progress / 100)
            overall_progress += stage_contribution * 100
            overall_animated = True
        else:
            overall_progress = 0
            overall_animated = True
    
    return (
        stage_message,  # classification-stage data
        stage_message,  # classification-stage-display
        overall_progress,  # overall progress value
        overall_animated,  # overall animated state
        True,  # overall striped
        *status.stage_values,  # Individual stage progress values
        *status.stage_animated,  # Individual stage animated states
        *status.stage_striped,  # Individual stage striped states
        result_display,  # Final result display
    )
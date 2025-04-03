# dash page for classification workflow

import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc

from dash import html, callback
from dash.dependencies import Output, Input

import base64
import logging

from src.utils.seq_utils import guess_seq_type
from src.config.settings import BLAST_DB_PATHS
from src.utils.classification_utils import classify_sequence
from src.components.callbacks import create_file_upload
from src.components.error_boundary import handle_callback_error
from src.utils.seq_utils import parse_fasta

dash.register_page(__name__)

logger = logging.getLogger(__name__)

layout = dmc.Container(
    size="md",
    children=[
        html.H1("Classification Workflow"),
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
            dbc.Progress(id="classification-progress",
                children=[
                    dbc.Progress(id="classification-progress-family", value=33, color="success", bar=True),
                    dbc.Progress(id="classification-progress-navis", value=33, color="warning", bar=True),
                    dbc.Progress(id="classification-progress-haplotype", value=33, color="danger", bar=True),
                ]
            )
        ], gap="md"),

        html.Div(id='classification-output')
    ]
)

@callback(
    Output("classification-fasta-sequence-upload", "children"),
    [
        Input("classification-fasta-upload", "contents"),
        Input("classification-fasta-upload", "filename"),
    ],
)
@handle_callback_error
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
            logger.error(e)
            return html.Div(["There was an error processing this file."])


@callback(
    [
        Output('classification-progress-family', 'value'),
        Output('classification-progress-navis', 'value'),
        Output('classification-progress-haplotype', 'value'),
        Output('classification-output', 'children')
    ],
    Input('classification-fasta-upload', 'contents'),
    prevent_initial_call=True
)
@handle_callback_error
def update_output(seq_content):
    from src.utils.seq_utils import parse_fasta_from_file, write_fasta
    import tempfile
    
    # Reset progress bars
    family_progress = navis_progress = haplotype_progress = 0
    
    def update_progress(stage, value):
        nonlocal family_progress, navis_progress, haplotype_progress
        if stage == 'family':
            family_progress = value
        elif stage == 'navis':
            navis_progress = value
        elif stage == 'haplotype':
            haplotype_progress = value

    try:
        # Parse FASTA
        # "," is the delimeter for splitting content_type from content_string
        # content_type, content_string = seq_content.split(",")  
        _, seq, fasta_error = parse_fasta_from_file(seq_content)
        
        if fasta_error:
            raise ValueError(f"FASTA parsing error: {fasta_error}")
            
        if not seq:
            raise ValueError("No sequence found in input")
            
        seq_type = guess_seq_type(seq)
        if not seq_type:
            raise ValueError("Could not determine sequence type")

        # Save sequence to temporary FASTA file
        tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        write_fasta({"query_sequence": seq}, tmp_fasta)

        family_dict, navis_dict, haplotype_dict = classify_sequence(
            fasta=tmp_fasta,
            seq_type=seq_type,
            blast_df=None,
            hmmer_dict=None,
            db_list=BLAST_DB_PATHS, 
            threads=1,
            progress_callback=update_progress
        )
        
        # print(f"Sequence {seq}: Family {family_dict}, Navis {navis_dict}, Haplotype {haplotype_dict}")
        
        return family_progress, navis_progress, haplotype_progress, None
        
    except Exception as e:
        # logger.error(f"Classification error: {str(e)}")
        return 0, 0, 0, dmc.Alert(
            title="Classification Error",
            children=str(e),
            color="red",
            variant="light"
        )

@callback(
    [Output("classification-progress", "value"), Output("classification-progress", "label")],
    [Input("classification-progress-interval", "n_intervals")],
)
def update_progress(n):
    # check progress of some background process, in this example we'll just
    # use n_intervals constrained to be in 0-100
    progress = min(n % 110, 100)
    # only add text after 5% progress to ensure text isn't squashed too much
    return progress, f"{progress} %" if progress >= 5 else ""
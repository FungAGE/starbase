import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import logging
import json
import tempfile
from pathlib import Path
import base64

from clinker.classes import find_files, parse_files
from clinker.align import align_clusters
from clinker.plot import plot_data, save_html
from src.utils.gff2gbk import main as gff2gb

logger = logging.getLogger(__name__)

def create_layout():
    """Create the Dash app layout"""
    return html.Div([
        dcc.Location(id='url', refresh=False),
        dcc.Store(id='cluster-data'),
        dbc.Container(
            fluid=True,
            children=[
                dbc.Card(
                    [
                        dbc.CardBody([
                            html.H1("Gene Cluster Comparison", className="mb-4"),
                            html.P(
                                "Interactive visualization of gene clusters using clustermap.js",
                                className="lead text-muted",
                            ),
                        ])
                    ],
                    className="mb-4",
                ),
                dbc.Row([
                    dbc.Col([
                        dbc.Card([
                            dbc.CardBody([
                                # Settings
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Identity Threshold"),
                                        dbc.Input(
                                            id='identity-threshold',
                                            type="number",
                                            value=0.3,
                                            min=0,
                                            max=1,
                                            step=0.1,
                                        ),
                                    ], md=6),
                                    dbc.Col([
                                        dbc.Switch(
                                            id='use-file-order',
                                            label="Use File Order",
                                            value=False,
                                            className="mt-4"
                                        )
                                    ], md=6),
                                ], className="mb-3"),
                                # Visualization container
                                html.Div(
                                    id="clustermap-viz",
                                    style={
                                        "width": "100%",
                                        "height": "800px",
                                        "position": "relative",
                                        "overflow": "hidden",
                                        "border": "1px solid #ddd"
                                    }
                                )
                            ])
                        ])
                    ])
                ])
            ]
        )
    ])

def process_local_files(gff_paths, fasta_paths):
    """Process local GFF and FASTA files using clinker"""
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_paths = []
            
            # Process each GFF file
            for gff_path in gff_paths:
                logger.info(f"Processing GFF file: {gff_path}")
                
                # Find corresponding FASTA file
                gff_stem = Path(gff_path).stem
                fasta_path = next((p for p in fasta_paths if Path(p).stem == gff_stem), None)
                
                if fasta_path is None:
                    logger.error(f"No corresponding FASTA file found for {gff_path}")
                    raise ValueError(f"Missing FASTA file for GFF: {gff_path}")
                
                # Convert GFF to GenBank
                logger.info(f"Converting {gff_path} with {fasta_path} to GenBank")
                gb_path = gff2gb(str(gff_path), str(fasta_path))
                temp_paths.append(Path(gb_path))
                logger.info(f"Successfully converted to GenBank: {gb_path}")
            
            # Parse files using clinker
            logger.info(f"Parsing {len(temp_paths)} files with clinker...")
            try:
                clusters = parse_files(temp_paths)
                logger.info(f"Successfully parsed {len(clusters)} clusters")
                return clusters
            except Exception as e:
                logger.error(f"Error parsing files with clinker: {str(e)}")
                raise
                
    except Exception as e:
        logger.error(f"Error in process_local_files: {str(e)}", exc_info=True)
        return None

def process_gbk_files(gbk_dir):
    """Process GenBank files directly, with a limit of 4 files"""
    try:
        # Get all GenBank files from the directory
        gbk_dir = Path(gbk_dir)
        gbk_files = list(gbk_dir.glob("*.gbk"))
        
        if not gbk_files:
            logger.error(f"No GenBank files found in {gbk_dir}")
            return None
            
        # Check if there are too many files
        if len(gbk_files) > 4:
            logger.warning(f"Found {len(gbk_files)} GenBank files, exceeding limit of 4")
            raise ValueError("Too many GenBank files. Please limit to 4 files for visualization.")
            
        logger.info(f"Found {len(gbk_files)} GenBank files")
        
        # Parse files using clinker
        logger.info("Parsing files with clinker...")
        try:
            clusters = parse_files(gbk_files)
            logger.info(f"Successfully parsed {len(clusters)} clusters")
            return clusters
        except Exception as e:
            logger.error(f"Error parsing files with clinker: {str(e)}")
            raise
                
    except Exception as e:
        logger.error(f"Error in process_gbk_files: {str(e)}", exc_info=True)
        raise  # Let the callback handle the error

def create_clustermap_data(globaligner, use_file_order=False):
    """Convert clinker Globaligner to clustermap.js format"""
    data = globaligner.to_data(use_file_order=use_file_order)
    
    # Transform data structure if needed
    clusters = []
    links = []
    groups = []
    
    # Process clusters
    for cluster in data['clusters']:
        formatted_cluster = {
            "uid": cluster['uid'],
            "name": cluster['name'],
            "loci": []
        }
        
        # Process loci
        for locus in cluster['loci']:
            formatted_locus = {
                "uid": locus['uid'],
                "name": locus['name'],
                "start": locus['start'],
                "end": locus['end'],
                "genes": []
            }
            
            # Process genes
            for gene in locus['genes']:
                formatted_gene = {
                    "uid": gene['uid'],
                    "name": gene.get('label', gene['uid']),
                    "start": gene['start'],
                    "end": gene['end'],
                    "strand": gene['strand']
                }
                formatted_locus['genes'].append(formatted_gene)
                
            formatted_cluster['loci'].append(formatted_locus)
            
        clusters.append(formatted_cluster)
    
    # Process links
    for link in data['links']:
        formatted_link = {
            "query": {
                "uid": link['query']['uid'],
                "name": link['query'].get('label', link['query']['uid'])
            },
            "target": {
                "uid": link['target']['uid'],
                "name": link['target'].get('label', link['target']['uid'])
            },
            "identity": link['identity']
        }
        links.append(formatted_link)
    
    # Process groups
    if 'groups' in data:
        for group in data['groups']:
            formatted_group = {
                "uid": group['uid'],
                "label": group['label'],
                "genes": group['genes'],
                "colour": group.get('colour', '#808080'),
                "hidden": group.get('hidden', False)
            }
            groups.append(formatted_group)
            
    return {
        "clusters": clusters,
        "links": links,
        "groups": groups
    }

def custom_save_html(data, output):
    """Customized version of clinker's save_html function"""
    directory = Path(__file__).resolve().parent / "assets"

    with (directory / "index.html").open() as fp:
        html = fp.read()

    # Remove the div-floater with instructions
    html = html.replace(
        '<div id="div-floater">',
        '<div id="div-floater" style="display: none;">'
    )

    # Load and inject the required dependencies
    css_string = '<link rel="stylesheet" href="style.css"></link>'
    d3_string = '<script src="d3.min.js"></script>'
    cl_string = '<script src="clinker.js"></script>'
    cm_string = '<script src="clustermap.min.js"></script>'

    with (directory / "style.css").open() as fp:
        css = fp.read()
        html = html.replace(css_string, f"<style>{css}</style>")

    with (directory / "d3.min.js").open() as fp:
        d3 = fp.read()
        html = html.replace(d3_string, f"<script>{d3}</script>")

    with (directory / "clustermap.min.js").open() as fp:
        cm = fp.read()
        html = html.replace(cm_string, f"<script>{cm}</script>")

    with (directory / "clinker.js").open() as fp:
        # Read the original clinker.js content
        cl = fp.read()
        
        # Modify the plot configuration to disable editing and show tooltips
        cl = cl.replace(
            'const chart = ClusterMap.ClusterMap()',
            '''const chart = ClusterMap.ClusterMap()
                .config({
                    plot: {
                        scaleFactor: 30,
                    },
                    cluster: {
                        spacing: 50,
                        alignLabels: true,
                    },
                    gene: {
                        tooltips: true,
                        label: {
                            show: false,
                        }
                    },
                    legend: {
                        show: false,
                    },
                    link: {
                        show: true,
                    }
                })'''
        )
        
        # Add the data and modified clinker.js content
        cl = f"const data={json.dumps(data)};" + cl
        html = html.replace(cl_string, f"<script>{cl}</script>")

    with open(output, "w") as fp:
        fp.write(html)

def create_app():
    """Create and configure the Dash application"""
    app = dash.Dash(
        __name__,
        external_stylesheets=[
            dbc.themes.BOOTSTRAP,
            "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.8.0/font/bootstrap-icons.css",
        ],
        external_scripts=[
            "https://d3js.org/d3.v6.min.js",
            "/assets/clustermap.min.js"
        ],
        suppress_callback_exceptions=True,
        assets_folder="assets",
        title="Clustermap"
    )
    
    app.layout = create_layout()
    init_callbacks(app)
    
    return app

def init_callbacks(app):
    """Initialize Dash callbacks"""
    
    @app.callback(
        [Output("clustermap-viz", "children"),
         Output("cluster-data", "data")],
        [Input("identity-threshold", "value"),
         Input("use-file-order", "value")]
    )
    def update_visualization(identity_threshold, use_file_order):
        try:
            # Path to the directory containing GenBank files
            gbk_dir = "/home/adrian/Systematics/Starship_Database/starbase/src/database/db/ships/gbks"
            
            # Process the GenBank files
            try:
                clusters = process_gbk_files(gbk_dir)
            except ValueError as e:
                return html.Div([
                    html.Div(
                        className="alert alert-warning",
                        children=[
                            html.I(className="bi bi-exclamation-triangle me-2"),
                            str(e)
                        ]
                    )
                ]), None
            
            if not clusters:
                return html.Div("Error processing files"), None
                
            globaligner = align_clusters(*clusters, cutoff=identity_threshold or 0.3)
            data = globaligner.to_data(use_file_order=use_file_order)
            
            # Generate visualization
            with tempfile.NamedTemporaryFile(suffix='.html', delete=False, mode='w') as f:
                logger.info(f"Creating visualization file: {f.name}")
                custom_save_html(data, f.name)
                
                with open(f.name, 'r') as html_file:
                    viz_html = html_file.read()
                
                Path(f.name).unlink()
            
            viz_container = html.Div([
                html.Div(
                    id="status-container",
                    children=[
                        html.P(f"Processed {len(clusters)} clusters"),
                    ],
                    style={"marginBottom": "10px"}
                ),
                html.Iframe(
                    srcDoc=viz_html,
                    style={
                        "width": "100%",
                        "height": "800px",
                        "border": "none",
                        "overflow": "hidden"
                    }
                )
            ])
            
            return viz_container, None
            
        except Exception as e:
            logger.error(f"Error in visualization: {str(e)}")
            return html.Div([
                html.Div(
                    className="alert alert-danger",
                    children=[
                        html.I(className="bi bi-x-circle me-2"),
                        f"Error: {str(e)}"
                    ]
                )
            ]), None

app = create_app()
server = app.server

if __name__ == "__main__":
    # Set up logging
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    app.run_server(debug=True, host='0.0.0.0', port=8080)
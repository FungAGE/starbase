import logging
import json
import tempfile
from pathlib import Path

from clinker.classes import parse_files
from src.utils.gff2gbk import main as gff2gb

logger = logging.getLogger(__name__)
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
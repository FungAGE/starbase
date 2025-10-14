import logging
import json
import tempfile
from pathlib import Path
import pandas as pd

from clinker.classes import parse_files
from clinker.align import align_clusters
from src.utils.gff2gbk import main as gff2gb
from src.database.sql_manager import fetch_accession_ship

logger = logging.getLogger(__name__)

def create_temp_gbk_from_gff(accession_tag, temp_dir):
    """Create a temporary GenBank file from GFF data in the database with comprehensive validation"""
    try:
        # Handle version tags - extract base accession for database lookup
        base_accession = accession_tag.split('.')[0] if '.' in accession_tag else accession_tag
        logger.info(f"Looking up base accession: {base_accession} for display accession: {accession_tag}")
        
        ship_data = fetch_accession_ship(base_accession)
        
        # Validate sequence data exists
        if ship_data["sequence"] is None or ship_data["sequence"].empty:
            raise ValueError(f"No sequence data found for accession: {base_accession}")
        
        # Validate GFF data exists
        if ship_data["gff"] is None or ship_data["gff"].empty:
            raise ValueError(f"No GFF annotation data found for accession: {base_accession}")
        
        # Extract sequence string from DataFrame
        sequence = ship_data["sequence"]["sequence"].iloc[0]
        
        # Validate sequence quality
        if not sequence or len(sequence) < 100:
            raise ValueError(f"Sequence too short ({len(sequence) if sequence else 0} bp) for accession: {base_accession}. Minimum 100 bp required.")
        
        # Validate GFF has minimum features
        gff_df = ship_data["gff"].copy()
        if len(gff_df) < 1:
            raise ValueError(f"No GFF features found for accession: {base_accession}")
        
        # Log data quality metrics
        logger.info(f"Ship {base_accession}: seq_len={len(sequence)} bp, total_features={len(gff_df)}")
        
        # Save to temporary FASTA file
        fasta_file = temp_dir / f"{accession_tag}.fasta"
        with open(fasta_file, "w") as f:
            f.write(f">{base_accession}\n{sequence}\n")

        # columns in database used for GFF format
        gff_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

        # Create temporary GFF file with proper GFF format and add seqid column (use base accession as seqid for GFF format to match FASTA)
        gff_file = temp_dir / f"{accession_tag}.gff"
        gff_df['seqid'] = base_accession
        
        # HACK: Convert "gene" features to "CDS" for clinker visualization. This is just temporary until we have more complete annotation data in the gff table.
        # Ensure phase is set for CDS features (default to 0 if null) and score has a value (use '.' if null, as per GFF3 standard)
        gff_df['type'] = gff_df['type'].replace('gene', 'CDS')        
        gff_df['phase'] = gff_df['phase'].fillna(0)        
        gff_df['score'] = gff_df['score'].fillna('.')
        
        # Validate we have CDS features for visualization
        cds_count = len(gff_df[gff_df['type'] == 'CDS'])
        if cds_count < 1:
            logger.warning(f"No CDS features found for {base_accession}. Found types: {gff_df['type'].unique().tolist()}")
            raise ValueError(f"No CDS features found for {base_accession}. At least 1 CDS required for synteny visualization.")
        
        logger.info(f"Ship {base_accession}: {cds_count} CDS features available for visualization")
                
        gff_df = gff_df[gff_columns]        
        gff_df.to_csv(gff_file, sep="\t", index=False, header=False)
        
        # Convert GFF to GenBank
        gbk_file = temp_dir / f"{accession_tag}.gbk"
        gff2gb(str(gff_file), str(fasta_file), str(gbk_file))
        
        logger.info(f"Successfully created GenBank file for {accession_tag}")
        return str(gbk_file)
        
    except ValueError as e:
        # Re-raise validation errors for better error reporting
        logger.error(f"Validation error for {accession_tag}: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error creating GenBank file for {accession_tag}: {str(e)}", exc_info=True)
        raise ValueError(f"Failed to create GenBank file for {accession_tag}: {str(e)}")

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

def process_gbk_files(gbk_files, accession_tags=None):
    """Process GenBank files directly, with a limit of 4 files.
    If accession_tags are provided, generate GenBank files on-the-fly from GFF data.
    
    Args:
        gbk_files: List of GenBank file paths or directory path
        accession_tags: List of accession tags to generate GenBank files from GFF data
    
    Raises:
        ValueError: If validation fails or processing errors occur
    """
    try:
        gbk_file_paths = []
        failed_accessions = []
        
        # If accession_tags are provided, generate GenBank files on-the-fly
        if accession_tags:
            logger.info(f"Generating GenBank files for {len(accession_tags)} accession tags")
            
            # Create temporary directory for generated files
            temp_dir = Path(tempfile.mkdtemp())
            
            for accession_tag in accession_tags:
                try:
                    gbk_file = create_temp_gbk_from_gff(accession_tag, temp_dir)
                    if gbk_file:
                        gbk_file_paths.append(gbk_file)
                        logger.info(f"Successfully generated GenBank file for {accession_tag}")
                except ValueError as e:
                    # Validation error - collect for reporting
                    failed_accessions.append(f"{accession_tag}: {str(e)}")
                    logger.warning(f"Failed to generate GenBank file for {accession_tag}: {str(e)}")
                except Exception as e:
                    failed_accessions.append(f"{accession_tag}: Unexpected error - {str(e)}")
                    logger.error(f"Unexpected error generating GenBank file for {accession_tag}: {str(e)}")
        else:
            # Handle both directory path and list of files
            if isinstance(gbk_files, (str, Path)):
                gbk_dir = Path(gbk_files)
                gbk_files = list(gbk_dir.glob("*.gbk"))
            else:
                gbk_files = [Path(f) for f in gbk_files]
            
            gbk_file_paths = [str(f) for f in gbk_files]
        
        # Report any failures
        if failed_accessions:
            error_msg = f"Failed to process {len(failed_accessions)} accession(s):\n" + "\n".join(f"  - {err}" for err in failed_accessions)
            if not gbk_file_paths:
                # All failed
                raise ValueError(f"All accessions failed validation. {error_msg}")
            else:
                # Some succeeded, log warnings
                logger.warning(f"Partial success: {len(gbk_file_paths)} succeeded, {len(failed_accessions)} failed")
        
        if not gbk_file_paths:
            raise ValueError("No GenBank files found or generated")
            
        # Check if there are too many files
        if len(gbk_file_paths) > 4:
            raise ValueError(f"Too many GenBank files ({len(gbk_file_paths)}). Maximum 4 files allowed for visualization.")
            
        logger.info(f"Processing {len(gbk_file_paths)} GenBank files successfully")
        
        # Parse files using clinker
        logger.info("Parsing files with clinker...")
        try:
            # parse_files returns a list of Cluster objects
            clusters = parse_files(gbk_file_paths)
            logger.info(f"Successfully parsed {len(clusters)} clusters")
            
            # Use align_clusters to create a Globaligner object
            # Set a lower cutoff (0.1) for alignment to handle potentially poor translations
            # This is more permissive than the default to accommodate lower quality annotations
            logger.info("Aligning clusters with cutoff=0.1...")
            globaligner = align_clusters(*clusters, cutoff=0.1)
            logger.info(f"Successfully created globaligner with {len(globaligner.clusters)} clusters")
            return globaligner
        except Exception as e:
            logger.error(f"Error during alignment with cutoff=0.1: {str(e)}")
            # If alignment fails, try with an even lower cutoff or no alignment
            try:
                logger.warning("Retrying alignment with cutoff=0.0 (no filtering)...")
                globaligner = align_clusters(*clusters, cutoff=0.0)
                logger.info(f"Successfully aligned with cutoff=0.0")
                return globaligner
            except Exception as e2:
                logger.error(f"Failed even with cutoff=0.0: {str(e2)}")
                raise ValueError(f"Failed to align clusters: {str(e)}. This may indicate issues with the sequence or annotation quality.")
                
    except ValueError:
        # Re-raise validation errors as-is
        raise
    except Exception as e:
        logger.error(f"Unexpected error in process_gbk_files: {str(e)}", exc_info=True)
        raise ValueError(f"Failed to process GenBank files: {str(e)}")

def create_clustermap_data(globaligner, use_file_order=False):
    """Convert clinker Globaligner to clustermap.js format.
    
    This function processes the gene data from clinker and extracts gene labels
    from the 'Alias' attribute in GFF files (if available). The label extraction
    follows this priority:
    1. Direct 'label' field from gene object
    2. 'gene', 'label', or 'locus_tag' from gene.names dictionary
    3. Falls back to gene UID if no label found
    """
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
                # Try to extract gene label from various sources
                gene_label = gene.get('label')  # First try direct label
                if not gene_label:
                    # Try to get from names dict if available
                    names = gene.get('names', {})
                    gene_label = names.get('gene') or names.get('label') or names.get('locus_tag')
                if not gene_label:
                    # Fallback to UID
                    gene_label = gene['uid']
                
                # Log gene label extraction for debugging
                if gene.get('names'):
                    logger.debug(f"Gene {gene['uid']}: extracted label '{gene_label}' from names {gene.get('names')}")
                
                formatted_gene = {
                    "uid": gene['uid'],
                    "name": gene_label,
                    "start": gene['start'],
                    "end": gene['end'],
                    "strand": gene['strand'],
                    "names": gene.get('names', {})  # Include gene names for labeling
                }
                formatted_locus['genes'].append(formatted_gene)
                
            formatted_cluster['loci'].append(formatted_locus)
            
        clusters.append(formatted_cluster)
    
    # Process links
    for link in data['links']:
        formatted_link = {
            "uid": f"link_{link['query']['uid']}_{link['target']['uid']}",
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

def get_clustermap_json(gbk_files, identity_threshold=0.3, use_file_order=False):
    """Process GenBank files and return JSON data for ClusterMap.js
    
    Args:
        gbk_files: List of GenBank file paths
        identity_threshold: Minimum identity threshold for links
        use_file_order: Whether to use file order for clustering
        
    Returns:
        dict: JSON-serializable data for ClusterMap.js
    """
    try:
        # Process the GenBank files
        globaligner = process_gbk_files(gbk_files)
        if not globaligner:
            raise ValueError("Failed to process GenBank files")
        
        # Create the clustermap data
        clustermap_data = create_clustermap_data(globaligner, use_file_order)
        
        # Filter links by identity threshold
        if identity_threshold > 0:
            clustermap_data['links'] = [
                link for link in clustermap_data['links'] 
                if link['identity'] >= identity_threshold
            ]
        
        return clustermap_data
        
    except Exception as e:
        logger.error(f"Error creating clustermap JSON: {str(e)}")
        raise

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
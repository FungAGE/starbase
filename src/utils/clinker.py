import logging
import json
import tempfile
import gc
import psutil
from pathlib import Path
import pandas as pd

from clinker.classes import parse_files
from clinker.align import align_clusters
from src.utils.gff2gbk import main as gff2gb
from src.database.sql_manager import fetch_accession_ship

logger = logging.getLogger(__name__)

# Performance and memory limits
MAX_SEQUENCE_LENGTH = 10_000_000  # 10MB per sequence
MAX_GENES_PER_CLUSTER = 10_000    # Max genes per cluster
MAX_TOTAL_GENES = 50_000         # Max total genes across all clusters
MAX_CLUSTERS = 4                 # Max number of clusters (already enforced)
MEMORY_THRESHOLD_MB = 500        # Memory usage threshold

def create_temp_gbk_from_gff(accession_tag, temp_dir):
    """Create a temporary GenBank file from GFF data in the database"""
    try:
        if not accession_tag or not isinstance(accession_tag, str):
            raise ValueError(f"Invalid accession tag: {accession_tag}")

        if not temp_dir or not temp_dir.exists():
            raise ValueError(f"Invalid temporary directory: {temp_dir}")

        base_accession = accession_tag.split('.')[0] if '.' in accession_tag else accession_tag
        logger.info(f"Looking up base accession: {base_accession} for display accession: {accession_tag}")

        try:
            ship_data = fetch_accession_ship(base_accession)
        except Exception as e:
            logger.error(f"Database error fetching data for {base_accession}: {str(e)}")
            raise RuntimeError(f"Failed to retrieve data from database for accession {base_accession}")

        # Validate data integrity
        if ship_data is None:
            raise ValueError(f"No data returned for accession: {base_accession}")

        sequence_data = ship_data.get("sequence")
        gff_data = ship_data.get("gff")

        if sequence_data is None or (hasattr(sequence_data, 'empty') and sequence_data.empty):
            raise ValueError(f"No sequence data found for accession: {base_accession}")

        if gff_data is None or (hasattr(gff_data, 'empty') and gff_data.empty):
            raise ValueError(f"No GFF annotation data found for accession: {base_accession}")

        # Validate sequence data structure
        try:
            if not hasattr(sequence_data, 'iloc'):
                raise ValueError("Sequence data is not in expected DataFrame format")

            sequence = sequence_data["sequence"].iloc[0]
            if not sequence or len(sequence.strip()) == 0:
                raise ValueError(f"Empty or invalid sequence data for accession: {base_accession}")

            # Check sequence length for performance
            seq_length = len(sequence)
            if seq_length > MAX_SEQUENCE_LENGTH:
                raise ValueError(f"Sequence too long ({seq_length:,} bp) for accession {base_accession}. Maximum allowed: {MAX_SEQUENCE_LENGTH:,} bp")

            # Basic sequence validation (check for valid nucleotides)
            valid_chars = set('ATCGNURYSWKMBDHVatcgnuryswkmbdhv-')
            invalid_chars = set(sequence.upper()) - valid_chars
            if invalid_chars:
                logger.warning(f"Sequence contains non-standard characters for {base_accession}: {invalid_chars}")

            # Memory check
            try:
                process = psutil.Process()
                memory_mb = process.memory_info().rss / 1024 / 1024
                if memory_mb > MEMORY_THRESHOLD_MB:
                    logger.warning(f"High memory usage detected: {memory_mb:.1f} MB")
                    gc.collect()  # Force garbage collection
            except Exception as e:
                logger.debug(f"Could not check memory usage: {e}")

        except (KeyError, IndexError) as e:
            raise ValueError(f"Invalid sequence data structure for accession {base_accession}: {str(e)}")

        # Extract sequence string from DataFrame and save to temporary FASTA file
        fasta_file = temp_dir / f"{accession_tag}.fasta"
        try:
            with open(fasta_file, "w") as f:
                f.write(f">{base_accession}\n{sequence}\n")
        except Exception as e:
            raise RuntimeError(f"Failed to write FASTA file for {accession_tag}: {str(e)}")

        # columns in database used for GFF format
        gff_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

        # Validate GFF data structure
        try:
            if not hasattr(gff_data, 'copy'):
                raise ValueError("GFF data is not in expected DataFrame format")

            gff_df = gff_data.copy()

            # Check for required columns
            missing_cols = set(gff_columns) - set(gff_df.columns)
            if missing_cols:
                # Try to add missing columns with defaults
                for col in missing_cols:
                    if col == 'seqid':
                        gff_df[col] = base_accession
                    elif col == 'source':
                        gff_df[col] = '.'
                    elif col == 'score':
                        gff_df[col] = '.'
                    elif col == 'strand':
                        gff_df[col] = '+'
                    elif col == 'phase':
                        gff_df[col] = 0
                    else:
                        gff_df[col] = '.'
                logger.warning(f"Added missing GFF columns for {base_accession}: {missing_cols}")

        except Exception as e:
            raise ValueError(f"Invalid GFF data structure for accession {base_accession}: {str(e)}")

        # Create temporary GFF file with proper GFF format and add seqid column
        gff_file = temp_dir / f"{accession_tag}.gff"
        gff_df['seqid'] = base_accession

        # HACK: Convert "gene" features to "CDS" for clinker visualization
        # Ensure phase is set for CDS features (default to 0 if null) and score has a value
        gff_df['type'] = gff_df['type'].replace('gene', 'CDS')
        gff_df['phase'] = gff_df['phase'].fillna(0)
        gff_df['score'] = gff_df['score'].fillna('.')

        # Validate coordinate data
        try:
            gff_df['start'] = pd.to_numeric(gff_df['start'], errors='coerce')
            gff_df['end'] = pd.to_numeric(gff_df['end'], errors='coerce')
            gff_df = gff_df.dropna(subset=['start', 'end'])

            if len(gff_df) == 0:
                raise ValueError("No valid GFF features with coordinates found")

            # Check coordinate ranges
            invalid_coords = gff_df[(gff_df['start'] < 1) | (gff_df['end'] < gff_df['start']) | (gff_df['end'] > len(sequence))]
            if len(invalid_coords) > 0:
                logger.warning(f"Found {len(invalid_coords)} GFF features with invalid coordinates for {base_accession}")
                gff_df = gff_df[(gff_df['start'] >= 1) & (gff_df['end'] >= gff_df['start']) & (gff_df['end'] <= len(sequence))]

        except Exception as e:
            raise ValueError(f"Invalid coordinate data in GFF for accession {base_accession}: {str(e)}")

        try:
            gff_df = gff_df[gff_columns]
            gff_df.to_csv(gff_file, sep="\t", index=False, header=False)
        except Exception as e:
            raise RuntimeError(f"Failed to write GFF file for {accession_tag}: {str(e)}")

        # Convert GFF to GenBank
        gbk_file = temp_dir / f"{accession_tag}.gbk"
        try:
            gff2gb(str(gff_file), str(fasta_file), str(gbk_file))

            # Verify the output file was created and has content
            if not gbk_file.exists() or gbk_file.stat().st_size == 0:
                raise RuntimeError("GenBank conversion produced empty or missing file")

        except Exception as e:
            raise RuntimeError(f"Failed to convert GFF to GenBank for {accession_tag}: {str(e)}")

        logger.info(f"Successfully created GenBank file for {accession_tag}")
        return str(gbk_file)

    except ValueError as e:
        logger.error(f"Validation error creating GenBank file for {accession_tag}: {str(e)}")
        raise  # Re-raise validation errors as they indicate data issues
    except RuntimeError as e:
        logger.error(f"Runtime error creating GenBank file for {accession_tag}: {str(e)}")
        raise  # Re-raise runtime errors as they indicate system issues
    except Exception as e:
        logger.error(f"Unexpected error creating GenBank file for {accession_tag}: {str(e)}")
        raise RuntimeError(f"Unexpected error processing accession {accession_tag}: {str(e)}")

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
    """
    try:
        # Input validation
        if accession_tags is None and gbk_files is None:
            raise ValueError("Either gbk_files or accession_tags must be provided")

        if accession_tags is not None and gbk_files is not None:
            raise ValueError("Cannot provide both gbk_files and accession_tags")

        gbk_file_paths = []

        # If accession_tags are provided, generate GenBank files on-the-fly
        if accession_tags:
            if not isinstance(accession_tags, (list, tuple)):
                raise ValueError("accession_tags must be a list or tuple")

            if len(accession_tags) == 0:
                raise ValueError("accession_tags list cannot be empty")

            if len(accession_tags) > 4:
                raise ValueError("Too many accession tags. Please limit to 4 for visualization.")

            logger.info(f"Generating GenBank files for {len(accession_tags)} accession tags")

            # Create temporary directory for generated files
            temp_dir = Path(tempfile.mkdtemp())
            failed_accessions = []

            for accession_tag in accession_tags:
                try:
                    gbk_file = create_temp_gbk_from_gff(accession_tag, temp_dir)
                    if gbk_file:
                        gbk_file_paths.append(gbk_file)
                        logger.info(f"Successfully generated GenBank file for {accession_tag}")
                    else:
                        failed_accessions.append(accession_tag)
                        logger.warning(f"Failed to generate GenBank file for {accession_tag}")
                except Exception as e:
                    failed_accessions.append(accession_tag)
                    logger.error(f"Error generating GenBank file for {accession_tag}: {str(e)}")

            if failed_accessions:
                if len(failed_accessions) == len(accession_tags):
                    raise RuntimeError(f"Failed to generate GenBank files for all accessions: {failed_accessions}")
                else:
                    logger.warning(f"Failed to generate GenBank files for {len(failed_accessions)} accessions: {failed_accessions}")

        else:
            # Handle both directory path and list of files
            if isinstance(gbk_files, (str, Path)):
                gbk_dir = Path(gbk_files)
                if not gbk_dir.exists():
                    raise ValueError(f"GenBank directory does not exist: {gbk_dir}")
                if not gbk_dir.is_dir():
                    raise ValueError(f"Path is not a directory: {gbk_dir}")

                gbk_files = list(gbk_dir.glob("*.gbk"))
                if not gbk_files:
                    raise ValueError(f"No .gbk files found in directory: {gbk_dir}")

            else:
                if not isinstance(gbk_files, (list, tuple)):
                    raise ValueError("gbk_files must be a list, tuple, or directory path")

                gbk_files = [Path(f) for f in gbk_files]

            # Validate file paths
            valid_files = []
            for gbk_file in gbk_files:
                if not gbk_file.exists():
                    logger.warning(f"GenBank file does not exist: {gbk_file}")
                    continue
                if not gbk_file.is_file():
                    logger.warning(f"Path is not a file: {gbk_file}")
                    continue
                if gbk_file.stat().st_size == 0:
                    logger.warning(f"GenBank file is empty: {gbk_file}")
                    continue
                valid_files.append(gbk_file)

            gbk_file_paths = [str(f) for f in valid_files]

        if not gbk_file_paths:
            raise ValueError("No valid GenBank files found or generated")

        # Check if there are too many files
        if len(gbk_file_paths) > 4:
            raise ValueError(f"Too many GenBank files ({len(gbk_file_paths)}). Please limit to 4 files for visualization.")

        if len(gbk_file_paths) < 2:
            raise ValueError("At least 2 GenBank files are required for synteny comparison")

        logger.info(f"Processing {len(gbk_file_paths)} GenBank files")

        # Parse files using clinker
        logger.info("Parsing files with clinker...")
        try:
            # parse_files returns a list of Cluster objects
            clusters = parse_files(gbk_file_paths)
            if not clusters or len(clusters) == 0:
                raise ValueError("No clusters found in GenBank files")

            logger.info(f"Successfully parsed {len(clusters)} clusters")

            # Validate clusters have genes and check performance limits
            total_genes = 0
            for i, cluster in enumerate(clusters):
                num_genes = len(cluster.genes)
                total_genes += num_genes

                if num_genes > MAX_GENES_PER_CLUSTER:
                    logger.warning(f"Cluster {i} has {num_genes} genes, exceeding recommended limit of {MAX_GENES_PER_CLUSTER}")
                elif num_genes == 0:
                    logger.warning(f"Cluster {i} has no genes")

            if total_genes == 0:
                raise ValueError("No genes found in any of the GenBank files")

            if total_genes > MAX_TOTAL_GENES:
                raise ValueError(f"Too many genes ({total_genes:,}) across all clusters. Maximum allowed: {MAX_TOTAL_GENES:,}")

            logger.info(f"Total genes across all clusters: {total_genes}")

            # Memory check before alignment
            try:
                process = psutil.Process()
                memory_mb = process.memory_info().rss / 1024 / 1024
                if memory_mb > MEMORY_THRESHOLD_MB:
                    logger.warning(f"High memory usage before alignment: {memory_mb:.1f} MB")
                    gc.collect()
            except Exception as e:
                logger.debug(f"Could not check memory usage: {e}")

            # Use align_clusters to create a Globaligner object
            # Set a lower cutoff for alignment to handle potentially poor translations
            logger.info("Creating sequence alignments...")
            globaligner = align_clusters(*clusters, cutoff=0.1)

            # Final memory cleanup
            try:
                gc.collect()
            except Exception:
                pass

            logger.info(f"Successfully created globaligner with {len(globaligner.clusters)} clusters")
            return globaligner

        except Exception as e:
            logger.error(f"Error parsing/aligning files with clinker: {str(e)}")
            # If alignment fails, try with an even lower cutoff or no alignment
            try:
                logger.warning("Retrying alignment with lower cutoff...")
                if 'clusters' in locals() and clusters:
                    globaligner = align_clusters(*clusters, cutoff=0.05)
                    logger.info("Successfully created globaligner with lower cutoff")
                    return globaligner
                else:
                    raise ValueError("No clusters available for alignment")
            except Exception as e2:
                logger.error(f"Failed with lower cutoff: {str(e2)}")
                try:
                    logger.warning("Retrying with no sequence alignment...")
                    if 'clusters' in locals() and clusters:
                        globaligner = align_clusters(*clusters, cutoff=0.0)
                        logger.info("Successfully created globaligner with no alignment")
                        return globaligner
                    else:
                        raise ValueError("No clusters available for alignment")
                except Exception as e3:
                    logger.error(f"Failed even with no alignment: {str(e3)}")
                    raise RuntimeError(f"Unable to create synteny visualization. Alignment failed for all cutoff values. Original error: {str(e)}")

    except ValueError as e:
        logger.error(f"Validation error in process_gbk_files: {str(e)}")
        raise  # Re-raise validation errors
    except RuntimeError as e:
        logger.error(f"Runtime error in process_gbk_files: {str(e)}")
        raise  # Re-raise runtime errors
    except Exception as e:
        logger.error(f"Unexpected error in process_gbk_files: {str(e)}", exc_info=True)
        raise RuntimeError(f"Unexpected error during GenBank processing: {str(e)}")

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
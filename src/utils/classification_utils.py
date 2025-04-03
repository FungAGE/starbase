import tempfile
import subprocess
import logging
from src.utils.seq_utils import write_combined_fasta, write_multi_fasta
from typing import Optional, Tuple, Dict, Any, Callable
from src.database.sql_manager import fetch_meta_data, fetch_ships, fetch_captains
import pandas as pd
import hashlib

logger = logging.getLogger(__name__)

accession_workflow = """
########################################################
# assigning accessions
########################################################
accession format:
- normal ship accession: SBS123456
- updated ship accession: SBS123456.1

workflow:
- first check for exact matches
- then check for contained within matches
- then check for almost identical matches
- if no matches, assign a new accession

if a sequence is a truncated version of a longer sequence, assign the longer sequence accession, flag for review

- fast method for checking for exact matches: md5 hash of sequence
- method for checking for contained/highly similar sequences: k-mers
"""

classify_sequence_workflow = """
classifcation pipeline that should be used:
- classify query sequences for blast page, display results
- classify submitted sequences from submission page, assign an accession, and input into submissionsdatabase
"""

def assign_accession(sequence: str, 
                    existing_ships: pd.DataFrame = None,
                    threshold: float = 0.95) -> Tuple[str, bool]:
    """Assign an accession to a new sequence.
    
    Args:
        sequence: New sequence to assign accession to
        existing_ships: DataFrame of existing ships (optional, will fetch if None)
        threshold: Similarity threshold for "almost identical" matches
    
    Returns:
        Tuple[str, bool]: (accession, needs_review)
            - accession: assigned accession number
            - needs_review: True if sequence needs manual review
    """
    if existing_ships is None:
        existing_ships = fetch_ships(curated=True)
    
    logger.info(f"Starting accession assignment process")
    
    logger.info("Step 1: Checking for exact matches using MD5 hash...")
    exact_match = check_exact_match(sequence, existing_ships)
    if exact_match:
        logger.info(f"Found exact match: {exact_match}")
        return exact_match, False
        
    logger.info("Step 2: Checking for contained matches...")
    container_match = check_contained_match(sequence, existing_ships)
    if container_match:
        logger.info(f"Found containing match: {container_match}")
        return container_match, True  # Flag for review since it's truncated
        
    logger.info(f"Step 3: Checking for similar matches (threshold={threshold})...")
    similar_match = check_similar_match(sequence, existing_ships, threshold)
    if similar_match:
        logger.info(f"Found similar match: {similar_match}")
        return similar_match, True  # Flag for review due to high similarity
        
    logger.info("No matches found, generating new accession...")
    new_accession = generate_new_accession(existing_ships)
    logger.info(f"Generated new accession: {new_accession}")
    return new_accession, False

def check_exact_match(sequence: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence exactly matches an existing sequence using MD5 hash."""
    sequence_hash = hashlib.md5(sequence.encode()).hexdigest()
    logger.debug(f"Query sequence hash: {sequence_hash}")
    
    # Calculate hashes for existing sequences, skipping None values
    existing_hashes = {}
    skipped_count = 0
    for acc, seq in zip(existing_ships['accession_tag'], existing_ships['sequence']):
        if seq is None:
            skipped_count += 1
            logger.warning(f"Skipping null sequence for accession {acc}")
            continue
        existing_hashes[hashlib.md5(seq.encode()).hexdigest()] = acc
    
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} sequences due to null values")
    
    logger.info(f"Compared against {len(existing_hashes)} valid sequences")
    match = existing_hashes.get(sequence_hash)
    if match:
        logger.info(f"Found exact hash match: {match}")
    return match

def check_contained_match(sequence: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence is contained within any existing sequences.
    Returns accession of the longest containing sequence."""
    containing_matches = []
    processed_count = 0
    
    for _, row in existing_ships.iterrows():
        processed_count += 1
        if processed_count % 1 == 0:
            logger.debug(f"Processed {processed_count}/{len(existing_ships)} sequences")
            
        if row['sequence'] is None:
            continue
            
        if sequence in row['sequence']:
            logger.info(f"Found containing match: {row['accession_tag']} (length: {len(row['sequence'])})")
            containing_matches.append((
                len(row['sequence']),  # length for sorting
                row['accession_tag']
            ))
    
    # Sort by length descending and return longest match if any
    if containing_matches:
        containing_matches.sort(reverse=True)  # Sort by length descending
        logger.info(f"Found {len(containing_matches)} containing matches")
        logger.info(f"Selected longest match: {containing_matches[0][1]} (length: {containing_matches[0][0]})")
        return containing_matches[0][1]  # Return accession of longest match
        
    logger.info("No containing matches found")
    return None

def check_similar_match(sequence: str, 
                       existing_ships: pd.DataFrame,
                       threshold: float) -> Optional[str]:
    """Check for sequences with similarity above threshold using k-mer comparison."""
    logger.info(f"Starting similarity comparison (threshold={threshold})")

    tmp_fasta_all_ships = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_multi_fasta(existing_ships, tmp_fasta_all_ships, sequence_col='sequence', id_col='accession_tag')

    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name

    # Create temporary FASTA with new and existing sequences
    write_combined_fasta(
        sequence,
        existing_ships,
        fasta_path=tmp_fasta,
        sequence_col='sequence',
        id_col='accession_tag'
    )
    logger.debug(f"Created temporary FASTA file: {tmp_fasta}")
    
    # Calculate similarities
    similarities = calculate_similarities(
        tmp_fasta,
        tmp_fasta_all_ships,
        seq_type='nucl'
    )
    
    # # Find best match above threshold
    # best_match = None
    # best_sim = 0
    
    # for seq_id, sim_dict in similarities.items():
    #     if seq_id == 'query_sequence':  # Skip self comparison
    #         continue
    #     sim = sim_dict.get('query_sequence', 0)
    #     if sim > threshold and sim > best_sim:
    #         best_match = seq_id
    #         best_sim = sim
            
    # return best_match
    logger.debug(f"Calculated similarities for {len(similarities)} sequences")


    # Check if we have any similarities above threshold
    if 'query_sequence' in similarities:
        for acc_id, sim in similarities['query_sequence'].items():
            logger.debug(f"Similarity to {acc_id}: {sim}")
            if sim >= threshold:
                logger.info(f"Found similar match: {acc_id} (similarity: {sim})")
                return acc_id
                
    logger.info("No similar matches found above threshold")
    return None

def generate_new_accession(existing_ships: pd.DataFrame) -> str:
    """Generate a new unique accession number."""
    # Extract existing accession numbers
    existing_nums = [
        int(acc.replace('SBS', '').split('.')[0])
        for acc in existing_ships['accession_tag']
        if acc.startswith('SBS')
    ]
    
    # Check if we have existing accessions
    if not existing_nums:
        error_msg = "Problem with loading existing ships. No existing SBS accessions found in database."
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    # Find next available number
    next_num = max(existing_nums) + 1
    logger.info(f"Last used accession number: SBS{max(existing_nums):06d}")
    logger.info(f"Assigning new accession number: SBS{next_num:06d}")
        
    return f"SBS{next_num:06d}"

########################################################
# classification pipeline
########################################################

def classify_sequence(
        sequence: str, 
        blast_df: pd.DataFrame,
        hmmer_dict: Dict[str, Any],
        pident_thresh: float = 90,
        db_list: Dict[str, Any] = None,
        threads: int = 1,
        progress_callback: Callable[[str, int], None] = None
) -> Tuple[Dict[str, str], str, str]:
    """Classify a new sequence/ship from BLAST results based on comparison to existing classified sequences.
    
    1. Family assignment via BLAST/HMMER results or run hmmsearch/diamond on new sequence
    2. Navis assignment via captain sequence clustering
    3. Haplotype assignment via full sequence similarity

    Args:
        sequence: The new sequence to classify
        blast_df: The BLAST results dataframe to use for family assignment
        hmmer_dict: The HMMER results dict to use for family assignment
        pident_thresh: Minimum percentage identity for captain gene
        db_list: Database configuration and paths
        threads: Number of CPU threads to use
        progress_callback: Optional callback function to report progress
    
    Returns:
        Tuple[Dict[str, str], str, str]: (family_dict, navis_dict, haplotype_dict)
    """
    from src.utils.blast_utils import make_captain_alert, process_captain_results
    from src.utils.blast_utils import create_no_matches_alert
    from src.utils.seq_utils import guess_seq_type

    def update_progress(stage: str, percent: int):
        if progress_callback:
            progress_callback(stage, percent)
            logger.debug(f"Progress update - {stage}: {percent}%")

    logger.info("Starting sequence classification pipeline...")
    update_progress('family', 0)
    update_progress('navis', 0)
    update_progress('haplotype', 0)

    # Get existing classified sequences from database
    logger.info("Fetching existing sequences from database...")
    existing_meta = fetch_meta_data(curated=True)
    existing_ships = fetch_ships(curated=True)
    existing_captains = fetch_captains()
    logger.info("Successfully loaded existing sequences")
    update_progress('family', 20)

    if sequence:
        logger.info("Determining sequence type...")
        seq_type = guess_seq_type(sequence)
        logger.info(f"Sequence type determined: {seq_type}")
        update_progress('family', 40)
    else:
        logger.warning("No sequence provided")
        seq_type = None

    logger.info("Starting family classification...")
    family_dict = classify_family(
        sequence=sequence, 
        seq_type=seq_type, 
        blast_df=blast_df, 
        hmmer_dict=hmmer_dict, 
        db_list=db_list, 
        pident_thresh=pident_thresh, 
        threads=threads
    )

    if family_dict is None:
        logger.warning("No family assignment possible")
        update_progress('family', 0)
        return family_dict, None, None
    else:           
        logger.info(f"Family assigned: {family_dict}")
        update_progress('family', 100)
        update_progress('navis', 20)
        logger.info("Starting navis classification...")
        
        navis_dict = classify_navis(
            sequence=sequence, 
            existing_captains=existing_captains, 
            existing_meta=existing_meta, 
            threads=threads
        )
        if navis_dict is None:
            logger.warning("No navis assignment possible")
            update_progress('navis', 0)
            return family_dict, None, None
        else:
            logger.info(f"Navis assigned: {navis_dict}")
            update_progress('navis', 100)
            update_progress('haplotype', 20)
            logger.info("Starting haplotype classification...")
            
            haplotype_dict = classify_haplotype(
                sequence=sequence, 
                existing_ships=existing_ships, 
                existing_meta=existing_meta
            )
            if haplotype_dict is None:
                logger.warning("No haplotype assignment possible")
                update_progress('haplotype', 0)
                return family_dict, navis_dict, None
            else:
                logger.info(f"Haplotype assigned: {haplotype_dict}")
                update_progress('haplotype', 100)
                logger.info("Classification pipeline completed successfully")
                return family_dict, navis_dict, haplotype_dict

def classify_family(sequence=None, seq_type=None, blast_df=None, hmmer_dict=None, db_list=None, pident_thresh=90, input_eval=0.001, threads=1):
    """Uses blast results, hmmsearch or diamond to assign family based on captain gene similarity"""
    # Part 1: Family Assignment via Captain Gene
    # - if given blast or hmmer results, use those to assign family
    # - if given a sequence, run hmmsearch or diamond to assign family
    # - compare captain genes to existing captain sequences
    # - Assign family based on closest match

    from src.utils.blast_utils import run_hmmer

    family_dict = None

    # make sure that inputs are valid
    inputs_provided = sum(x is not None for x in [sequence, blast_df, hmmer_dict])
    
    if inputs_provided == 0:
        raise ValueError("Must provide one of: sequence, blast_df, or hmmer_dict")
    if inputs_provided > 1:
        raise ValueError("Can only provide one of: sequence, blast_df, or hmmer_dict")
        
    # Additional validation for sequence input
    if sequence is not None and seq_type is None:
        raise ValueError("Sequence provided but no sequence type")

    # Simple selection of top hit using blast or hmmer results
    if blast_df is not None and len(blast_df) > 0:
        # Sort by evalue (ascending) and pident (descending) to get best hits
        blast_df = blast_df.sort_values(['evalue', 'pident'], ascending=[True, False])
        top_hit = blast_df.iloc[0]
        logger.info(f"Top hit: {top_hit}")

        top_evalue = float(top_hit['evalue'])
        top_aln_length = int(top_hit['aln_length'])
        top_pident = float(top_hit['pident'])

        if top_pident >= pident_thresh:
            # look up family name from accession tag
            query_accession = top_hit['query_id']
            top_family = fetch_meta_data(accession_tag=query_accession)['familyName']
            family_dict = {"family": top_family, "aln_length": top_aln_length, "evalue": top_evalue}

    if sequence:        
        hmmer_dict = run_hmmer(
            db_list=db_list,
            query_type=seq_type,
            input_genes="tyr",
            input_eval=0.01,
            query_fasta=sequence,
            threads=2,
        )

    if hmmer_dict is not None and len(hmmer_dict) > 0:
        # Handle case where captain_results_dict is a list of results
        if isinstance(hmmer_dict, list):
            # Get the first result if available
            family_name = hmmer_dict[0].get('family_name')
            aln_length = hmmer_dict[0].get('aln_length')
            evalue = hmmer_dict[0].get('evalue')
        else:
            # Original behavior for dict
            family_name = hmmer_dict.get('family_name')
            aln_length = hmmer_dict.get('aln_length')
            evalue = hmmer_dict.get('evalue')
        family_dict = {"family": family_name, "aln_length": aln_length, "evalue": evalue}    

    return family_dict

def classify_navis(sequence: str,
                  existing_captains: pd.DataFrame, 
                  existing_meta: pd.DataFrame,
                  threads: int = 1) -> str:
    """Assign navis based on captain sequence clustering."""
    # Part 2: Navis Assignment
    # - Compare captain sequence to existing classified captains
    # - Use mmseqs clustering to group with existing navis

    # Create temporary FASTA with:
    # - Captain sequence from new sequence
    # - All existing classified captain sequences
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_combined_fasta(sequence, existing_captains, fasta_path=tmp_fasta, sequence_col='sequence', id_col='accession_tag')
        
    # Run mmseqs clustering
    clusters = mmseqs_easy_cluster(
        tmp_fasta,
        min_seq_id=0.5,
        coverage=0.25,
        threads=threads
    )

    # Find which cluster has the member "query_sequence"
    cluster_rep = None
    for cluster, members in clusters.items():
        if 'query_sequence' in members:
            cluster_rep = cluster
            break
    
    # id of new sequence is 'query_sequence'
    # look up the navis name from the cluster representative
    # ! how can we retain navis names if we are re-clustering?
    # HACK: use the existing meta data to get the navis names of captain sequences in the cluster
    navis_name = existing_meta[existing_meta['captainID'] == cluster_rep]['starship_navis'].values[0]
    if navis_name is None:
        logger.warning(f"No navis name, but sequence clusters with captainID: {cluster_rep}")

    return navis_name

def classify_haplotype(sequence: str,
                      existing_ships: pd.DataFrame,
                      existing_meta: pd.DataFrame) -> str:
    """Assign haplotype based on full sequence similarity."""
    # Part 3: Haplotype Assignment
    # - Compare full sequence to existing classified sequences
    # - Use sourmash for similarity calculation
    # - Use MCL clustering to group with existing haplotypes

    # Create temporary FASTA with:
    # - New sequence
    # - All existing classified sequences
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_combined_fasta(sequence, existing_ships, fasta_path=tmp_fasta, sequence_col='sequence', id_col='accession_tag')
        
        # Calculate similarities
    similarities = calculate_similarities(
        tmp_fasta,
        seq_type='nucl'
    )
        
    # Cluster sequences
    groups, _, _ = cluster_sequences(
        similarities,
        group_prefix="HAP",
        inflation=1.5,
        threshold=0.05
    )
    
    # groups looks like this:
    # 'sequence_id': {
    #     'group_id': 'HAP000001',
    #     'is_representative': True
    # }

    # Find which group has the member "query_sequence"
    rep_group_id = None
    for group, members in groups.items():
        if 'query_sequence' in members:
            rep_group_id = group['group_id']
            break
    
    # Return haplotype based on existing classifications in that group
    seq_ids_in_rep_group = []
    for group, members in groups.items():
        if group['group_id'] == rep_group_id:
            seq_ids_in_rep_group.append(group['sequence_id'])

    # just get the first sequence id in the group
    haplotype_name = existing_meta[existing_meta['captainID'] == seq_ids_in_rep_group[0]]['starship_haplotype'].values[0]

    if haplotype_name is None:
        logger.warning(f"No haplotype name, but sequence clusters with captainID: {seq_ids_in_rep_group[0]}")

    return groups

def mmseqs_easy_cluster(fasta_file, min_seq_id=0.5, coverage=0.25, threads=1):
    """Run mmseqs easy-cluster on input sequences.
    
    Args:
        fasta_file (str): Path to input FASTA file
        min_seq_id (float): Minimum sequence identity threshold (0-1)
        coverage (float): Minimum coverage threshold (0-1) 
        threads (int): Number of threads to use
        
    Returns:
        dict: Parsed clustering results mapping sequence IDs to cluster assignments
    """
    # Create temporary directory for mmseqs files
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Build command with all parameters
        base_name = os.path.join(tmp_dir, "results")
        cmd = [
            "mmseqs",
            "easy-cluster",
            fasta_file,
            base_name,
            tmp_dir,
            "--threads", str(threads),
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--alignment-mode", "3",
            "--cov-mode", "0", 
            "--cluster-reassign"
        ]
        
        try:
            # Run mmseqs
            logger.info(f"Running MMseqs2: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Parse results from the _cluster.tsv file
            cluster_file = f"{base_name}_cluster.tsv"
            return parse_mmseqs_results(cluster_file)
            
        except subprocess.CalledProcessError as e:
            logger.error(f"MMseqs2 clustering failed: {e.stderr}")
            raise

def parse_mmseqs_results(sequence_db, results_db):
    """Create a TSV file from MMseqs2 clustering results.
    
    Args:
        sequence_db (str): Path to input sequence database
        results_db (str): Path to results database

    Returns:
        dict: Mapping of sequence IDs to their cluster assignments
    """
    try:
        results_tsv = results_db + ".tsv" 
        cmd = [
            "mmseqs",
            "createtsv",
            sequence_db,
            sequence_db,
            results_db,
            results_tsv
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except Exception as e:
        logger.error(f"Error during conversion of MMseqs2 clustering results to TSV: {str(e)}")
        raise

    try:
        clusters = {}
        with open(results_tsv) as f:
            for line in f:
                # The resultsDB_clu.tsv file follows the following format:
                # #cluster-representative 	cluster-member
                # Q0KJ32	Q0KJ32
                # Q0KJ32	C0W539
                # Q0KJ32	D6KVP9
                # E3HQM9	E3HQM9
                # E3HQM9	F0YHT8
                cluster, member = line.strip().split('\t')
                clusters[member] = cluster
        return clusters
        
    except Exception as e:
        logger.error(f"Error parsing MMseqs2 results: {str(e)}")
        raise

def sourmash_sketch(sequence_file, sig_file, kmer_size, scaled, sketch_type='dna'):
    sketch_cmd = [
        "sourmash", "sketch",
        f"{sketch_type}",
        "--singleton",
        "--output", sig_file,
        "-p", f"k={kmer_size},scaled={scaled},noabund",
        sequence_file
    ]
    
    logger.info(f"Creating sourmash signatures: {' '.join(sketch_cmd)}")
    try:
        subprocess.run(sketch_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Sourmash sketch failed: {e.stderr}")
        raise

def sourmash_compare(ship_sketch, new_sketch, matrix_file, kmer_size, sketch_type='dna'):
    compare_cmd = [
        "sourmash", "compare",
        ship_sketch,
        new_sketch,
        f"--{sketch_type}",  # dna or protein
        "--csv", matrix_file,
        "-k", str(kmer_size)
    ]
    
    logger.info(f"Calculating pairwise similarities: {' '.join(compare_cmd)}")
    try:
        subprocess.run(compare_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Sourmash compare failed: {e.stderr}")
        raise


def calculate_similarities(ship_sequence_file, new_sequence_file, seq_type='nucl', threads=1):
    """Calculate k-mer similarity between sequences using sourmash.
    
    Args:
        sequence_file (str): Path to input FASTA file
        seq_type (str): Sequence type to compare - 'nucl' or 'prot'
        threads (int): Number of threads to use
        
    Returns:
        dict: Nested dictionary of pairwise similarities {seq1: {seq2: similarity}}
    """
    # Set sourmash parameters based on sequence type
    if seq_type == 'nucl':
        kmer_size = 510  # Default from starfish sig
        scaled = 100
        sketch_type = 'dna'
    elif seq_type == 'prot':
        kmer_size = 17
        scaled = 20
        sketch_type = 'protein'
    else:
        raise ValueError(f"Invalid sequence type: {seq_type}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        # 0. create sketch of ships
        ship_sketch_file = os.path.join(tmp_dir, "ships.sig")
        sourmash_sketch(ship_sequence_file, ship_sketch_file, kmer_size, scaled, sketch_type)

        # 1. Create signature file
        sig_file = os.path.join(tmp_dir, "sequences.sig")
        sourmash_sketch(new_sequence_file, sig_file, kmer_size, scaled, sketch_type)

        # 2. Calculate pairwise similarities
        matrix_file = os.path.join(tmp_dir, "similarity.csv")
        sourmash_compare(ship_sketch_file, sig_file, matrix_file, kmer_size, sketch_type)

        # 3. Parse similarity matrix
        return parse_similarity_matrix(matrix_file)

def parse_similarity_matrix(matrix_file):
    """Parse sourmash CSV output into pairwise similarities.
    
    Args:
        matrix_file (str): Path to sourmash CSV output
        
    Returns:
        dict: Nested dictionary of pairwise similarities
    """
    similarities = {}
    
    with open(matrix_file) as f:
        # First line contains sequence IDs
        headers = next(f).strip().split(',')
        seq_ids = headers[1:]  # First column is empty
        
        # Read similarity values
        for i, line in enumerate(f):
            values = line.strip().split(',')
            seq_id = values[0]
            similarities[seq_id] = {}
            
            # Convert similarity values to float and store
            for j, val in enumerate(values[1:]):
                target_id = seq_ids[j]
                if seq_id != target_id:  # Skip self-comparisons
                    sim = float(val)
                    if sim > 0:  # Only store non-zero similarities
                        similarities[seq_id][target_id] = sim

    return similarities

def write_similarity_file(similarities, output_file):
    """Write similarities to output file in starfish format.
    
    Args:
        similarities (dict): Nested dictionary of pairwise similarities
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        for seq1 in sorted(similarities):
            for seq2 in sorted(similarities[seq1]):
                # Format: seq1\tseq2\tsimilarity
                sim = similarities[seq1][seq2]
                f.write(f"{seq1}\t{seq2}\t{sim:.3f}\n")

def cluster_sequences(similarity_file, group_prefix="HAP", inflation=1.5, threshold=0.0, threads=1):
    """Group sequences into clusters using MCL algorithm.
    
    Args:
        similarity_file (str): Path to similarity file (from calculate_similarities)
        group_prefix (str): Prefix for naming groups
        inflation (float): MCL inflation parameter
        threshold (float): Minimum similarity threshold (0-1)
        threads (int): Number of threads to use
        
    Returns:
        tuple: (groups, node_data, edge_data) where:
            groups: dict mapping sequence IDs to group assignments
            node_data: dict of node metadata
            edge_data: dict of edge weights
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        # 1. Remap similarity values if threshold > 0
        if threshold > 0:
            remapped_sim = os.path.join(tmp_dir, "remapped.sim")
            remap_similarities(similarity_file, remapped_sim, threshold)
            input_file = remapped_sim
        else:
            input_file = similarity_file

        # 2. Run MCL clustering
        mcl_output = os.path.join(tmp_dir, "mcl_clusters.txt")
        run_mcl_clustering(input_file, mcl_output, inflation, threads)
        
        # 3. Parse results and assign group names
        groups = parse_and_name_groups(mcl_output, group_prefix)
        
        # 4. Generate node and edge data for visualization
        node_data = generate_node_data(groups)
        edge_data = generate_edge_data(similarity_file)
        
        return groups, node_data, edge_data

def remap_similarities(input_file, output_file, threshold):
    """Remap similarity values based on threshold.
    
    Following MCL manual recommendation to increase contrast between edge weights.
    Values are remapped by subtracting the threshold.
    """
    with open(input_file) as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            ref, query, sim = line.strip().split('\t')
            sim = float(sim)
            if sim < threshold:
                remapped = 0
            else:
                remapped = sim - threshold
            f_out.write(f"{ref}\t{query}\t{remapped:.3f}\n")

def run_mcl_clustering(input_file, output_file, inflation, threads):
    """Run MCL clustering algorithm."""
    cmd = [
        "mcl",
        input_file,
        "-I", str(inflation),
        "--abc",  # Input format is label pairs with weight
        "-te", str(threads),
        "-o", output_file
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"MCL clustering failed: {e.stderr}")
        raise

def parse_and_name_groups(mcl_output, prefix):
    """Parse MCL output and assign group names.
    
    Returns dict mapping sequence IDs to their group assignments.
    """
    groups = {}
    group_count = 0
    
    with open(mcl_output) as f:
        for line in f:
            members = line.strip().split('\t')
            if len(members) > 0:
                group_count += 1
                group_name = f"{prefix}{group_count:04d}"
                
                # First member is considered the representative
                for member in members:
                    groups[member] = {
                        'group_id': group_name,
                        'is_representative': member == members[0]
                    }
    
    logger.info(f"Grouped sequences into {group_count} clusters")
    return groups

def generate_node_data(groups):
    """Generate node metadata for visualization."""
    node_data = {}
    
    for seq_id, data in groups.items():
        node_data[seq_id] = {
            'id': seq_id,
            'group': data['group_id']
        }
    
    return node_data

def generate_edge_data(similarity_file):
    """Generate edge data with average similarities for visualization."""
    edges = {}
    
    with open(similarity_file) as f:
        for line in f:
            node1, node2, weight = line.strip().split('\t')
            weight = float(weight)
            
            # Sort nodes to ensure consistent edge keys
            key = tuple(sorted([node1, node2]))
            
            if key not in edges:
                edges[key] = {'sum': weight, 'count': 1}
            else:
                edges[key]['sum'] += weight
                edges[key]['count'] += 1
    
    # Calculate averages
    edge_data = {}
    for (node1, node2), data in edges.items():
        avg_weight = data['sum'] / data['count']
        edge_data[(node1, node2)] = avg_weight
    
    return edge_data

def write_cluster_files(groups, node_data, edge_data, output_prefix):
    """Write clustering results to files."""
    # Write main clustering results
    with open(f"{output_prefix}.mcl", 'w') as f:
        current_group = None
        for seq_id, data in sorted(groups.items(), key=lambda x: (x[1]['group_id'], not x[1]['is_representative'])):
            if data['group_id'] != current_group:
                if current_group is not None:
                    f.write('\n')
                f.write(f"{data['group_id']}\t{seq_id}")
                current_group = data['group_id']
            else:
                f.write(f"\t{seq_id}")
    
    # Write node data
    with open(f"{output_prefix}.nodes.txt", 'w') as f:
        f.write("id\tgroup\n")
        for node_id, data in sorted(node_data.items()):
            f.write(f"{node_id}\t{data['group']}\n")
    
    # Write edge data
    with open(f"{output_prefix}.edges.txt", 'w') as f:
        f.write("from\tto\tweight\n")
        for (node1, node2), weight in sorted(edge_data.items()):
            f.write(f"{node1}\t{node2}\t{weight:.3f}\n")
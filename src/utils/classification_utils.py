import tempfile
import subprocess
import logging
from src.utils.blast_utils import hmmsearch, mmseqs_easy_cluster, parse_hmmer, calculate_similarities, write_similarity_file, cluster_sequences, write_cluster_files, extract_gene_from_hmmer
from src.utils.seq_utils import write_combined_fasta, write_multi_fasta
from typing import Optional, Tuple, Dict, Any
from src.database.sql_manager import fetch_meta_data, fetch_ships, fetch_all_captains
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
        threads: int = 1) -> Tuple[Dict[str, str, str], str, str]:
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
    
    Returns:
        Tuple[Dict[str, str, str], str, str]: (family, navis, haplotype) assignments
    """
    from src.utils.blast_utils import make_captain_alert, process_captain_results
    from src.utils.blast_utils import create_no_matches_alert
    from src.utils.seq_utils import guess_seq_type

    logger.info("Starting sequence classification pipeline...")

    # Get existing classified sequences from database
    logger.info("Fetching existing sequences from database...")
    existing_meta = fetch_meta_data(curated=True)
    existing_ships = fetch_ships(curated=True)
    existing_captains = fetch_all_captains()
    logger.info("Successfully loaded existing sequences")

    family_dict = None
    navis_dict = None
    haplotype_dict = None

    if sequence:
        logger.info("Determining sequence type...")
        seq_type = guess_seq_type(sequence)
        logger.info(f"Sequence type determined: {seq_type}")
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
        return family_dict, None, None
    else:           
        logger.info(f"Family assigned: {family_dict}")
        logger.info("Starting navis classification...")
        
        navis_dict = classify_navis(
            sequence=sequence, 
            existing_captains=existing_captains, 
            existing_meta=existing_meta, 
            threads=threads
        )
        if navis_dict is None:
            logger.warning("No navis assignment possible")
            return family_dict, None, None
        else:
            logger.info(f"Navis assigned: {navis_dict}")
            logger.info("Starting haplotype classification...")
            
            haplotype_dict = classify_haplotype(
                sequence=sequence, 
                existing_ships=existing_ships, 
                existing_meta=existing_meta
            )
            if haplotype_dict is None:
                logger.warning("No haplotype assignment possible")
                return family_dict, navis_dict, None
            else:
                logger.info(f"Haplotype assigned: {haplotype_dict}")
                logger.info("Classification pipeline completed successfully")
                return family_dict, navis_dict, haplotype_dict

def classify_family(sequence, seq_type, blast_df, hmmer_dict, db_list, pident_thresh=90, input_eval=0.001, threads=1):
    """Uses blast results, hmmsearch or diamond to assign family based on captain gene similarity"""
    # Part 1: Family Assignment via Captain Gene
    # - if given blast or hmmer results, use those to assign family
    # - if given a sequence, run hmmsearch or diamond to assign family
    # - compare captain genes to existing captain sequences
    # - Assign family based on closest match

    from src.utils.blast_utils import run_hmmer

    family_dict = None

    # make sure that inputs are valid
    if not sequence and not blast_df and not hmmer_dict:
        raise ValueError("No sequence, blast_df, or hmmer_dict provided")
    if sequence and seq_type is None:
        raise ValueError("Sequence provided but no sequence type")
    if blast_df and hmmer_dict:
        raise ValueError("Both blast_df and hmmer_dict provided")
    if (blast_df or hmmer_dict) and sequence:
            raise ValueError("blast_df or hmmer_dict and sequence provided")

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
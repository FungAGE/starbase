import tempfile
import subprocess
import logging
from src.utils.blast_utils import hmmsearch, mmseqs_easy_cluster, parse_hmmer, calculate_similarities, write_similarity_file, cluster_sequences, write_cluster_files, extract_gene_from_hmmer
from src.utils.seq_utils import write_combined_fasta
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
    
    # Step 1: Check for exact matches using MD5 hash
    exact_match = check_exact_match(sequence, existing_ships)
    if exact_match:
        return exact_match, False
        
    # Step 2: Check if sequence is contained within existing sequence
    container_match = check_contained_match(sequence, existing_ships)
    if container_match:
        return container_match, True  # Flag for review since it's truncated
        
    # Step 3: Check for highly similar sequences
    similar_match = check_similar_match(sequence, existing_ships, threshold)
    if similar_match:
        return similar_match, True  # Flag for review due to high similarity
        
    # Step 4: If no matches, assign new accession
    new_accession = generate_new_accession(existing_ships)
    return new_accession, False


def check_exact_match(sequence: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence exactly matches an existing sequence using MD5 hash."""
    sequence_hash = hashlib.md5(sequence.encode()).hexdigest()
    
    # Calculate hashes for existing sequences
    existing_hashes = {
        hashlib.md5(seq.encode()).hexdigest(): acc 
        for acc, seq in zip(existing_ships['accession_tag'], existing_ships['sequence'])
    }
    
    return existing_hashes.get(sequence_hash)

def check_contained_match(sequence: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence is contained within any existing sequences."""
    for _, row in existing_ships.iterrows():
        if sequence in row['sequence']:
            return row['accession_tag']
    return None

def check_similar_match(sequence: str, 
                       existing_ships: pd.DataFrame,
                       threshold: float) -> Optional[str]:
    """Check for sequences with similarity above threshold using k-mer comparison."""
    # Create temporary FASTA with new and existing sequences
    tmp_fasta = write_combined_fasta(
        sequence,
        existing_ships,
        sequence_col='sequence',
        id_col='accession_tag'
    )
    
    # Calculate similarities
    similarities = calculate_similarities(
        tmp_fasta,
        mode='element',
        seq_type='nucl'
    )
    
    # Find best match above threshold
    best_match = None
    best_sim = 0
    
    for seq_id, sim_dict in similarities.items():
        if seq_id == 'query_sequence':  # Skip self comparison
            continue
        sim = sim_dict.get('query_sequence', 0)
        if sim > threshold and sim > best_sim:
            best_match = seq_id
            best_sim = sim
            
    return best_match

def generate_new_accession(existing_ships: pd.DataFrame) -> str:
    """Generate a new unique accession number."""
    # Extract existing accession numbers
    existing_nums = [
        int(acc.replace('SBS', '').split('.')[0])
        for acc in existing_ships['accession_tag']
        if acc.startswith('SBS')
    ]
    
    # Find next available number
    if existing_nums:
        next_num = max(existing_nums) + 1
    else:
        next_num = 1
        
    return f"SBS{next_num:06d}"

########################################################
# classification pipeline
########################################################

def classify_sequence(sequence: str, db_list: Dict[str, Any], threads: int = 1) -> Tuple[str, str, str]:
    """Classify a new sequence based on comparison to existing classified sequences.
    
    Args:
        sequence: The new sequence to classify
        db_list: Database configuration and paths
        threads: Number of CPU threads to use
    
    Returns:
        Tuple[str, str, str]: (family, navis, haplotype) assignments
    """
    # Get existing classified sequences from database
    existing_meta = fetch_meta_data(curated=True)
    existing_ships = fetch_ships(curated=True)
    existing_captains = fetch_all_captains()
    
    # Part 1: Family Assignment via Captain Gene
    # - First identify captain gene in new sequence via hmmsearch
    # - Compare to existing captain sequences
    # - Assign family based on closest match
    family = classify_family(sequence, db_list, existing_captains, threads)
    
    # Part 2: Navis Assignment
    # - Compare captain sequence to existing classified captains
    # - Use mmseqs clustering to group with existing navis
    navis = classify_navis(sequence, existing_captains, existing_meta, threads)
    
    # Part 3: Haplotype Assignment
    # - Compare full sequence to existing classified sequences
    # - Use sourmash for similarity calculation
    # - Use MCL clustering to group with existing haplotypes
    haplotype = classify_haplotype(sequence, existing_ships, existing_meta)
    
    return family, navis, haplotype

def classify_family(sequence, db_list, input_eval=0.001, threads=1):
    """Uses hmmsearch to assign family based on captain gene similarity"""
    # first, assign all Starships to a family by searching the candidate captain sequences against the reference captain HMM profile database and identifying the best profile hit to each sequence:
    tmp_hmmer='elementFinder/macpha6_tyr_vs_YRsuperfams.out'
    hmmer_db = db_list["gene"]["captain"]["hmm"]["tyr"]
    # TODO: add this to hmmsearch method?
    hmmer_cmd = f"hmmsearch --noali --notextw --max -o {tmp_hmmer} --cpu {threads} -E {input_eval} {hmmer_db} {sequence}"
    logger.info(f"Running HMMER search: {hmmer_cmd}")
    subprocess.run(hmmer_cmd, shell=True)
    hmmer_results = hmmsearch(db_list=db_list, 
                             query_type="captain",
                             threads=threads,
                             query_fasta=sequence)
    # Parse results and return family assignment
    # perl -p -e 's/ +/\t/g' elementFinder/macpha6_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v '#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/macpha6_tyr_vs_YRsuperfams_besthits.txt
    # TODO: make sure this parses correctly, or just use a simpler method for grabbing the family assignment from the hmmsearch output
    family_assignment = parse_hmmer(hmmer_results)

    result = extract_gene_from_hmmer(hmmer_results)
    if not result:
        return None
        
    header, captain_seq = result

    # TODO: determine what the output of this function should be

    return family_assignment

def classify_navis(sequence: str,
                  existing_captains: pd.DataFrame, 
                  existing_meta: pd.DataFrame,
                  threads: int = 1) -> str:
    """Assign navis based on captain sequence clustering."""
    # Create temporary FASTA with:
    # - Captain sequence from new sequence
    # - All existing classified captain sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta') as tmp_fasta:
        write_combined_fasta(tmp_fasta.name, sequence, existing_captains)
        
        # Run mmseqs clustering
        clusters = mmseqs_easy_cluster(
            tmp_fasta.name,
            min_seq_id=0.5,
            coverage=0.25,
            threads=threads
        )
        
        # clusters will be a dict like:
        # {
        #    'seq1': 'cluster1_rep',
        #    'seq2': 'cluster1_rep',
        #    'seq3': 'cluster2_rep',
        #    ...
        # }

        # Find which cluster contains our sequence
        # Return navis assignment based on existing classifications in that cluster

        # TODO: determine what the output of this function should be

        return navis_assignment

def classify_haplotype(sequence: str,
                      existing_ships: pd.DataFrame,
                      existing_meta: pd.DataFrame) -> str:
    """Assign haplotype based on full sequence similarity."""
    # Create temporary FASTA with:
    # - New sequence
    # - All existing classified sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta') as tmp_fasta:
        write_combined_fasta(tmp_fasta.name, sequence, existing_ships)
        
        # Calculate similarities
        similarities = calculate_similarities(
            tmp_fasta.name,
            mode='element',
            seq_type='nucl'
        )
        
        # Cluster sequences
        groups, _, _ = cluster_sequences(
            similarities,
            group_prefix="HAP",
            inflation=1.5,
            threshold=0.05
        )
        
        # Find which group contains our sequence
        # Return haplotype based on existing classifications in that group

        # TODO: determine what the output of this function should be
        return haplotype_assignment
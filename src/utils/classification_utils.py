import tempfile
import subprocess
import pandas as pd
import hashlib
import os
import glob
import signal
import shutil
import screed
from sourmash import (
    SourmashSignature,
    MinHash,
    save_signatures,
    load_file_as_signatures,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

from src.utils.seq_utils import (
    write_combined_fasta,
    write_multi_fasta,
    create_tmp_fasta_dir,
    load_fasta_to_dict,
    clean_sequence,
    revcomp,
)
from src.database.sql_manager import fetch_ships
from Bio import SeqIO
from Bio.Seq import Seq
import networkx as nx

from typing import Optional, Tuple, Dict, Any
from src.database.sql_manager import fetch_ships, fetch_captains
from src.utils.blast_data import WorkflowState, ClassificationData

accession_workflow = """
########################################################
# assigning accessions
########################################################
accession format:
- normal ship accession: SSA123456
- updated ship accession: SSA123456.1

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

# Define our workflow stages with their colors
WORKFLOW_STAGES = [
    {"id": "exact", "label": "Checking for exact matches", "color": "red"},
    {"id": "contained", "label": "Checking for contained matches", "color": "orange"},
    {"id": "similar", "label": "Checking for similar matches", "color": "yellow"},
    # {"id": "denovo", "label": "Running denovo annotation", "color": "pink"},
    {"id": "family", "label": "Running family classification", "color": "green"},
    {"id": "navis", "label": "Running navis classification", "color": "blue"},
    {"id": "haplotype", "label": "Running haplotype classification", "color": "violet"},
]

# Define badge colors based on source
source_colors = {
    "exact": "green",
    "contained": "teal",
    "similar": "cyan",
    "blast_hit": "blue",
    "hmmsearch": "purple",
    "captain clustering": "orange",
    "k-mer similarity": "red",
    "unknown": "gray",
}

# Define icon and color based on confidence level
confidence_colors = {
    "High": "green",
    "Medium": "yellow",
    "Low": "red",
    "Unknown": "gray",
}

confidence_icons = {
    "High": "mdi:shield-check",
    "Medium": "mdi:shield-half-full",
    "Low": "mdi:shield-outline",
}


def assign_accession(
    sequence: str, existing_ships: pd.DataFrame = None, threshold: float = 0.95,
    precomputed_sig_path: str = None
) -> Tuple[str, bool]:
    """Assign an accession to a new sequence.

    Args:
        sequence: New sequence to assign accession to
        existing_ships: DataFrame of existing ships (optional, will fetch if None)
        threshold: Similarity threshold for "almost identical" matches
        precomputed_sig_path: Path to pre-computed sourmash signature file (optional, for efficiency)

    Returns:
        Tuple[str, bool]: (accession, needs_review)
            - accession: assigned accession number
            - needs_review: True if sequence needs manual review
    """

    logger.debug("Starting accession assignment process")

    logger.debug("Step 1: Checking for exact matches using MD5 hash...")
    exact_match = check_exact_match(sequence, existing_ships)
    if exact_match:
        logger.debug(f"Found exact match: {exact_match}")
        return exact_match, False

    logger.debug("Step 2: Checking for contained matches...")
    container_result = check_contained_match(
        fasta=sequence,
        existing_ships=existing_ships,
        min_coverage=0.95,
        min_identity=0.95,
    )
    if container_result:
        container_match, is_perfect_match = container_result
        logger.debug(f"Found containing match: {container_match} (perfect: {is_perfect_match})")

        # If it's a perfect match (coverage=1.0, identity=1.0), treat as exact match
        if is_perfect_match:
            logger.debug("Perfect match found - treating as exact match")
            return container_match, False  # No review needed for perfect matches
        else:
            logger.debug("Imperfect match found - requires review")
            return container_match, True  # Flag for review since it's truncated or imperfect

    logger.debug(f"Step 3: Checking for similar matches (threshold={threshold})...")
    similar_match, similarities = check_similar_match(
        sequence, existing_ships, threshold, precomputed_sig_path=precomputed_sig_path
    )
    if similar_match:
        logger.debug(f"Found similar match: {similar_match}")
        return similar_match, True  # Flag for review due to high similarity

    logger.debug("No matches found, generating new accession...")
    new_accession = generate_new_accession(existing_ships)
    logger.debug(f"Generated new accession: {new_accession}")
    return new_accession, False


def generate_new_accession(existing_ships: pd.DataFrame) -> str:
    """Generate a new unique accession number."""
    # Extract existing accession numbers
    existing_nums = [
        int(acc.replace("SSA", "").split(".")[0])
        for acc in existing_ships["accession_tag"]
        if acc is not None and acc.startswith("SSA")
    ]

    # Check if we have existing accessions
    if not existing_nums:
        logger.info("No existing SSA accessions found - starting fresh accession numbering from 1")
        next_num = 1
        logger.debug(f"Assigning new accession number: SSA{next_num:06d}")
    else:
        # Find next available number
        next_num = max(existing_nums) + 1
        logger.debug(f"Last used accession number: SSA{max(existing_nums):06d}")
        logger.debug(f"Assigning new accession number: SSA{next_num:06d}")

    return f"SSA{next_num:06d}"


def get_version_sort_key(version_tag):
    """Convert version_tag to sortable format. Handle None/empty as 0."""
    if not version_tag:
        return 0
    try:
        return int(version_tag)
    except (ValueError, TypeError):
        return 0


########################################################
# sequence matching functions
########################################################


# sequence can be a string, Seq object, or SeqRecord object
def generate_md5_hash(sequence):
    if sequence is None:
        raise ValueError("Sequence is None, can't generate MD5 hash")
    if type(sequence) == SeqIO.SeqRecord:
        sequence = str(sequence.seq)
    if type(sequence) == Seq:  # Handle Seq objects
        sequence = str(sequence)
    return hashlib.md5(sequence.encode()).hexdigest()


def check_exact_match(fasta: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence exactly matches an existing sequence using MD5 hash."""
    sequence = None

    # Handle both file paths and direct sequence strings
    if os.path.exists(fasta):
        logger.debug(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None
        sequence_list = list(sequences.values())
    else:
        logger.debug(f"Treating input as sequence string, length: {len(fasta)}")
        sequence_list = [fasta]

    if not sequence_list:
        logger.error("No sequence provided")
        return None

    # generate md5 for query sequence(s)
    # create nested dicts for query and query_revcomp
    query_md5 = {}
    for seq in sequence_list:
        if seq is None:
            logger.warning("Skipping None sequence")
            continue
        clean_seq = clean_sequence(seq)
        if clean_seq is None:
            logger.warning(f"clean_sequence returned None for sequence: {seq[:50]}...")
            continue
        md5_hash = generate_md5_hash(clean_seq)
        md5_hash_revcomp = generate_md5_hash(revcomp(clean_seq))
        if md5_hash is None or md5_hash_revcomp is None:
            logger.warning(
                f"generate_md5_hash returned None for sequence: {seq[:50]}..."
            )
            continue
        query_md5[seq] = md5_hash
        # Also store reverse complement hash for checking
        query_md5[f"{seq}_revcomp"] = md5_hash_revcomp

    # collect existing md5 in dict
    # Create reverse mapping: md5 -> accession_display (includes version)
    existing_hashes = {}
    sequences_without_md5 = []

    stored_md5_count = 0
    stored_rev_comp_count = 0
    calculated_md5_count = 0
    
    for _, row in existing_ships.iterrows():
        # Use accession_display if available, otherwise fall back to accession_tag
        accession_display = row.get("accession_display", row.get("accession_tag"))

        if row.get("md5") is not None and accession_display is not None:
            existing_hashes[row["md5"]] = accession_display
            stored_md5_count += 1
        # Also check reverse complement MD5 if available in database
        if row.get("rev_comp_md5") is not None and accession_display is not None:
            existing_hashes[row["rev_comp_md5"]] = accession_display
            stored_rev_comp_count += 1
        # Calculate MD5 on the fly for sequences without stored MD5 (both md5 and rev_comp_md5 are None)
        if (row.get("md5") is None and row.get("rev_comp_md5") is None and 
            row.get("sequence") is not None and accession_display is not None):
            # Calculate MD5 on the fly for sequences without stored MD5
            clean_seq = clean_sequence(row["sequence"])
            if clean_seq:
                calculated_md5 = generate_md5_hash(clean_seq)
                calculated_md5_revcomp = generate_md5_hash(revcomp(clean_seq))
                if calculated_md5:
                    existing_hashes[calculated_md5] = accession_display
                    sequences_without_md5.append(accession_display)
                    calculated_md5_count += 1
                if calculated_md5_revcomp:
                    existing_hashes[calculated_md5_revcomp] = accession_display
                    calculated_md5_count += 1
    
    logger.debug(f"MD5 hash sources: stored={stored_md5_count}, stored_rev_comp={stored_rev_comp_count}, calculated={calculated_md5_count}")

    logger.debug(f"Found {len(existing_hashes)} existing MD5 hashes in database")
    if sequences_without_md5:
        logger.debug(
            f"Calculated MD5 for {len(sequences_without_md5)} sequences on-the-fly"
        )
    logger.debug(f"Query MD5 hashes: {list(query_md5.values())}")

    # Check if query hash exists in database
    for seq, md5 in query_md5.items():
        match = existing_hashes.get(md5)
        if match:
            logger.info(f"Found exact hash match: {match}")
            return match

    logger.debug("No exact match found")
    return None


def check_contained_match(
    fasta: str,
    existing_ships: pd.DataFrame,
    min_coverage: float = 0.95,
    min_identity: float = 0.95,
) -> Optional[Tuple[str, bool]]:
    """Check if sequence is contained within any existing sequences.

    Args:
        sequence: Query sequence to check
        existing_ships: DataFrame containing existing sequences
        min_coverage: Minimum coverage of query sequence required (default: 0.95)
        min_identity: Minimum sequence identity required (default: 0.95)

    Returns:
        Tuple of (accession_tag, is_perfect_match) of the best containing match,
        or None if no match found. is_perfect_match is True if coverage=1.0 and identity=1.0.
    """
    containing_matches = []

    if os.path.exists(fasta):
        logger.debug(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None
        sequence = list(sequences.values())[0]
    else:
        # If fasta is not a file path, treat it as a sequence string
        logger.debug(f"Treating input as sequence string, length: {len(fasta)}")
        sequence = fasta

    if not sequence:
        logger.error("No sequence provided")
        return None

    query_len = len(sequence)

    logger.debug(f"Checking for contained matches (query length: {query_len})")

    # ruff: noqa
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as query_file:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as ref_file:
            logger.debug("Writing query sequence to temporary file")
            query_file.write(f">query\n{sequence}\n")
            query_file.flush()

            ref_count = 0
            logger.debug("Writing reference sequences to temporary file")
            for _, row in existing_ships.iterrows():
                if row["sequence"] is not None and len(row["sequence"]) >= query_len:
                    # Use accession_display if available, otherwise fall back to accession_tag
                    accession_display = row.get(
                        "accession_display", row.get("accession_tag")
                    )
                    ref_file.write(f">{accession_display}\n{row['sequence']}\n")
                    ref_count += 1
            ref_file.flush()
            logger.debug(f"Written {ref_count} reference sequences for comparison")

            logger.debug("Running minimap2 alignment")
            cmd = f"minimap2 -c --cs -t 1 {ref_file.name} {query_file.name}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.returncode != 0:
                logger.error(f"minimap2 failed with error: {result.stderr}")
                return None

            alignment_count = 0
            for line in result.stdout.splitlines():
                alignment_count += 1
                fields = line.split("\t")
                if len(fields) < 10:
                    continue

                ref_name = fields[5]
                matches = int(fields[9])
                align_length = int(fields[10])

                coverage = align_length / query_len
                identity = matches / align_length

                if coverage >= min_coverage and identity >= min_identity:
                    logger.debug(
                        f"Found containing match: {ref_name} "
                        f"(coverage: {coverage:.2f}, identity: {identity:.2f})"
                    )

                    # Find the sequence length using either accession_display or accession_tag
                    matching_rows = existing_ships[
                        (
                            existing_ships.get(
                                "accession_display", existing_ships.get("accession_tag")
                            )
                            == ref_name
                        )
                    ]
                    if matching_rows.empty:
                        # Fallback to accession_tag only
                        matching_rows = existing_ships[
                            existing_ships["accession_tag"] == ref_name
                        ]

                    if not matching_rows.empty:
                        sequence_length = len(matching_rows["sequence"].iloc[0])
                    # Check if this is a perfect match (coverage=1.0 and identity=1.0)
                    is_perfect_match = (coverage >= 1.0 and identity >= 1.0)

                    containing_matches.append(
                        (
                            identity * coverage,  # score for sorting
                            sequence_length,  # length for tiebreaking
                            ref_name,
                            is_perfect_match,  # whether this is a perfect match
                        )
                    )

            logger.debug(f"Processed {alignment_count} alignments from minimap2")

    # Sort by score descending, then by length descending
    if containing_matches:
        containing_matches.sort(reverse=True)
        logger.debug(f"Found {len(containing_matches)} containing matches")
        logger.debug(
            f"Selected best match: {containing_matches[0][2]} "
            f"(score: {containing_matches[0][0]:.2f}, length: {containing_matches[0][1]}, "
            f"perfect_match: {containing_matches[0][3]})"
        )
        return containing_matches[0][2], containing_matches[0][3]  # Return (accession, is_perfect_match)

    logger.debug("No containing matches found")
    return None


def check_similar_match(
    fasta: str, existing_ships: pd.DataFrame, threshold: float, 
    precomputed_sig_path: str = None
) -> Tuple[Optional[str], Any]:
    """
    Check for sequences with similarity above threshold using k-mer comparison.
    
    Args:
        fasta: Query sequence or path to FASTA file
        existing_ships: DataFrame of existing ships
        threshold: Similarity threshold
        precomputed_sig_path: Path to pre-computed signature file for existing ships
        
    Returns:
        Tuple of (match_accession, similarities_dict) or (None, None)
    """
    logger.debug(f"Starting similarity comparison (threshold={threshold})")

    # Handle sequence input
    if os.path.exists(fasta):
        logger.debug(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None, None
        sequence = list(sequences.values())[0]
    else:
        # If fasta is not a file path, treat it as a sequence string
        logger.debug(f"Treating input as sequence string, length: {len(fasta)}")
        sequence = fasta

    if not sequence:
        logger.error("No sequence provided")
        return None, None

    # Create temp FASTA with ONLY the query sequence
    tmp_query_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    with open(tmp_query_fasta, 'w') as f:
        f.write(f">query_sequence\n{sequence}\n")
    
    # Try to load pre-computed signatures for existing ships
    existing_signatures = None
    if precomputed_sig_path:
        logger.debug(f"Attempting to load pre-computed signatures from {precomputed_sig_path}")
        try:
            from src.database.blastdb import load_sourmash_signatures
            existing_signatures = load_sourmash_signatures(precomputed_sig_path.replace('.sig', ''))
            if existing_signatures:
                logger.info(f"Loaded {len(existing_signatures)} pre-computed signatures")
        except Exception as e:
            logger.warning(f"Failed to load pre-computed signatures: {e}, will compute on-the-fly")
    
    # If no pre-computed signatures, create them from existing ships
    if existing_signatures is None:
        logger.debug("No pre-computed signatures available, creating from existing ships DataFrame")
        tmp_existing_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
        write_multi_fasta(
            existing_ships,
            tmp_existing_fasta,
            sequence_col="sequence",
            id_col="accession_display" if "accession_display" in existing_ships.columns else "accession_tag",
        )
        existing_signatures = sourmash_sketch(tmp_existing_fasta, "nucl")
        os.unlink(tmp_existing_fasta)
    
    # Calculate similarities using pre-computed or newly created signatures
    similarities = calculate_similarities_with_precomputed(
        query_fasta=tmp_query_fasta,
        existing_signatures=existing_signatures,
        seq_type="nucl"
    )
    
    # Clean up temp file
    os.unlink(tmp_query_fasta)
    
    # Find matches above threshold
    query_id = "query_sequence"
    
    if isinstance(similarities, dict) and query_id in similarities:
        for match_id, sim in similarities[query_id].items():
            logger.debug(f"Similarity to {match_id}: {sim}")
            if sim >= threshold:
                logger.debug(f"Found similar match: {match_id} (similarity: {sim})")
                return match_id, similarities
    
    logger.debug("No similar matches found above threshold")
    return None, None


########################################################
# classification pipeline
########################################################


def classify_family(
    fasta=None,
    seq_type=None,
    meta_dict=None,
    pident_thresh=90,
    input_eval=0.001,
    threads=1,
):
    """Uses blast results, hmmsearch or diamond to assign family based on captain gene similarity
    Part 1: Family Assignment via Captain Gene
    - if given blast or hmmer results, use those to assign family
    - if given a sequence, run hmmsearch or diamond to assign family
    - compare captain genes to existing captain sequences
    - Assign family based on closest match
    Returns:
        - For nucleotide input (seq_type=="nucl"): (family_dict, protein_file)
        - For protein input (seq_type=="prot"): (family_dict, None)
    """
    from src.utils.blast_utils import run_hmmer, select_ship_family

    family_dict = None
    tmp_protein_filename = None
    hmmer_dict = None

    if meta_dict is not None and isinstance(meta_dict, list):
        meta_df = pd.DataFrame(meta_dict)

    if os.path.exists(fasta):
        logger.debug(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None, None
    else:
        # If fasta is not a file path, treat it as a sequence string
        logger.debug(f"Treating input as sequence string, length: {len(fasta)}")
        # For classify_family, we need to create a temporary file if we have a sequence string
        # but since run_hmmer expects a file path, this will be handled by run_hmmer itself
        pass

    hmmer_dict, tmp_protein_filename = run_hmmer(
        query_type=seq_type,
        input_gene="tyr",
        input_eval=0.01,
        query_fasta=fasta,
        threads=2,
    )

    if hmmer_dict is not None:
        hmmer_df = pd.DataFrame(hmmer_dict)
        if len(hmmer_df) > 0:
            family_name, family_aln_length, family_evalue = select_ship_family(hmmer_df)
            if family_name:
                family_dict = {
                    "n_hits": len(hmmer_df),
                    "family": family_name,
                    "aln_length": family_aln_length,
                    "evalue": family_evalue,
                }

    # Return based on sequence type
    if seq_type == "nucl":
        return family_dict, tmp_protein_filename
    else:
        return family_dict, None


def classify_navis(
    fasta: str, existing_captains: pd.DataFrame, threads: int = 1
) -> str:
    """Assign navis based on captain sequence clustering.
    - Compare captain sequence to existing classified captains
    - Use mmseqs clustering to group with existing navis

    Create temporary dir with FASTAs from:
    - Captain sequence from new sequence
    - All existing classified captain sequences
    """

    logger.debug("Starting navis classification")
    tmp_fasta_dir = create_tmp_fasta_dir(fasta, existing_captains)
    logger.debug(f"Created temporary FASTA directory: {tmp_fasta_dir}")

    # Run mmseqs clustering
    clusters = mmseqs_easy_cluster(
        tmp_fasta_dir,
        output_dir=os.path.join(tmp_fasta_dir, "clusters"),
        min_seq_id=0.5,
        coverage=0.25,
        threads=threads,
    )
    logger.debug(f"Clustering complete, got {len(clusters)} results")

    # The clusters dict maps member sequences to their representatives
    # So we can directly look up the query sequence
    if "query_sequence" not in clusters:
        logger.warning("Query sequence not found in clustering results")
        return None

    # Get the representative sequence for our query's cluster
    cluster_rep = clusters["query_sequence"]
    logger.debug(
        f"Query sequence belongs to cluster with representative: {cluster_rep}"
    )

    if cluster_rep == "query_sequence":
        logger.warning("Query sequence is the cluster representative")
        return None
    else:
        # Look up the navis name from the cluster representative
        # Try captainID first, then accession_display, then accession_tag
        matching_captains = existing_captains[
            existing_captains["captainID"] == cluster_rep
        ]

        # If not found and accession_display column exists, try that too
        if matching_captains.empty and "accession_display" in existing_captains.columns:
            matching_captains = existing_captains[
                existing_captains["accession_display"] == cluster_rep
            ]

        # If still not found and accession_tag column exists, try that too
        if matching_captains.empty and "accession_tag" in existing_captains.columns:
            matching_captains = existing_captains[
                existing_captains["accession_tag"] == cluster_rep
            ]

        if matching_captains.empty:
            logger.warning(
                f"No matching captain found for cluster representative: {cluster_rep}"
            )
            return None

        navis_name = matching_captains["navis_name"].iloc[0]
        logger.debug(f"Found navis name: {navis_name}")

        return navis_name


def classify_haplotype(fasta, existing_ships, navis=None, similarities=None):
    """
    Classify a sequence by haplotype based on sequence similarity.

    Args:
        fasta (str): Path to the FASTA file containing the query sequence.
        existing_ships (DataFrame or list): DataFrame or list of existing ships data.
        navis (str, list, or dict): The navis value to filter ships by. If None, all ships are used.

    Returns:
        dict: Classification result with haplotype information.
    """
    try:
        # Convert to DataFrame if needed
        if isinstance(existing_ships, list):
            logger.debug("Converting existing_ships from list to DataFrame")
            existing_ships = pd.DataFrame(existing_ships)

        # Print the type and value of navis for debugging
        logger.debug(f"navis type: {type(navis)}, value: {navis}")

        # Handle different types of navis parameter
        if navis is not None:
            if isinstance(navis, list):
                if len(navis) > 0:
                    navis = str(navis[0])
                    logger.debug(f"Using first value from navis list: {navis}")
                else:
                    navis = None
                    logger.warning("Empty navis list provided, using all ships")
            elif isinstance(navis, dict):
                if len(navis) > 0:
                    first_key = next(iter(navis))
                    navis = str(navis[first_key])
                    logger.debug(f"Using first value from navis dict: {navis}")
                else:
                    navis = None
                    logger.warning("Empty navis dict provided, using all ships")
            else:
                navis = str(navis)
                logger.debug(f"Using navis value: {navis}")

        if existing_ships.empty:
            logger.warning("No existing ships data provided")
            return {
                "haplotype_name": "Unknown",
                "confidence": 0,
                "note": "No ships data available",
            }

        # Check if navis_name column exists
        if navis is not None and "navis_name" in existing_ships.columns:
            filtered_ships = existing_ships[existing_ships["navis_name"] == navis]
            logger.debug(f"Filtered to {len(filtered_ships)} ships with navis={navis}")
            if filtered_ships.empty:
                logger.warning(
                    f"No ships found with navis={navis}, using all available ships"
                )
                filtered_ships = existing_ships
        else:
            if navis is not None:
                logger.warning(
                    "No navis_name column in ships DataFrame, using all ships"
                )
            filtered_ships = existing_ships

        # Check if necessary columns exist
        if "sequence" not in filtered_ships.columns:
            logger.error("No 'sequence' column in ships DataFrame")
            return {
                "haplotype_name": "Unknown",
                "confidence": 0,
                "note": "Missing sequence data",
            }

        required_columns = ["sequence", "captainID", "haplotype_name"]
        missing_columns = [
            col for col in required_columns if col not in filtered_ships.columns
        ]
        if missing_columns:
            logger.warning(f"Missing columns in ships DataFrame: {missing_columns}")
            return {
                "haplotype_name": "Unknown",
                "confidence": 0,
                "note": f"Missing required columns: {missing_columns}",
            }

        # Cluster sequences
        if similarities is None:
            logger.warning("No similarities provided, cannot perform clustering")
            return {
                "haplotype_name": "Unknown",
                "confidence": 0,
                "note": "No similarity data available for clustering",
            }
        
        try:
            groups = cluster_sequences(similarities, threshold=0.95)
            logger.debug(f"Clustered sequences into {len(groups)} groups")
        except Exception as e:
            logger.error(f"Error clustering sequences: {e}")
            return {
                "haplotype_name": "Error",
                "confidence": 0,
                "note": f"Clustering error: {str(e)}",
            }

        # Find the group containing the query sequence
        query_group = None
        query_seq = SeqIO.read(fasta, "fasta")
        for group in groups:
            if query_seq.id in group:
                query_group = group
                break

        if query_group is None:
            logger.warning("Query sequence not found in any cluster")
            return {
                "haplotype_name": "Novel",
                "confidence": 0,
                "note": "Query did not cluster with any existing sequence",
            }

        # Get ship IDs in the same group
        ship_ids = [seq_id for seq_id in query_group if seq_id != query_seq.id]
        logger.debug(f"Found {len(ship_ids)} ships in the same cluster as query")

        if not ship_ids:
            logger.warning("No ships in the same cluster as query")
            return {
                "haplotype_name": "Novel",
                "confidence": 0,
                "note": "Query formed its own cluster",
            }

        # Get haplotypes for these ships
        # Match by accession_display since ship_ids are sequence IDs
        if "accession_display" in filtered_ships.columns:
            matching_ships = filtered_ships[
                filtered_ships["accession_display"].isin(ship_ids)
            ]
            ship_haplotypes = matching_ships["haplotype_name"]
        else:
            # Fallback to accession_tag if accession_display not available
            matching_ships = filtered_ships[
                filtered_ships["accession_tag"].isin(ship_ids)
            ]
            ship_haplotypes = matching_ships["haplotype_name"]
        

        if ship_haplotypes.empty:
            logger.warning("No haplotype information for clustered ships")
            return {
                "haplotype_name": "Novel",
                "confidence": 0,
                "note": "No haplotype information available",
            }

        # Count haplotype occurrences
        haplotype_counts = ship_haplotypes.value_counts()
        logger.debug(f"Haplotype counts: {haplotype_counts.to_dict()}")

        if haplotype_counts.empty:
            logger.warning("No haplotype counts found")
            return {
                "haplotype_name": "Novel",
                "confidence": 0,
                "note": "No haplotype information available",
            }

        # Get the most common haplotype
        most_common_haplotype = haplotype_counts.index[0]
        confidence = haplotype_counts.iloc[0] / haplotype_counts.sum()

        logger.debug(
            f"Most common haplotype: {most_common_haplotype} with confidence {confidence:.2f}"
        )

        return {
            "haplotype_name": most_common_haplotype,
            "confidence": float(confidence),
            "counts": haplotype_counts.to_dict(),
            "cluster_size": len(query_group) - 1,  # Excluding the query itself
        }
    except Exception as e:
        logger.error(f"Unexpected error in classify_haplotype: {e}")
        logger.exception("Full traceback:")
        return {
            "haplotype_name": "Error",
            "confidence": 0,
            "note": f"Classification error: {str(e)}",
        }


def mmseqs_easy_cluster(
    fasta_dir: str, output_dir: str, min_seq_id=0.5, coverage=0.25, threads=1
) -> Dict[str, str]:
    """Run mmseqs easy-cluster on input sequences.

    Args:
        fasta_dir: Directory containing FASTA files
        output_dir: Directory for output files
        min_seq_id: Minimum sequence identity threshold (0-1)
        coverage: Minimum coverage threshold (0-1)
        threads: Number of threads to use

    Returns:
        dict: Parsed clustering results mapping sequence IDs to cluster assignments
    """
    logger.debug("Starting sequence clustering with MMseqs2")
    logger.debug(f"Input directory: {fasta_dir}")
    logger.debug(f"Output directory: {output_dir}")

    # Get all FASTA files in directory
    fasta_files = glob.glob(os.path.join(fasta_dir, "*.fa"))
    if not fasta_files:
        raise ValueError(f"No .fa files found in {fasta_dir}")
    logger.debug(f"Found {len(fasta_files)} FASTA files")

    # Create directories for MMseqs2 output
    clusters_dir = os.path.join(output_dir, "clusters")
    tmp_dir = os.path.join(clusters_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # Run easy-cluster directly on FASTA files
    results_prefix = os.path.join(clusters_dir, "results")
    cluster_cmd = [
        "mmseqs",
        "easy-cluster",
        *fasta_files,  # Pass all FASTA files directly
        results_prefix,  # This is just the prefix, MMseqs2 will add _cluster.tsv
        tmp_dir,  # Temporary directory for MMseqs2
        "--threads",
        str(threads),
        "--min-seq-id",
        str(min_seq_id),
        "-c",
        str(coverage),
        "--alignment-mode",
        "3",
        "--cov-mode",
        "0",
        "--cluster-reassign",
        "--createdb-mode",
        "0",  # Add this to ensure proper database creation
    ]

    process = None
    try:
        # Start process with a new process group
        process = subprocess.Popen(
            cluster_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,  # This ensures the process has its own session
        )

        # Wait for process with timeout
        stdout, stderr = process.communicate(timeout=300)  # 5 minute timeout

        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                process.returncode, cluster_cmd, stdout, stderr
            )

        if stderr:
            logger.debug(f"STDERR: {stderr}")

        # Parse clustering results from the TSV file
        cluster_file = f"{results_prefix}_cluster.tsv"
        logger.debug(f"Looking for cluster results in: {cluster_file}")

        if not os.path.exists(cluster_file):
            raise FileNotFoundError(f"Cluster results file not found: {cluster_file}")

        clusters = {}
        with open(cluster_file) as f:
            for line in f:
                rep_seq, member_seq = line.strip().split("\t")
                clusters[member_seq] = rep_seq

        logger.debug(f"Found {len(set(clusters.values()))} clusters")
        return clusters

    except subprocess.TimeoutExpired:
        logger.error("MMseqs2 process timed out")
        if process:
            # Kill the entire process group
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                # Give it a moment to terminate gracefully
                process.wait(timeout=5)
            except Exception as e:
                logger.warning(f"Error while terminating process gracefully: {str(e)}")
                # If it doesn't terminate gracefully, force kill
                try:
                    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                except Exception as e:
                    logger.warning(f"Error while force killing process: {str(e)}")
        raise

    except subprocess.CalledProcessError as e:
        logger.error(f"MMseqs2 command failed: {e.stderr}")
        raise

    except Exception as e:
        logger.error(f"Error during sequence clustering: {str(e)}")
        logger.exception("Full traceback:")
        raise

    finally:
        # Cleanup process if it's still running
        if process and process.poll() is None:
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGTERM)
                process.wait(timeout=5)
            except Exception as e:
                logger.warning(f"Error while terminating process gracefully: {str(e)}")
                try:
                    os.killpg(os.getpgid(process.pid), signal.SIGKILL)
                except Exception as e:
                    logger.warning(f"Error while force killing process: {str(e)}")
                    pass


def sourmash_sketch(fasta_file, seq_type="nucl"):
    """
    Create sourmash signatures directly from a FASTA file without intermediate files.

    Args:
        fasta_file (str): Path to the FASTA file containing sequences
        seq_type (str): Type of sequences, either 'nucl' or 'prot'

    Returns:
        list: List of (sequence_id, signature) tuples
    """
    try:
        # Set parameters based on sequence type
        if seq_type == "nucl":
            k = 21
            scaled = 1000
            is_protein = False
        else:
            k = 7
            scaled = 100
            is_protein = True

        # Create signatures for each sequence
        signatures = []

        for record in screed.open(fasta_file):
            # Create MinHash object
            mh = MinHash(n=0, ksize=k, scaled=scaled, is_protein=is_protein)

            # Add sequence data
            if is_protein:
                mh.add_protein(record.sequence)
            else:
                mh.add_sequence(record.sequence, force=True)

            # Create signature
            sig = SourmashSignature(mh, name=record.name)
            signatures.append((record.name, sig))

        return signatures

    except Exception as e:
        logger.error(f"Error creating sourmash signatures: {e}")
        return []


def calculate_similarities(fasta_file, seq_type="nucl", restricted_comparisons=None):
    """
    Calculate pairwise similarities between sequences in a FASTA file directly
    using the sourmash API without creating intermediate files.

    Args:
        fasta_file (str): Path to the FASTA file containing sequences
        seq_type (str): Type of sequences, either 'nucl' or 'prot'
        restricted_comparisons (dict): Dictionary of form {seq_id1: {seq_id2: True}}
                                       for comparisons to skip

    Returns:
        dict: Nested dictionary of similarities {seq_id1: {seq_id2: similarity}}
    """
    logger.debug(f"Directly calculating similarities for {fasta_file}, type={seq_type}")

    # Initialize restricted comparisons if not provided
    if restricted_comparisons is None:
        restricted_comparisons = {}

    try:
        # Get signatures directly
        signatures = sourmash_sketch(fasta_file, seq_type)

        if not signatures:
            logger.error("Failed to create signatures")
            return {}

        # Extract all sequence IDs
        all_seq_ids = [seq_id for seq_id, _ in signatures]

        # Initialize similarity dictionary with zeros
        similarities = {}
        for seq_id1 in all_seq_ids:
            similarities[seq_id1] = {}
            for seq_id2 in all_seq_ids:
                if seq_id1 != seq_id2:
                    similarities[seq_id1][seq_id2] = 0.0

        # Calculate pairwise similarities
        observed_comparisons = set()

        for i, (seq_id1, sig1) in enumerate(signatures):
            for j, (seq_id2, sig2) in enumerate(signatures[i + 1 :], i + 1):
                # Skip self-comparisons and restricted comparisons
                if seq_id1 == seq_id2:
                    continue

                # Check if this comparison is restricted
                if (
                    seq_id1 in restricted_comparisons
                    and seq_id2 in restricted_comparisons[seq_id1]
                ) or (
                    seq_id2 in restricted_comparisons
                    and seq_id1 in restricted_comparisons[seq_id2]
                ):
                    logger.debug(
                        f"Skipping restricted comparison between {seq_id1} and {seq_id2}"
                    )
                    continue

                # Calculate Jaccard similarity directly from signature objects
                similarity = sig1.jaccard(sig2)

                # Store sorted to ensure consistent keys (like in the Perl version)
                seq1, seq2 = sorted([seq_id1, seq_id2])
                similarities[seq1][seq2] = similarity
                similarities[seq2][seq1] = similarity  # Store symmetrically

                # Mark as observed
                observed_comparisons.add((seq1, seq2))
                observed_comparisons.add((seq2, seq1))

        # Ensure all valid comparisons have an entry (already initialized to 0.0)
        logger.debug(
            f"Calculated {len(observed_comparisons) / 2} pairwise similarities"
        )
        return similarities

    except Exception as e:
        logger.error(f"Error in direct similarity calculation: {e}")
        return {}


def calculate_similarities_with_precomputed(
    query_fasta, existing_signatures, seq_type="nucl", restricted_comparisons=None
):
    """
    Calculate similarities between query sequences and pre-computed signatures.
    
    This is much more efficient than calculate_similarities() when you have
    pre-computed signatures for existing sequences.
    
    Args:
        query_fasta: Path to FASTA file with query sequence(s)
        existing_signatures: List of (seq_id, signature) tuples for existing sequences
        seq_type: 'nucl' or 'prot'
        restricted_comparisons: Dict of comparisons to skip
        
    Returns:
        Dict of similarities {seq_id1: {seq_id2: similarity}}
    """
    logger.debug(f"Calculating similarities with pre-computed signatures")
    
    if restricted_comparisons is None:
        restricted_comparisons = {}
    
    try:
        # Create signatures for query sequences only
        query_signatures = sourmash_sketch(query_fasta, seq_type)
        
        if not query_signatures:
            logger.error("Failed to create signatures for query")
            return {}
        
        # Combine signatures
        all_signatures = query_signatures + existing_signatures
        all_seq_ids = [seq_id for seq_id, _ in all_signatures]
        
        logger.debug(f"Comparing {len(query_signatures)} query vs {len(existing_signatures)} existing signatures")
        
        # Initialize similarity dictionary
        similarities = {}
        for seq_id1 in all_seq_ids:
            similarities[seq_id1] = {}
            for seq_id2 in all_seq_ids:
                if seq_id1 != seq_id2:
                    similarities[seq_id1][seq_id2] = 0.0
        
        # Calculate pairwise similarities
        observed_comparisons = set()
        
        for i, (seq_id1, sig1) in enumerate(all_signatures):
            for j, (seq_id2, sig2) in enumerate(all_signatures[i + 1:], i + 1):
                if seq_id1 == seq_id2:
                    continue
                
                # Check if comparison is restricted
                if ((seq_id1 in restricted_comparisons and 
                     seq_id2 in restricted_comparisons[seq_id1]) or
                    (seq_id2 in restricted_comparisons and 
                     seq_id1 in restricted_comparisons[seq_id2])):
                    logger.debug(f"Skipping restricted comparison: {seq_id1} vs {seq_id2}")
                    continue
                
                # Calculate Jaccard similarity
                similarity = sig1.jaccard(sig2)
                
                # Store symmetrically
                seq1, seq2 = sorted([seq_id1, seq_id2])
                similarities[seq1][seq2] = similarity
                similarities[seq2][seq1] = similarity
                
                observed_comparisons.add((seq1, seq2))
                observed_comparisons.add((seq2, seq1))
        
        logger.debug(f"Calculated {len(observed_comparisons) / 2} pairwise similarities")
        return similarities
        
    except Exception as e:
        logger.error(f"Error in similarity calculation with pre-computed signatures: {e}")
        return {}


def cluster_sequences(similarities, threshold=0.95):
    """
    Cluster sequences based on pairwise similarities.

    Args:
        similarities (list): List of tuples (seq_id1, seq_id2, similarity)
        threshold (float): Similarity threshold for clustering

    Returns:
        list: List of sets, where each set contains sequence IDs in the same cluster
    """
    logger.debug(f"Clustering sequences with threshold {threshold}")

    try:
        # Filter similarities by threshold
        filtered_similarities = [
            (id1, id2) for id1, id2, sim in similarities if sim >= threshold
        ]
        logger.debug(
            f"Using {len(filtered_similarities)} out of {len(similarities)} similarities above threshold"
        )

        if not filtered_similarities:
            logger.warning("No similarities above threshold")
            return []

        # Create a graph of connected sequences
        G = nx.Graph()

        # Add all sequence IDs to the graph
        all_seq_ids = set()
        for id1, id2, _ in similarities:
            all_seq_ids.add(id1)
            all_seq_ids.add(id2)

        G.add_nodes_from(all_seq_ids)

        # Add edges for pairs with similarity >= threshold
        G.add_edges_from(filtered_similarities)

        # Find connected components (clusters)
        clusters = list(nx.connected_components(G))
        logger.debug(f"Found {len(clusters)} clusters")

        return clusters

    except Exception as e:
        logger.error(f"Error in cluster_sequences: {e}")
        return []


def remap_similarities(input_file, output_file, threshold):
    """Remap similarity values based on threshold.

    Following MCL manual recommendation to increase contrast between edge weights.
    Values are remapped by subtracting the threshold.
    """
    with open(input_file) as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            ref, query, sim = line.strip().split("\t")
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
        "-I",
        str(inflation),
        "--abc",  # Input format is label pairs with weight
        "-te",
        str(threads),
        "-o",
        output_file,
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
            members = line.strip().split("\t")
            if len(members) > 0:
                group_count += 1
                group_name = f"{prefix}{group_count:04d}"

                # First member is considered the representative
                for member in members:
                    groups[member] = {
                        "group_id": group_name,
                        "is_representative": member == members[0],
                    }

    logger.debug(f"Grouped sequences into {group_count} clusters")
    return groups


def generate_node_data(groups):
    """Generate node metadata for visualization."""
    node_data = {}

    for seq_id, data in groups.items():
        node_data[seq_id] = {"id": seq_id, "group": data["group_id"]}

    return node_data


def generate_edge_data(similarity_file):
    """Generate edge data with average similarities for visualization."""
    edges = {}

    with open(similarity_file) as f:
        for line in f:
            node1, node2, weight = line.strip().split("\t")
            weight = float(weight)

            # Sort nodes to ensure consistent edge keys
            key = tuple(sorted([node1, node2]))

            if key not in edges:
                edges[key] = {"sum": weight, "count": 1}
            else:
                edges[key]["sum"] += weight
                edges[key]["count"] += 1

    # Calculate averages
    edge_data = {}
    for (node1, node2), data in edges.items():
        avg_weight = data["sum"] / data["count"]
        edge_data[(node1, node2)] = avg_weight

    return edge_data


def write_cluster_files(groups, node_data, edge_data, output_prefix):
    """Write clustering results to files."""
    # Write main clustering results
    with open(f"{output_prefix}.mcl", "w") as f:
        current_group = None
        for seq_id, data in sorted(
            groups.items(),
            key=lambda x: (x[1]["group_id"], not x[1]["is_representative"]),
        ):
            if data["group_id"] != current_group:
                if current_group is not None:
                    f.write("\n")
                f.write(f"{data['group_id']}\t{seq_id}")
                current_group = data["group_id"]
            else:
                f.write(f"\t{seq_id}")

    # Write node data
    with open(f"{output_prefix}.nodes.txt", "w") as f:
        f.write("id\tgroup\n")
        for node_id, data in sorted(node_data.items()):
            f.write(f"{node_id}\t{data['group']}\n")

    # Write edge data
    with open(f"{output_prefix}.edges.txt", "w") as f:
        f.write("from\tto\tweight\n")
        for (node1, node2), weight in sorted(edge_data.items()):
            f.write(f"{node1}\t{node2}\t{weight:.3f}\n")

def metaeuk_createdb(query_fasta):
    """Run MetaEuk createdb for reference database."""
    
    tmp_dir = tempfile.mkdtemp()
    try:
        cmd = [
            "metaeuk",
            "createdb",
            os.path.abspath(query_fasta),
            os.path.join(tmp_dir, "db"),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return os.path.join(tmp_dir, "db"), tmp_dir
    except Exception:
        # Clean up on failure
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise

def metaeuk_easy_predict(query_fasta, threads=20):
    """Run MetaEuk easy-predict for de novo annotation.

    Args:
        query_fasta: Path to input FASTA file
        threads: Number of threads to use
    
    Returns:
        Tuple of (codon_fasta_path, fasta_path, gff_path) or (None, None, None) on error
    """

    if not os.path.exists(query_fasta):
        logger.warning("Query FASTA file does not exist")
        return None, None, None

    if not isinstance(query_fasta, str) or query_fasta == "" or query_fasta is None:
        logger.warning("query_fasta should be a path to a fasta file")
        return None, None, None  

    ref_db = None
    ref_db_tmp_dir = None
    
    try:
        ref_db, ref_db_tmp_dir = metaeuk_createdb(query_fasta)
    except Exception as e:
        logger.error(f"Error in metaeuk_createdb: {e}")
        return None, None, None

    try:
        with tempfile.TemporaryDirectory() as output_dir:
            # Output prefix for MetaEuk files
            output_prefix = os.path.join(output_dir, "prediction")
            
            # Run MetaEuk with explicit paths
            cmd = [
                "metaeuk",
                "easy-predict",
                os.path.abspath(query_fasta),
                os.path.abspath(ref_db),
                output_prefix,
                "--metaeuk-eval",
                "0.0001",
                "-e",
                "100",
                "--max-seqs",
                "1",
                "--min-length",
                "40",
                "--search-type",
                "3",
                "--threads",
                str(threads),
            ]

            # Run command and capture output
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.debug(f"MetaEuk command completed successfully")

            # Check if output files were created with correct paths
            codon_fasta = f"{output_prefix}.codon.fas"
            fasta = f"{output_prefix}.fas"
            gff = f"{output_prefix}.gff"

            missing_files = []
            for file_path in [codon_fasta, fasta, gff]:
                if not os.path.exists(file_path):
                    missing_files.append(file_path)

            if missing_files:
                logger.error(f"Expected output files not created: {missing_files}")
                logger.error(f"MetaEuk stdout: {result.stdout}")
                logger.error(f"MetaEuk stderr: {result.stderr}")
                raise FileNotFoundError(f"Expected output files not created: {missing_files}")

            # Copy files to a permanent location outside the temp directory
            permanent_dir = tempfile.mkdtemp(prefix="metaeuk_output_")
            
            permanent_codon = os.path.join(permanent_dir, "prediction.codon.fas")
            permanent_fasta = os.path.join(permanent_dir, "prediction.fas") 
            permanent_gff = os.path.join(permanent_dir, "prediction.gff")
            
            shutil.copy2(codon_fasta, permanent_codon)
            shutil.copy2(fasta, permanent_fasta)
            shutil.copy2(gff, permanent_gff)

            return permanent_codon, permanent_fasta, permanent_gff

    except subprocess.CalledProcessError as e:
        logger.error(f"MetaEuk easy-predict failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error during MetaEuk easy-predict: {str(e)}")
        logger.exception("Full traceback:")
        raise
    finally:
        # Clean up reference database temporary directory
        if ref_db_tmp_dir and os.path.exists(ref_db_tmp_dir):
            shutil.rmtree(ref_db_tmp_dir, ignore_errors=True)


def run_classification_workflow(workflow_state, blast_data, classification_data, meta_dict=None):
    """Run the classification workflow and return results."""
    import pandas as pd
    
    # Initialize similarities to None to avoid NameError in haplotype stage
    similarities = None

    try:
        # Parse parameters for database fetches
        ships_df = fetch_ships(
            curated=workflow_state.fetch_ship_params.curated,
            with_sequence=workflow_state.fetch_ship_params.with_sequence,
            dereplicate=workflow_state.fetch_ship_params.dereplicate,
        )
        captains_df = fetch_captains(
            curated=workflow_state.fetch_captain_params.curated,
            with_sequence=workflow_state.fetch_captain_params.with_sequence,
        )

        fasta_path = blast_data.fasta_file
        if isinstance(fasta_path, dict) and "content" in fasta_path:
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta_path["content"])
            blast_data.fasta_file = tmp_file
            fasta_path = tmp_file

        # Make sure blast_df is serializable
        if blast_data.blast_df and isinstance(blast_data.blast_df, pd.DataFrame):
            blast_data.blast_df = blast_data.blast_df.to_dict("records")

        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]

            # Update state to show progress
            workflow_state.current_stage = stage_id
            workflow_state.current_stage_idx = i

            workflow_state.stages[stage_id]["progress"] = 10
            workflow_state.stages[stage_id]["status"] = "running"

            logger.debug(f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}")

            if stage_id == "exact":
                logger.debug("Running exact match check")

                result = check_exact_match(
                    fasta=blast_data.fasta_file, existing_ships=ships_df
                )

                if result:
                    logger.debug(f"Found exact match: {result}")
                    workflow_state.stages[stage_id]["progress"] = 100
                    workflow_state.stages[stage_id]["status"] = "complete"

                    classification_data.source = "exact"
                    classification_data.closest_match = result
                    classification_data.confidence = "High"
                    classification_data.match_details = f"Exact sequence match to {result}"
                    
                    workflow_state.set_classification(classification_data)
                    workflow_state.complete = True
                    logger.debug(
                        f"Exact match found - stage: {workflow_state.match_stage}, result: {classification_data}"
                    )
                    return workflow_state

            if stage_id == "contained":
                logger.debug("Running contained match check")
                result = check_contained_match(
                    fasta=blast_data.fasta_file,
                    existing_ships=ships_df,
                    min_coverage=0.95,
                    min_identity=0.95,
                )

                if result:
                    match_accession, is_perfect_match = result
                    logger.debug(f"Found contained match: {match_accession} (perfect: {is_perfect_match})")
                    workflow_state.stages[stage_id]["progress"] = 30
                    workflow_state.stages[stage_id]["status"] = "complete"

                    # Set source to "exact" for perfect matches, "contained" for imperfect matches
                    classification_data.source = "exact" if is_perfect_match else "contained"
                    classification_data.closest_match = match_accession
                    classification_data.confidence = "High" if is_perfect_match else "Medium"

                    # Update workflow state to reflect the correct stage
                    workflow_state.match_stage = classification_data.source
                    if is_perfect_match:
                        classification_data.match_details = f"Exact sequence match to {match_accession}"
                    else:
                        classification_data.match_details = f"Query sequence contained within {match_accession}"
                    
                    workflow_state.set_classification(classification_data)
                    workflow_state.complete = True
                    logger.debug(
                        f"Contained match found - stage: {workflow_state.match_stage}, result: {classification_data}"
                    )
                    return workflow_state

            if stage_id == "similar":
                logger.debug("Running similarity match check")
                result, similarities = check_similar_match(
                    fasta=blast_data.fasta_file,
                    existing_ships=ships_df,
                    threshold=0.95,
                )

                if result:
                    logger.debug(f"Found similar match: {result}")
                    workflow_state.stages[stage_id]["progress"] = 50
                    workflow_state.stages[stage_id]["status"] = "complete"
                    
                    classification_data.source = "similar"
                    classification_data.closest_match = result
                    classification_data.confidence = "High"
                    classification_data.match_details = f"High similarity to {result}"
                    
                    workflow_state.set_classification(classification_data)
                    workflow_state.complete = True
                    # Don't store similarities in workflow_state to avoid serialization issues
                    logger.debug(
                        f"Similar match found - stage: {workflow_state.match_stage}, result: {classification_data}, similarities: {len(similarities) if similarities else 0} entries"
                    )
                    return workflow_state

            if stage_id == "family":
                logger.debug("Running family classification")

                family_dict, protein_file = classify_family(
                    fasta=blast_data.fasta_file,
                    seq_type=blast_data.seq_type,
                    meta_dict=meta_dict,
                    pident_thresh=90,
                    input_eval=0.001,
                    threads=1,
                )

                if family_dict:
                    family_name = family_dict["family"]
                    logger.debug(f"Found family classification: {family_name}")
                    workflow_state.stages[stage_id]["progress"] = 70
                    workflow_state.stages[stage_id]["status"] = "complete"
                    workflow_state.complete = True
                    workflow_state.found_match = True
                    workflow_state.match_stage = "family"

                    # Simplify the result before setting it in the workflow state
                    if (
                        isinstance(family_dict, tuple)
                        and len(family_dict) > 0
                        and isinstance(family_dict[0], dict)
                        and "family" in family_dict[0]
                    ):
                        family=family_dict[0]["family"]
                    elif isinstance(family_dict, dict) and "family" in family_dict:
                        family = family_dict["family"]
                    else:
                        family = family_dict   

                    classification_data.source = "family"
                    classification_data.family = family
                    classification_data.confidence = "Medium"

                    # Check hits count
                    n_hits = family_dict.get("n_hits", 0)
                    if n_hits == 0:
                        logger.debug("Family matched but with 0 hits (unusual state)")

                    workflow_state.set_classification(classification_data)
                    logger.debug(
                        f"Family classification result: {classification_data}"
                    )
                    return workflow_state
                else:
                    # Handle the case where no HMMER hits were found
                    logger.debug(
                        "No family classification found (HMMER returned no hits)"
                    )
                    workflow_state.stages[stage_id]["progress"] = 100
                    workflow_state.stages[stage_id]["status"] = "complete"
                    # Mark this as a meaningful "no match" result rather than incomplete
                    workflow_state.found_match = False
                    workflow_state.match_stage = "family"
                    workflow_state.error = "No hits found"
                    workflow_state.complete = True
                    return workflow_state

            if stage_id == "navis":
                logger.debug("Running navis classification")
                if captains_df.empty:
                    logger.warning("No captain data available for navis classification")
                    workflow_state.stages[stage_id]["progress"] = 80
                    workflow_state.stages[stage_id]["status"] = "skipped"
                else:
                    result = classify_navis(
                        fasta=blast_data.fasta_file,
                        existing_captains=captains_df,
                        threads=1,
                    )

                if result:
                    logger.debug(f"Found navis classification: {result}")
                    workflow_state.stages[stage_id]["progress"] = 90
                    workflow_state.stages[stage_id]["status"] = "complete"
                    workflow_state.found_match = True
                    workflow_state.match_stage = "navis"
                    classification_data.source = "navis"
                    classification_data.navis = result
                    classification_data.confidence = "Medium"
                    workflow_state.set_classification(classification_data)
                    workflow_state.complete = True
                    return workflow_state

            if stage_id == "haplotype":
                logger.debug("Running haplotype classification")
                if captains_df.empty or ships_df.empty:
                    logger.warning("Missing data for haplotype classification")
                    workflow_state.stages[stage_id]["progress"] = 90
                    workflow_state.stages[stage_id]["status"] = "skipped"
                else:
                    # Extract navis value from the first captain record
                    navis_value = None
                    try:
                        if (
                            not captains_df.empty
                            and "navis_name" in captains_df.columns
                        ):
                            navis_values = captains_df["navis_name"].dropna().unique()
                            if len(navis_values) > 0:
                                navis_value = navis_values[0]
                                logger.debug(f"Using navis value: {navis_value}")
                            else:
                                logger.warning("No non-null navis values found")
                        else:
                            logger.warning("No navis column in captains data")

                        if navis_value is None and "navis_name" in ships_df.columns:
                            navis_values = ships_df["navis_name"].dropna().unique()
                            if len(navis_values) > 0:
                                navis_value = navis_values[0]
                                logger.debug(
                                    f"Using navis value from ships_df: {navis_value}"
                                )
                    except Exception as e:
                        logger.error(f"Error extracting navis value: {e}")

                if navis_value is None:
                    logger.warning("No navis value found, using fallback value 'UNK'")
                    navis_value = "UNK"

                try:
                    result = classify_haplotype(
                        fasta=blast_data.fasta_file,
                        existing_ships=ships_df,
                        navis=navis_value,
                        similarities=similarities,
                    )

                    if result:
                        logger.debug(f"Found haplotype classification: {result}")
                        workflow_state.stages[stage_id]["progress"] = 100
                        workflow_state.stages[stage_id]["status"] = "complete"
                        workflow_state.complete = True
                        workflow_state.found_match = True
                        workflow_state.match_stage = "haplotype"
                        classification_data.source = "haplotype"
                        classification_data.haplotype = result
                        classification_data.confidence = "High"
                        workflow_state.set_classification(classification_data)
                        return workflow_state
                except Exception as e:
                    logger.error(f"Error in haplotype classification: {e}")
                    workflow_state.stages[stage_id]["status"] = "error"
                    workflow_state.error = (
                        f"Haplotype classification error: {str(e)}"
                    )

            # Mark this stage as complete
            workflow_state.stages[stage_id]["progress"] = 100
            workflow_state.stages[stage_id]["status"] = "complete"

        # If no classification found through our methods, try BLAST results as a final fallback
        if not workflow_state.found_match == False and blast_data.blast_df:
            logger.debug(
                "No classification found through workflow methods, trying BLAST fallback"
            )
            try:
                # Convert blast_df back to DataFrame if it's a list of records
                if isinstance(blast_data.blast_df, list):
                    blast_df = pd.DataFrame(blast_data.blast_df)
                else:
                    blast_df = blast_data.blast_df

                if not blast_df.empty:
                    # Sort by evalue (ascending) and pident (descending) to get best hits
                    blast_df = blast_df.sort_values(
                        ["evalue", "pident"], ascending=[True, False]
                    )
                    top_hit = blast_df.iloc[0]

                    top_evalue = float(top_hit["evalue"])
                    top_aln_length = int(top_hit["aln_length"])
                    top_pident = float(top_hit["pident"])

                    if top_pident >= 90:
                        hit_IDs = top_hit["hit_IDs"]
                        # Convert hit_IDs to a list if it's a string
                        hit_IDs_list = (
                            [hit_IDs] if isinstance(hit_IDs, str) else hit_IDs
                        )

                        # Look up metadata for this hit
                        if meta_dict:
                            meta_df = pd.DataFrame(meta_dict)
                            meta_df_sub = meta_df[
                                meta_df["accession_tag"].isin(hit_IDs_list)
                            ]

                            if not meta_df_sub.empty:
                                top_family = meta_df_sub["familyName"].iloc[0]

                                workflow_state.found_match = True
                                workflow_state.match_stage = "blast_hit"
                                classification_data.source = "blast_hit"
                                classification_data.family = top_family
                                classification_data.closest_match=hit_IDs
                                classification_data.match_details=f"BLAST hit with {top_pident:.1f}% identity (length {top_aln_length}bp, E-value: {top_evalue})"
                                classification_data.confidence="High" if top_pident >= 90 and top_aln_length > 1000 else "Medium" if top_pident >= 70 else "Low"
                                
                                workflow_state.set_classification(classification_data)
                                logger.debug(
                                    f"Found BLAST-based classification: {top_family}"
                                )
            except Exception as e:
                logger.error(f"Error processing BLAST fallback: {e}")

        # Mark workflow as complete even if no matches were found
        workflow_state.complete = True
        return workflow_state

    except Exception as e:
        error_message = str(e)
        logger.error(f"Error in classification workflow: {error_message}")
        logger.exception("Full traceback:")

        workflow_state.error = error_message
        workflow_state.complete = True
        return workflow_state


def create_classification_card(classification_data):
    """
    Create a card displaying classification results from multiple sources

    Args:
        classification_data: Dictionary containing classification information

    Returns:
        dmc.Paper component with classification information
    """
    from dash_iconify import DashIconify
    import dash_mantine_components as dmc
    from src.config.logging import get_logger
    
    logger = get_logger(__name__)

    if not classification_data:
        return None

    # Debug: Log the classification data being passed to the card
    logger.debug(f"create_classification_card received data: {classification_data}")

    source = classification_data.get("source")
    family = classification_data.get("family")
    navis = classification_data.get("navis")
    haplotype = classification_data.get("haplotype")
    # TODO: update this to be the accession_display (accession_tag and version_tag)
    closest_match = classification_data.get("closest_match")
    match_details = classification_data.get("match_details")
    confidence = classification_data.get("confidence")
    
    logger.debug(f"Extracted values - source: {source}, family: {family}, navis: {navis}, haplotype: {haplotype}, confidence: {confidence}")
    
    # Check if this is an empty classification (all important fields are None)
    if not any([source, family, navis, haplotype, closest_match, confidence]):
        logger.debug("Classification data is empty (all fields are None), returning None")
        return None

    source_color = source_colors.get(source, "gray")
    confidence_color = confidence_colors.get(confidence, "gray")
    confidence_icon = confidence_icons.get(confidence, "mdi:shield-outline")

    # Create list of classification details
    details = []

    if confidence:
        details.append(
            dmc.Group(
                [
                    dmc.Text(f"Confidence: {confidence}", fw=700, size="lg"),
                    dmc.ThemeIcon(
                        DashIconify(icon=confidence_icon, width=16),
                        size="md",
                        variant="light",
                        color=confidence_color,
                    ),
                ],
                pos="right",
            ),
        )

    if match_details:
        details.append(
            dmc.Group(
                [
                    dmc.Badge(source.replace("_", " ").title(), color=source_color),
                    dmc.Text(match_details, c="dimmed", size="sm"),
                ],
                pos="apart",
            )
        )

    if closest_match and confidence != "Low":
        if source == "exact":
            closest_match_text = f"Exact match to {closest_match}"
        elif source == "contained":
            # Check if this is a perfect match (High confidence + perfect match in details)
            if (confidence == "High" and match_details and
                ("perfect match" in match_details.lower() or "perfect" in match_details.lower())):
                closest_match_text = f"Perfect match to {closest_match}"
            else:
                closest_match_text = f"Contained in {closest_match}"
        elif source == "similar":
            closest_match_text = f"Similar to {closest_match}"
        else:
            closest_match_text = f"Closest match: {closest_match}"
        details.append(
            dmc.Group(
                [
                    dmc.Text(closest_match_text, fw=700, size="lg"),
                ],
                pos="apart",
            )
        )

    if family:
        if isinstance(family, dict) and "family" in family:
            family = family["family"]

        details.append(
            dmc.Group(
                [
                    dmc.Badge(source.replace("_", " ").title(), color=source_color)
                    if source
                    not in ["exact", "contained", "similar", "blast_hit", "hmmsearch"]
                    else None,
                    dmc.Text("Family:", fw=700),
                    dmc.Text(family),
                ],
                pos="apart",
            )
        )

    if navis and confidence != "Low":
        details.append(
            dmc.Group(
                [
                    dmc.Badge(source.replace("_", " ").title(), color=source_color)
                    if source
                    not in ["exact", "contained", "similar", "blast_hit", "hmmsearch"]
                    else None,
                    dmc.Text("Navis:", fw=700),
                    dmc.Text(navis),
                ],
                pos="apart",
            )
        )

    if haplotype and confidence == "High":
        details.append(
            dmc.Group(
                [
                    dmc.Badge(source.replace("_", " ").title(), color=source_color)
                    if source
                    not in ["exact", "contained", "similar", "blast_hit", "hmmsearch"]
                    else None,
                    dmc.Text("Haplotype:", fw=700),
                    dmc.Text(haplotype),
                ],
                pos="apart",
            )
        )

    # Create card
    return dmc.Paper(
        children=[
            *details,
        ],
        p="md",
        withBorder=True,
        shadow="sm",
        radius="md",
        style={
            "marginBottom": "1rem",
            "width": "fit-content",
            "maxWidth": "100%",
            "display": "inline-block",
        },
    )


def create_classification_output(workflow_state=None, classification_data=None):
    """Create the classification output component"""
    import dash_html_components as html
    import dash_mantine_components as dmc

    # Extract classification data
    classification_title = dmc.Title(
        "Classification Results",
        order=2,
        style={"marginTop": "15px", "marginBottom": "20px"},
    )

    if classification_data:
        classification_card = create_classification_card(classification_data)
        if classification_card:
            return html.Div(
                [
                    classification_title,
                    classification_card,
                ]
            )
        else:
            # Classification data exists but is empty - treat as no classification
            classification_data = None
    
    # No valid classification data - check workflow state
    # Handle workflow_state as either object or dictionary
    workflow_complete = False
    if workflow_state:
        if hasattr(workflow_state, 'complete'):
            workflow_complete = workflow_state.complete
        elif isinstance(workflow_state, dict):
            workflow_complete = workflow_state.get('complete', False)
        
    # Check if workflow is still running
    if workflow_state and not workflow_complete:
        # Calculate progress from workflow state
        progress = 0
        current_stage_text = "Starting classification..."

        # Get current stage index
        current_stage_idx = None
        if hasattr(workflow_state, 'current_stage_idx'):
            current_stage_idx = workflow_state.current_stage_idx
        elif isinstance(workflow_state, dict):
            current_stage_idx = workflow_state.get('current_stage_idx')
                
        if current_stage_idx is not None:
            try:
                stage_idx = int(current_stage_idx)
                total_stages = 6  # Number of workflow stages

                # Get current stage progress
                current_stage = None
                stages_dict = None
                if hasattr(workflow_state, 'current_stage'):
                    current_stage = workflow_state.current_stage
                    stages_dict = workflow_state.stages
                elif isinstance(workflow_state, dict):
                    current_stage = workflow_state.get('current_stage')
                    stages_dict = workflow_state.get('stages', {})

                if current_stage and current_stage in stages_dict:
                    stage_data = stages_dict[current_stage]
                    stage_progress = (
                        stage_data.get("progress", 0)
                        if isinstance(stage_data, dict)
                        else 0
                    )
                else:
                    stage_progress = 0

                # Calculate overall progress
                progress = int(
                    (stage_idx / total_stages) * 100
                    + (stage_progress / total_stages)
                )
                progress = max(0, min(100, progress))
                # get stage label from WORKFLOW_STAGES
                stage_labels = {
                    stage["id"]: stage["label"] for stage in WORKFLOW_STAGES
                }
                current_stage_text = stage_labels.get(
                    current_stage,
                    f"Processing {current_stage}"
                    if current_stage
                    else "Running classification...",
                )

            except (ValueError, ZeroDivisionError, TypeError):
                progress = 0

        # Handle case where workflow started but current_stage info is missing
        if not current_stage_text or current_stage_text == "Processing None":
            current_stage_text = "Running comprehensive classification..."

        # Workflow is still running - show progress bar and loader
        return html.Div(
            [
                classification_title,
                dmc.Stack(
                    [
                        dmc.Group(
                                [
                                    dmc.Loader(size="sm", color="blue"),
                                    dmc.Text(
                                        "Classification In Progress",
                                        size="lg",
                                        fw=500,
                                        c="blue",
                                    ),
                                ],
                                gap="md",
                                align="center",
                            ),
                            dmc.Progress(
                                value=progress,
                                color="blue",
                                size="lg",
                                animated=True,
                                striped=True,
                                style={"width": "100%"},
                            ),
                            dmc.Text(
                                current_stage_text, size="sm", c="dimmed", ta="center"
                            ),
                        ],
                        gap="sm",
                    ),
                ]
            )
    else:
        # Workflow is complete but no classification available
        return html.Div(
            [
                classification_title,
                dmc.Alert(
                    title="No Classification Available",
                    children="Could not classify this sequence with any available method.",
                    color="yellow",
                    variant="light",
                ),
                ]
            )

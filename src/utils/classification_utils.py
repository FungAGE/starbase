import tempfile
import subprocess
import pandas as pd
import hashlib
import os
import glob
import signal
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
    generate_md5_hash,
)
from src.database.sql_manager import fetch_ships
from Bio import SeqIO
import networkx as nx

from typing import Optional, Tuple, Dict, Any
from src.database.sql_manager import fetch_ships, fetch_captains
from src.utils.blast_data import WorkflowState

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
    sequence: str, existing_ships: pd.DataFrame = None, threshold: float = 0.95
) -> Tuple[str, bool]:
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

    logger.debug("Starting accession assignment process")

    logger.debug("Step 1: Checking for exact matches using MD5 hash...")
    exact_match = check_exact_match(sequence, existing_ships)
    if exact_match:
        logger.debug(f"Found exact match: {exact_match}")
        return exact_match, False

    logger.debug("Step 2: Checking for contained matches...")
    container_match = check_contained_match(
        fasta=sequence,
        existing_ships=existing_ships,
        min_coverage=0.95,
        min_identity=0.95,
    )
    if container_match:
        logger.debug(f"Found containing match: {container_match}")
        return container_match, True  # Flag for review since it's truncated

    logger.debug(f"Step 3: Checking for similar matches (threshold={threshold})...")
    similar_match = check_similar_match(sequence, existing_ships, threshold)
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
        int(acc.replace("SBS", "").split(".")[0])
        for acc in existing_ships["accession_tag"]
        if acc.startswith("SBS")
    ]

    # Check if we have existing accessions
    if not existing_nums:
        error_msg = "Problem with loading existing ships. No existing SBS accessions found in database."
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Find next available number
    next_num = max(existing_nums) + 1
    logger.debug(f"Last used accession number: SBS{max(existing_nums):06d}")
    logger.debug(f"Assigning new accession number: SBS{next_num:06d}")

    return f"SBS{next_num:06d}"


########################################################
# sequence matching functions
########################################################


def generate_md5_hash(sequence: str) -> str:
    """Generate an MD5 hash of a sequence."""
    if sequence is None:
        logger.error("Cannot generate MD5 hash for None sequence")
        return None
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
        sequence = list(sequences.values())[0]

    # Normalize sequence by removing whitespace and making uppercase
    clean_sequence = "".join(sequence.upper().split())
    sequence_hash = generate_md5_hash(clean_sequence)
    logger.debug(
        f"Query sequence hash: {sequence_hash} (sequence length: {len(clean_sequence)})"
    )

    # Calculate hashes for existing sequences, skipping None values
    existing_hashes = {}
    skipped_count = 0
    for acc, seq in zip(existing_ships["accession_tag"], existing_ships["sequence"]):
        if seq is None:
            skipped_count += 1
            logger.warning(f"Skipping null sequence for accession {acc}")
            continue

        # Normalize sequence the same way
        clean_db_seq = "".join(seq.upper().split())
        existing_hashes[generate_md5_hash(clean_db_seq)] = acc

    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} sequences due to null values")

    logger.debug(f"Compared against {len(existing_hashes)} valid sequences")
    match = existing_hashes.get(sequence_hash)
    if match:
        logger.debug(f"Found exact hash match: {match}")
    else:
        # If fasta is not a file path, treat it as a sequence string
        logger.debug(f"Treating input as sequence string, length: {len(fasta)}")
        sequence_list = [fasta]

    if not sequence_list:
        logger.error("No sequence provided")
        return None

    # generate md5 for query sequence(s)
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
        if md5_hash is None:
            logger.warning(
                f"generate_md5_hash returned None for sequence: {seq[:50]}..."
            )
            continue
        query_md5[seq] = md5_hash

    # collect existing md5 in dict
    # Create reverse mapping: md5 -> accession_display (includes version)
    existing_hashes = {}
    sequences_without_md5 = []

    for _, row in existing_ships.iterrows():
        # Use accession_display if available, otherwise fall back to accession_tag
        accession_display = row.get("accession_display", row.get("accession_tag"))

        if row.get("md5") is not None and accession_display is not None:
            existing_hashes[row["md5"]] = accession_display
        elif row.get("sequence") is not None and accession_display is not None:
            # Calculate MD5 on the fly for sequences without stored MD5
            clean_seq = clean_sequence(row["sequence"])
            if clean_seq:
                calculated_md5 = generate_md5_hash(clean_seq)
                if calculated_md5:
                    existing_hashes[calculated_md5] = accession_display
                    sequences_without_md5.append(accession_display)

    logger.debug(f"Found {len(existing_hashes)} existing MD5 hashes in database")
    if sequences_without_md5:
        logger.debug(
            f"Calculated MD5 for {len(sequences_without_md5)} sequences on-the-fly"
        )
    logger.debug(f"Query MD5 hashes: {list(query_md5.values())}")

    # Check if query hash exists in database
    for seq, md5 in query_md5.items():
        logger.debug(f"Checking query MD5: {md5}")
        match = existing_hashes.get(md5)
        if match:
            logger.debug(f"Found exact hash match: {match}")
            return match

    logger.debug("No exact match found")
    return None


def check_contained_match(
    fasta: str,
    existing_ships: pd.DataFrame,
    min_coverage: float = 0.95,
    min_identity: float = 0.95,
) -> Optional[str]:
    """Check if sequence is contained within any existing sequences.

    Args:
        sequence: Query sequence to check
        existing_ships: DataFrame containing existing sequences
        min_coverage: Minimum coverage of query sequence required (default: 0.95)
        min_identity: Minimum sequence identity required (default: 0.95)

    Returns:
        accession_tag of the best containing match, or None if no match found
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
                        containing_matches.append(
                            (
                                identity * coverage,  # score for sorting
                                sequence_length,  # length for tiebreaking
                                ref_name,
                            )
                        )

            logger.debug(f"Processed {alignment_count} alignments from minimap2")

    # Sort by score descending, then by length descending
    if containing_matches:
        containing_matches.sort(reverse=True)
        logger.debug(f"Found {len(containing_matches)} containing matches")
        logger.debug(
            f"Selected best match: {containing_matches[0][2]} "
            f"(score: {containing_matches[0][0]:.2f}, length: {containing_matches[0][1]})"
        )
        return containing_matches[0][2]  # Return accession of best match

    logger.debug("No containing matches found")
    return None


def check_similar_match(
    fasta: str, existing_ships: pd.DataFrame, threshold: float
) -> Tuple[Optional[str], Any]:
    """Check for sequences with similarity above threshold using k-mer comparison."""
    logger.debug(f"Starting similarity comparison (threshold={threshold})")

    tmp_fasta_all_ships = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_multi_fasta(
        existing_ships,
        tmp_fasta_all_ships,
        sequence_col="sequence",
        id_col="accession_display"
        if "accession_display" in existing_ships.columns
        else "accession_tag",
    )

    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name

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

    # Create temporary FASTA with new and existing sequences
    write_combined_fasta(
        sequence,
        existing_ships,
        fasta_path=tmp_fasta,
        sequence_col="sequence",
        id_col="accession_display"
        if "accession_display" in existing_ships.columns
        else "accession_tag",
    )
    logger.debug(f"Created temporary FASTA file: {tmp_fasta}")

    # Calculate similarities
    similarities = calculate_similarities(
        fasta_file=tmp_fasta,
        seq_type="nucl",
    )

    # Convert similarities dictionary to a list of triplets for processing
    similarity_triplets = []
    query_id = "query_sequence"  # This is the ID used in write_combined_fasta

    # Check if similarities is a dictionary (new format) or a list (old format)
    if isinstance(similarities, dict):
        for seq_id1 in similarities:
            for seq_id2, sim in similarities[seq_id1].items():
                if seq_id1 == query_id or seq_id2 == query_id:
                    # The other ID is the match
                    match_id = seq_id2 if seq_id1 == query_id else seq_id1
                    logger.debug(f"Similarity to {match_id}: {sim}")
                    if sim >= threshold:
                        logger.debug(
                            f"Found similar match: {match_id} (similarity: {sim})"
                        )
                        return match_id, similarities
    else:
        # Handle the case where similarities is already a list of triplets
        for similarity_tuple in similarities:
            if len(similarity_tuple) == 3:
                seq_id1, seq_id2, sim = similarity_tuple
                if seq_id1 == query_id or seq_id2 == query_id:
                    match_id = seq_id2 if seq_id1 == query_id else seq_id1
                    logger.debug(f"Similarity to {match_id}: {sim}")
                    if sim >= threshold:
                        logger.debug(
                            f"Found similar match: {match_id} (similarity: {sim})"
                        )
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

    protein = list(load_fasta_to_dict(fasta).values())[0]

    logger.debug("Starting navis classification")
    tmp_fasta_dir = create_tmp_fasta_dir(protein, existing_captains)
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
        # Try both captainID and accession_display if available
        matching_captains = existing_captains[
            existing_captains["captainID"] == cluster_rep
        ]

        # If not found and accession_display column exists, try that too
        if matching_captains.empty and "accession_display" in existing_captains.columns:
            matching_captains = existing_captains[
                existing_captains["accession_display"] == cluster_rep
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

        logger.debug(
            f"existing_ships shape: {existing_ships.shape if not existing_ships.empty else 'empty'}"
        )
        if not existing_ships.empty:
            logger.debug(f"existing_ships columns: {existing_ships.columns.tolist()}")

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

            # Cluster sequences
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
            ship_haplotypes = filtered_ships[
                filtered_ships["captainID"].isin(ship_ids)
            ]["haplotype_name"].dropna()

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
        # Handle None similarities
        if similarities is None:
            logger.warning("No similarities provided (None)")
            return []

        # Handle empty similarities
        if not similarities:
            logger.warning("No similarities provided (empty list)")
            return []

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


# TODO: make sure `ref_db` is a fasta file (do we need to `createdb` for this fasta file?)
# TODO: make sure that `output_prefix` is a temp directory
def metaeuk_easy_predict(query_fasta, ref_db, output_prefix, threads=20):
    """Run MetaEuk easy-predict for de novo annotation.

    Args:
        query_fasta: Path to input FASTA file
        ref_db: Path to reference database
        output_prefix: Prefix for output files
        threads: Number of threads to use
    """

    try:
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Run MetaEuk with explicit paths
            cmd = [
                "metaeuk",
                "easy-predict",
                os.path.abspath(query_fasta),
                os.path.abspath(ref_db),
                os.path.abspath(output_prefix),
                tmp_dir,
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
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Check if output files were created
            codon_fasta = f"{output_prefix}.codon.fas"
            fasta = f"{output_prefix}.fas"
            gff = f"{output_prefix}.gff"

            for file_path in [codon_fasta, fasta, gff]:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(
                        f"Expected output file not created: {file_path}"
                    )

            return codon_fasta, fasta, gff

    except subprocess.CalledProcessError as e:
        logger.error(f"MetaEuk easy-predict failed with return code {e.returncode}")
        logger.error(f"stdout: {e.stdout}")
        logger.error(f"stderr: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"Error during MetaEuk easy-predict: {str(e)}")
        logger.exception("Full traceback:")
        raise


def run_classification_workflow(upload_data, meta_dict=None):
    """Run the classification workflow and return results."""
    import pandas as pd

    # Initialize workflow state as dictionary for Dash compatibility
    workflow_state = {
        "complete": False,
        "error": None,
        "found_match": False,
        "match_stage": None,
        "match_result": None,
        "stages": {
            stage["id"]: {"progress": 0, "status": "pending"}
            for stage in WORKFLOW_STAGES
        },
        "task_id": "",
        "status": "initialized",
        "workflow_started": True,
        "current_stage": None,
        "current_stage_idx": 0,
        "start_time": 0.0,
        "class_dict": {},
    }

    # Initialize similarities to None to avoid NameError in haplotype stage
    similarities = None

    try:
        # Parse parameters for database fetches
        ships_df = fetch_ships(
            curated=upload_data.fetch_ship_params.curated,
            with_sequence=upload_data.fetch_ship_params.with_sequence,
            dereplicate=upload_data.fetch_ship_params.dereplicate,
        )
        captains_df = fetch_captains(
            curated=upload_data.fetch_captain_params.curated,
            with_sequence=upload_data.fetch_captain_params.with_sequence,
        )

        fasta_path = upload_data.fasta_file
        if isinstance(fasta_path, dict) and "content" in fasta_path:
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta_path["content"])
            upload_data.fasta_file = tmp_file
            fasta_path = tmp_file

        # Make sure blast_df is serializable
        if upload_data.blast_df and isinstance(upload_data.blast_df, pd.DataFrame):
            upload_data.blast_df = upload_data.blast_df.to_dict("records")

        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]

            # Update state to show progress
            workflow_state["current_stage"] = stage_id
            workflow_state["current_stage_idx"] = i
            workflow_state["stages"][stage_id]["progress"] = 10
            workflow_state["stages"][stage_id]["status"] = "running"

            logger.debug(f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}")

            if stage_id == "exact":
                logger.debug("Running exact match check")

                result = check_exact_match(
                    fasta=upload_data.fasta_file, existing_ships=ships_df
                )

                if result:
                    logger.debug(f"Found exact match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "exact"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    logger.debug(
                        f"Exact match found - stage: {workflow_state['match_stage']}, result: {workflow_state['match_result']}"
                    )
                    return workflow_state

            if stage_id == "contained":
                logger.debug("Running contained match check")
                result = check_contained_match(
                    fasta=upload_data.fasta_file,
                    existing_ships=ships_df,
                    min_coverage=0.95,
                    min_identity=0.95,
                )

                if result:
                    logger.debug(f"Found contained match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 30
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "contained"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    logger.debug(
                        f"Contained match found - stage: {workflow_state['match_stage']}, result: {workflow_state['match_result']}"
                    )
                    return workflow_state

            if stage_id == "similar":
                logger.debug("Running similarity match check")
                result, similarities = check_similar_match(
                    fasta=upload_data.fasta_file,
                    existing_ships=ships_df,
                    threshold=0.95,
                )

                if result:
                    logger.debug(f"Found similar match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 50
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "similar"
                    workflow_state["match_result"] = result
                    # Don't store similarities in workflow_state to avoid serialization issues
                    logger.debug(
                        f"Similar match found - stage: {workflow_state['match_stage']}, result: {workflow_state['match_result']}, similarities: {len(similarities) if similarities else 0} entries"
                    )
                    return workflow_state

            if stage_id == "family":
                logger.debug("Running family classification")

                family_dict, protein_file = classify_family(
                    fasta=upload_data.fasta_file,
                    seq_type=upload_data.seq_type,
                    meta_dict=meta_dict,
                    pident_thresh=90,
                    input_eval=0.001,
                    threads=1,
                )

                if family_dict:
                    family_name = family_dict["family"]
                    logger.debug(f"Found family classification: {family_name}")
                    workflow_state["stages"][stage_id]["progress"] = 70
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["complete"] = True
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "family"

                    # Simplify the result before setting it in the workflow state
                    if (
                        isinstance(family_dict, tuple)
                        and len(family_dict) > 0
                        and isinstance(family_dict[0], dict)
                        and "family" in family_dict[0]
                    ):
                        workflow_state["match_result"] = family_dict[0]["family"]
                    elif isinstance(family_dict, dict) and "family" in family_dict:
                        workflow_state["match_result"] = family_dict["family"]
                    else:
                        workflow_state["match_result"] = family_dict

                    # Check hits count
                    n_hits = family_dict.get("n_hits", 0)
                    if n_hits == 0:
                        logger.debug("Family matched but with 0 hits (unusual state)")

                    logger.debug(
                        f"Family classification result: {workflow_state['match_result']}"
                    )
                    return workflow_state
                else:
                    # Handle the case where no HMMER hits were found
                    logger.debug(
                        "No family classification found (HMMER returned no hits)"
                    )
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    # Mark this as a meaningful "no match" result rather than incomplete
                    workflow_state["found_match"] = False
                    workflow_state["match_stage"] = "family"
                    workflow_state["match_result"] = "No hits found"
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "navis":
                logger.debug("Running navis classification")
                if captains_df.empty:
                    logger.warning("No captain data available for navis classification")
                    workflow_state["stages"][stage_id]["progress"] = 80
                    workflow_state["stages"][stage_id]["status"] = "skipped"
                else:
                    result = classify_navis(
                        fasta=upload_data.fasta_file,
                        existing_captains=captains_df,
                        threads=1,
                    )

                if result:
                    logger.debug(f"Found navis classification: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 90
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "navis"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "haplotype":
                logger.debug("Running haplotype classification")
                if captains_df.empty or ships_df.empty:
                    logger.warning("Missing data for haplotype classification")
                    workflow_state["stages"][stage_id]["progress"] = 90
                    workflow_state["stages"][stage_id]["status"] = "skipped"
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
                        fasta=upload_data.fasta_file,
                        existing_ships=ships_df,
                        navis=navis_value,
                        similarities=similarities,
                    )

                    if result:
                        logger.debug(f"Found haplotype classification: {result}")
                        workflow_state["stages"][stage_id]["progress"] = 100
                        workflow_state["stages"][stage_id]["status"] = "complete"
                        workflow_state["complete"] = True
                        workflow_state["found_match"] = True
                        workflow_state["match_stage"] = "haplotype"
                        workflow_state["match_result"] = result
                        return workflow_state
                except Exception as e:
                    logger.error(f"Error in haplotype classification: {e}")
                    workflow_state["stages"][stage_id]["status"] = "error"
                    workflow_state["error"] = (
                        f"Haplotype classification error: {str(e)}"
                    )

            # Mark this stage as complete
            workflow_state["stages"][stage_id]["progress"] = 100
            workflow_state["stages"][stage_id]["status"] = "complete"

        # If no classification found through our methods, try BLAST results as a final fallback
        if not workflow_state.get("found_match", False) and upload_data.blast_df:
            logger.debug(
                "No classification found through workflow methods, trying BLAST fallback"
            )
            try:
                # Convert blast_df back to DataFrame if it's a list of records
                if isinstance(upload_data.blast_df, list):
                    import pandas as pd

                    blast_df = pd.DataFrame(upload_data.blast_df)
                else:
                    blast_df = upload_data.blast_df

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

                                workflow_state["found_match"] = True
                                workflow_state["match_stage"] = "blast_hit"
                                workflow_state["match_result"] = {
                                    "source": "blast_hit",
                                    "family": top_family,
                                    "closest_match": hit_IDs,
                                    "match_details": f"BLAST hit with {top_pident:.1f}% identity (length {top_aln_length}bp, E-value: {top_evalue})",
                                    "confidence": "High"
                                    if top_pident >= 90 and top_aln_length > 1000
                                    else "Medium"
                                    if top_pident >= 70
                                    else "Low",
                                }

                                logger.debug(
                                    f"Found BLAST-based classification: {top_family}"
                                )
            except Exception as e:
                logger.error(f"Error processing BLAST fallback: {e}")

        # Mark workflow as complete even if no matches were found
        workflow_state["complete"] = True
        return workflow_state

    except Exception as e:
        error_message = str(e)
        logger.error(f"Error in classification workflow: {error_message}")
        logger.exception("Full traceback:")

        workflow_state["error"] = error_message
        workflow_state["complete"] = True
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

    if not classification_data:
        return None

    source = classification_data.get("source", "Unknown")
    family = classification_data.get("family")
    navis = classification_data.get("navis")
    haplotype = classification_data.get("haplotype")
    # TODO: update this to be the accession_display (accession_tag and version_tag)
    closest_match = classification_data.get("closest_match")
    match_details = classification_data.get("match_details")
    confidence = classification_data.get("confidence", "Low")

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


def create_classification_output(sequence_results, workflow_state=None):
    """Create the classification output component"""
    import dash_html_components as html
    import dash_mantine_components as dmc

    # Extract classification data
    classification_data = sequence_results.get("classification")
    classification_title = dmc.Title(
        "Classification Results",
        order=2,
        style={"marginTop": "15px", "marginBottom": "20px"},
    )

    if classification_data:
        return html.Div(
            [
                classification_title,
                create_classification_card(classification_data),
            ]
        )
    else:
        # Check if workflow is still running
        if workflow_state and not workflow_state.get("complete", False):
            # Calculate progress from workflow state
            progress = 0
            current_stage_text = "Starting classification..."

            if workflow_state.get("current_stage_idx") is not None:
                try:
                    stage_idx = int(workflow_state.get("current_stage_idx", 0))
                    total_stages = 6  # Number of workflow stages

                    # Get current stage progress
                    current_stage = workflow_state.get("current_stage")
                    stages_dict = workflow_state.get("stages", {})

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

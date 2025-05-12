import tempfile
import subprocess
import logging
import pandas as pd
import hashlib
import os
import glob
import signal
from src.utils.seq_utils import (
    write_combined_fasta,
    write_multi_fasta,
    create_tmp_fasta_dir,
    load_fasta_to_dict,
)
from typing import Optional, Tuple, Dict, Any, Callable
from src.database.sql_manager import fetch_meta_data, fetch_ships, fetch_captains

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

# Define our workflow stages with their colors
WORKFLOW_STAGES = [
    {"id": "exact", "label": "Checking for exact matches", "color": "green"},
    {"id": "contained", "label": "Checking for contained matches", "color": "blue"},
    {"id": "similar", "label": "Checking for similar matches", "color": "violet"},
    # {"id": "denovo", "label": "Running denovo annotation", "color": "orange"},
    {"id": "family", "label": "Running family classification", "color": "green"},
    {"id": "navis", "label": "Running navis classification", "color": "blue"},
    {"id": "haplotype", "label": "Running haplotype classification", "color": "violet"},
]


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
    if existing_ships is None:
        existing_ships = fetch_ships(curated=True)

    logger.info("Starting accession assignment process")

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
    logger.info(f"Last used accession number: SBS{max(existing_nums):06d}")
    logger.info(f"Assigning new accession number: SBS{next_num:06d}")

    return f"SBS{next_num:06d}"


########################################################
# sequence matching functions
########################################################


def check_exact_match(fasta: str, existing_ships: pd.DataFrame) -> Optional[str]:
    """Check if sequence exactly matches an existing sequence using MD5 hash."""
    # fasta should be a file path
    if os.path.exists(fasta):
        logger.info(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None
        sequence = list(sequences.values())[0]

    # Normalize sequence by removing whitespace and making uppercase
    clean_sequence = "".join(sequence.upper().split())
    sequence_hash = hashlib.md5(clean_sequence.encode()).hexdigest()
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
        existing_hashes[hashlib.md5(clean_db_seq.encode()).hexdigest()] = acc

    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} sequences due to null values")

    logger.info(f"Compared against {len(existing_hashes)} valid sequences")
    match = existing_hashes.get(sequence_hash)
    if match:
        logger.info(f"Found exact hash match: {match}")
    else:
        logger.info("No exact match found")
    return match


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
        logger.info(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None
        sequence = list(sequences.values())[0]

    query_len = len(sequence)

    logger.info(f"Checking for contained matches (query length: {query_len})")

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
                    ref_file.write(f">{row['accession_tag']}\n{row['sequence']}\n")
                    ref_count += 1
            ref_file.flush()
            logger.info(f"Written {ref_count} reference sequences for comparison")

            logger.info("Running minimap2 alignment")
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
                    logger.info(
                        f"Found containing match: {ref_name} "
                        f"(coverage: {coverage:.2f}, identity: {identity:.2f})"
                    )
                    containing_matches.append(
                        (
                            identity * coverage,  # score for sorting
                            len(
                                existing_ships[
                                    existing_ships["accession_tag"] == ref_name
                                ]["sequence"].iloc[0]
                            ),  # length for tiebreaking
                            ref_name,
                        )
                    )

            logger.info(f"Processed {alignment_count} alignments from minimap2")

    # Sort by score descending, then by length descending
    if containing_matches:
        containing_matches.sort(reverse=True)
        logger.info(f"Found {len(containing_matches)} containing matches")
        logger.info(
            f"Selected best match: {containing_matches[0][2]} "
            f"(score: {containing_matches[0][0]:.2f}, length: {containing_matches[0][1]})"
        )
        return containing_matches[0][2]  # Return accession of best match

    logger.info("No containing matches found")
    return None


def check_similar_match(
    fasta: str, existing_ships: pd.DataFrame, threshold: float
) -> Optional[str]:
    """Check for sequences with similarity above threshold using k-mer comparison."""
    logger.info(f"Starting similarity comparison (threshold={threshold})")

    tmp_fasta_all_ships = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_multi_fasta(
        existing_ships,
        tmp_fasta_all_ships,
        sequence_col="sequence",
        id_col="accession_tag",
    )

    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name

    if os.path.exists(fasta):
        logger.info(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None
        sequence = list(sequences.values())[0]

    # Create temporary FASTA with new and existing sequences
    write_combined_fasta(
        sequence,
        existing_ships,
        fasta_path=tmp_fasta,
        sequence_col="sequence",
        id_col="accession_tag",
    )
    logger.debug(f"Created temporary FASTA file: {tmp_fasta}")

    # Calculate similarities
    similarities = calculate_similarities(
        ship_sequence_file=tmp_fasta,
        new_sequence_file=tmp_fasta_all_ships,
        seq_type="nucl",
    )

    # return best_match
    logger.debug(f"Calculated similarities for {len(similarities)} sequences")

    # Check if we have any similarities above threshold
    if "query_sequence" in similarities:
        for acc_id, sim in similarities["query_sequence"].items():
            logger.debug(f"Similarity to {acc_id}: {sim}")
            if sim >= threshold:
                logger.info(f"Found similar match: {acc_id} (similarity: {sim})")
                return acc_id

    logger.info("No similar matches found above threshold")
    return None


########################################################
# classification pipeline
########################################################


def classify_family(
    fasta=None,
    seq_type=None,
    blast_df=None,
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

    if os.path.exists(fasta):
        logger.info(f"Loading sequence from file: {fasta}")
        sequences = load_fasta_to_dict(fasta)
        if not sequences:
            logger.error("No sequences found in FASTA file")
            return None

    # Simple selection of top hit using blast or hmmer results
    if blast_df is not None and len(blast_df) > 0:
        # Sort by evalue (ascending) and pident (descending) to get best hits
        blast_df = blast_df.sort_values(["evalue", "pident"], ascending=[True, False])
        top_hit = blast_df.iloc[0]
        logger.info(f"Top hit: {top_hit}")

        top_evalue = float(top_hit["evalue"])
        top_aln_length = int(top_hit["aln_length"])
        top_pident = float(top_hit["pident"])

        if top_pident >= pident_thresh:
            # look up family name from accession tag
            query_accession = top_hit["query_id"]
            top_family = fetch_meta_data(accession_tag=query_accession)["familyName"]
            family_dict = {
                "family": top_family,
                "aln_length": top_aln_length,
                "evalue": top_evalue,
            }

        hmmer_dict, tmp_protein_filename = run_hmmer(
            query_type=seq_type,
            input_gene="tyr",
            input_eval=0.01,
            query_fasta=fasta,
            threads=2,
        )

    if hmmer_dict is not None and len(hmmer_dict) > 0:
        hmmer_df = pd.DataFrame(hmmer_dict)
        if len(hmmer_df) > 0:
            family_name, family_aln_length, family_evalue = select_ship_family(hmmer_df)
            if family_name:
                family_dict = {
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
    """Assign navis based on captain sequence clustering."""
    # Part 2: Navis Assignment
    # - Compare captain sequence to existing classified captains
    # - Use mmseqs clustering to group with existing navis

    # Create temporary dir with FASTAs from:
    # - Captain sequence from new sequence
    # - All existing classified captain sequences

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
        matching_captains = existing_captains[
            existing_captains["captainID"] == cluster_rep
        ]
        if matching_captains.empty:
            logger.warning(
                f"No matching captain found for cluster representative: {cluster_rep}"
            )
            return None

        navis_name = matching_captains["starship_navis"].iloc[0]
        logger.debug(f"Found navis name: {navis_name}")

        return navis_name


def classify_haplotype(fasta: str, existing_ships: pd.DataFrame, navis: str) -> str:
    """Assign haplotype based on full sequence similarity."""
    # Part 3: Haplotype Assignment
    # - Compare full sequence to existing classified sequences
    # - Use sourmash for similarity calculation
    # - Use MCL clustering to group with existing haplotypes

    # Create temporary FASTA with:
    # - New sequence
    # - All existing classified sequences for this navis

    existing_ships_navis = existing_ships[existing_ships["starship_navis"] == navis]
    tmp_fasta = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
    write_combined_fasta(
        fasta,
        existing_ships_navis,
        fasta_path=tmp_fasta,
        sequence_col="sequence",
        id_col="accession_tag",
    )

    # Calculate similarities
    similarities = calculate_similarities(
        ship_sequence_file=tmp_fasta, new_sequence_file=tmp_fasta, seq_type="nucl"
    )

    # Cluster sequences
    groups, _, _ = cluster_sequences(
        similarities, group_prefix="HAP", inflation=1.5, threshold=0.05
    )

    # groups looks like this:
    # 'sequence_id': {
    #     'group_id': 'HAP000001',
    #     'is_representative': True
    # }

    # Find which group has the member "query_sequence"
    rep_group_id = None
    for group, members in groups.items():
        if "query_sequence" in members:
            rep_group_id = group["group_id"]
            break

    # Return haplotype based on existing classifications in that group
    seq_ids_in_rep_group = []
    for group, members in groups.items():
        if group["group_id"] == rep_group_id:
            seq_ids_in_rep_group.append(group["sequence_id"])

    # just get the first sequence id in the group
    haplotype_name = existing_ships[
        existing_ships["captainID"] == seq_ids_in_rep_group[0]
    ]["starship_haplotype"].values[0]

    if haplotype_name is None:
        logger.warning(
            f"No haplotype name, but sequence clusters with captainID: {seq_ids_in_rep_group[0]}"
        )

    return haplotype_name


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


def sourmash_sketch(sequence_file, sig_file, kmer_size, scaled, sketch_type="dna"):
    sketch_cmd = [
        "sourmash",
        "sketch",
        f"{sketch_type}",
        "--singleton",
        "--output",
        sig_file,
        "-p",
        f"k={kmer_size},scaled={scaled},noabund",
        sequence_file,
    ]

    logger.info(f"Creating sourmash signatures: {' '.join(sketch_cmd)}")
    try:
        subprocess.run(sketch_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Sourmash sketch failed: {e.stderr}")
        raise


def sourmash_compare(
    ship_sketch, new_sketch, matrix_file, kmer_size, sketch_type="dna"
):
    compare_cmd = [
        "sourmash",
        "compare",
        ship_sketch,
        new_sketch,
        f"--{sketch_type}",  # dna or protein
        "--csv",
        matrix_file,
        "-k",
        str(kmer_size),
    ]

    logger.info(f"Calculating pairwise similarities: {' '.join(compare_cmd)}")
    try:
        subprocess.run(compare_cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Sourmash compare failed: {e.stderr}")
        raise


def calculate_similarities(
    ship_sequence_file=None, new_sequence_file=None, seq_type="nucl", threads=1
):
    """Calculate k-mer similarity between sequences using sourmash.

    Args:
        sequence_file (str): Path to input FASTA file
        seq_type (str): Sequence type to compare - 'nucl' or 'prot'
        threads (int): Number of threads to use

    Returns:
        dict: Nested dictionary of pairwise similarities {seq1: {seq2: similarity}}
    """
    # Set sourmash parameters based on sequence type
    if seq_type == "nucl":
        kmer_size = 510  # Default from starfish sig
        scaled = 100
        sketch_type = "dna"
    elif seq_type == "prot":
        kmer_size = 17
        scaled = 20
        sketch_type = "protein"
    else:
        raise ValueError(f"Invalid sequence type: {seq_type}")

    with tempfile.TemporaryDirectory() as tmp_dir:
        # 0. create sketch of ships
        ship_sketch_file = os.path.join(tmp_dir, "ships.sig")
        sourmash_sketch(
            ship_sequence_file, ship_sketch_file, kmer_size, scaled, sketch_type
        )

        # 1. Create signature file
        sig_file = os.path.join(tmp_dir, "sequences.sig")
        sourmash_sketch(new_sequence_file, sig_file, kmer_size, scaled, sketch_type)

        # 2. Calculate pairwise similarities
        matrix_file = os.path.join(tmp_dir, "similarity.csv")
        sourmash_compare(
            ship_sketch_file, sig_file, matrix_file, kmer_size, sketch_type
        )

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
        headers = next(f).strip().split(",")
        seq_ids = headers[1:]  # First column is empty

        # Read similarity values
        for i, line in enumerate(f):
            values = line.strip().split(",")
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
    with open(output_file, "w") as f:
        for seq1 in sorted(similarities):
            for seq2 in sorted(similarities[seq1]):
                # Format: seq1\tseq2\tsimilarity
                sim = similarities[seq1][seq2]
                f.write(f"{seq1}\t{seq2}\t{sim:.3f}\n")


def cluster_sequences(
    similarity_file, group_prefix="HAP", inflation=1.5, threshold=0.0, threads=1
):
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

    logger.info(f"Grouped sequences into {group_count} clusters")
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


# denova annotation using metaeuk easy-predict
# TODO: make sure `ref_db` is a fasta file (do we need to `createdb` for this fasta file?)
# TODO: make sure that `output_prefix` is a temp directory
def metaeuk_easy_predict(query_fasta, ref_db, output_prefix, threads=20):
    """Run MetaEuk easy-predict with proper directory handling.

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
                os.path.abspath(query_fasta),  # Use absolute path
                os.path.abspath(ref_db),  # Use absolute path
                os.path.abspath(output_prefix),  # Use absolute path
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


def run_classification_workflow(upload_data):
    """Run the classification workflow and return results."""
    # Initialize state object
    workflow_state = {
        "complete": False,
        "error": None,
        "found_match": False,
        "match_stage": None,
        "match_result": None,
        "stages": {
            stage["id"]: {"progress": 0, "complete": False} for stage in WORKFLOW_STAGES
        },
    }

    try:
        # Process file content if provided instead of file path
        fasta_path = upload_data["fasta"]
        if isinstance(fasta_path, dict) and "content" in fasta_path:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta_path["content"])
            # Update the upload_data to use the temporary file path
            upload_data["fasta"] = tmp_file
            fasta_path = tmp_file

        # Process each stage in sequence
        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]

            # Update state to show progress
            workflow_state["current_stage"] = stage_id
            workflow_state["current_stage_idx"] = i
            workflow_state["stages"][stage_id]["progress"] = 10
            workflow_state["stages"][stage_id]["status"] = "running"

            logger.info(f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}")

            if stage_id == "exact":
                from src.database.sql_manager import fetch_ships

                ships_df = fetch_ships(**upload_data["fetch_ship_params"])
                result = check_exact_match(
                    fasta=upload_data["fasta"], existing_ships=ships_df
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "exact"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "contained":
                result = check_contained_match(
                    fasta=upload_data["fasta"],
                    existing_ships=fetch_ships(**upload_data["fetch_ship_params"]),
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 30
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "contained"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "similar":
                result = check_similar_match(
                    fasta=upload_data["fasta"],
                    existing_ships=fetch_ships(**upload_data["fetch_ship_params"]),
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 50
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "similar"
                    workflow_state["match_result"] = result

            if stage_id == "family":
                result = classify_family(
                    fasta=upload_data["fasta"],
                    seq_type=upload_data["seq_type"],
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 70
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "navis":
                result = classify_navis(
                    fasta=upload_data["fasta"],
                    existing_captains=fetch_captains(
                        **upload_data["fetch_captain_params"]
                    ).to_dict("records"),
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 90
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["complete"] = True
                    return workflow_state

            if stage_id == "haplotype":
                result = classify_haplotype(
                    fasta=upload_data["fasta"],
                    existing_ships=fetch_ships(
                        **upload_data["fetch_ship_params"]
                    ).to_dict("records"),
                    navis=fetch_captains(**upload_data["fetch_captain_params"]).to_dict(
                        "records"
                    ),
                )

                if result:
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["complete"] = True
                    return workflow_state

            # Mark this stage as complete
            workflow_state["stages"][stage_id]["progress"] = 100
            workflow_state["stages"][stage_id]["status"] = "complete"

        # Mark workflow as complete
        workflow_state["complete"] = True
        return workflow_state

    except Exception as e:
        # Handle errors
        error_message = f"Error during classification: {str(e)}"
        logger.error(error_message)
        logger.exception("Full traceback:")

        workflow_state["error"] = error_message
        workflow_state["complete"] = True
        return workflow_state

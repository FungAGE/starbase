import tempfile
import pysam
import pandas as pd
from collections import defaultdict
import os
import subprocess
from typing import Optional
import logging

logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


def process_minimap2_alignments(sam_file, min_aln_length=1000, flank_size=5000):
    """
    Process minimap2 alignments and extract significant matches in flanking regions

    Parameters:
    -----------
    sam_file : str
        Path to the SAM file
    min_aln_length : int
        Minimum alignment length to consider (default: 1000 bp)
    flank_size : int
        Size of flanking regions to consider at start and end (default: 5000 bp)
    """
    alignments = defaultdict(list)

    # Open SAM/BAM file
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam:
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Filter by alignment length
            if read.query_alignment_length < min_aln_length:
                continue

            # Check if alignment overlaps with flanking regions
            query_start = read.query_alignment_start
            query_end = read.query_alignment_end
            query_length = read.query_length

            # Check if alignment is in first or last flank_size bases
            in_first_flank = query_start < flank_size
            in_last_flank = query_end > (query_length - flank_size)

            if not (in_first_flank or in_last_flank):
                continue

            # Store which flank this alignment belongs to
            flank_type = "first" if in_first_flank else "last"

            alignments[read.query_name].append(
                {
                    "query_name": read.query_name,
                    "target_name": read.reference_name,
                    "target_start": read.reference_start,
                    "target_end": read.reference_end,
                    "query_length": query_length,
                    "aligned_length": read.query_alignment_length,
                    "is_reverse": read.is_reverse,
                    "flank_type": flank_type,
                    "query_start": query_start,
                    "query_end": query_end,
                }
            )

    # Convert to DataFrame
    results = []
    for query, hits in alignments.items():
        results.extend(hits)

    df = pd.DataFrame(results)
    if not df.empty:
        logger.info(
            f"\nFound {len(df)} alignments in flanking regions after filtering:"
        )
        logger.info(f"  Minimum alignment length: {min_aln_length}")
        logger.info(f"  Flank size: {flank_size}bp")

        # Print summary of alignments by query and flank type
        for query in df["query_name"].unique():
            query_hits = df[df["query_name"] == query]
            logger.info(f"\n{query}:")
            for flank in ["first", "last"]:
                flank_hits = query_hits[query_hits["flank_type"] == flank]
                if not flank_hits.empty:
                    logger.info(f"  {flank}:")
                    for _, hit in flank_hits.iterrows():
                        logger.info(
                            f"    Match in {hit['target_name']}: "
                            f"positions {hit['target_start']}-{hit['target_end']} "
                            f"(length: {hit['aligned_length']}bp)"
                        )
    else:
        logger.info("\nNo alignments found meeting the filtering criteria.")

    return df


def analyze_flank_pairs(
    df: pd.DataFrame,
    genome_file: str,
    max_gap: Optional[int] = None,
    flank_size: int = 5000,
) -> pd.DataFrame:
    """
    Analyze pairs of flank hits and full-length matches to find potential complete elements
    Only returns the best match for each query in the genome.
    Collapses hits to the same genomic region, grouping queries that match equally well.

    Best matches are determined by:
    1. Full-length matches take precedence over partial matches
    2. For partial matches, prioritize smallest insert size
    3. For equal insert sizes, prioritize longest total alignment length

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing alignment results
    genome_file : str
        Name of the genome file being processed
    max_gap : int | None
        Maximum allowed gap between flanking hits. If None, uses query length
    flank_size : int
        Size of flanking regions to consider

    Returns:
    --------
    pd.DataFrame
        DataFrame containing paired flank information and full-length matches
    """
    all_matches = []

    # First, find all possible matches as before
    for query in df["query_name"].unique():
        query_hits = df[df["query_name"] == query]
        query_length = query_hits.iloc[0]["query_length"]
        best_match = None
        best_score = float("-inf")

        # Skip if query_length is 0
        if query_length == 0:
            logger.warning(f"Query {query} has length 0, skipping")
            continue

        # If max_gap is None, use query length
        effective_max_gap = max_gap if max_gap is not None else query_length

        # First check for full-length alignments
        for target in query_hits["target_name"].unique():
            target_hits = query_hits[query_hits["target_name"] == target]

            for _, hit in target_hits.iterrows():
                alignment_coverage = hit["aligned_length"] / query_length
                if alignment_coverage > 0.9:
                    # Ensure coordinates are positive and properly ordered
                    first_start = max(0, hit["target_start"])
                    first_end = first_start + flank_size
                    last_start = max(0, hit["target_end"] - flank_size)
                    last_end = hit["target_end"]

                    if hit["is_reverse"]:
                        # Swap for reverse orientation
                        first_start, last_start = last_start, first_start
                        first_end, last_end = last_end, first_end

                    # Score full-length matches higher than any partial match
                    score = (
                        1000000 + hit["aligned_length"]
                    )  # Ensures full-length matches always win
                    match = {
                        "genome_file": genome_file,
                        "query_name": query,
                        "target_name": target,
                        "first_start": first_start,
                        "first_end": first_end,
                        "last_start": last_start,
                        "last_end": last_end,
                        "insert_size": 0,
                        "total_length": hit["aligned_length"],
                        "orientation": "reverse" if hit["is_reverse"] else "forward",
                    }
                    if score > best_score:
                        best_score = score
                        best_match = match

        # Only process partial matches if no full-length match was found
        if best_match is None:
            # Process separate flank hits
            first_hits = query_hits[query_hits["flank_type"] == "first"]
            last_hits = query_hits[query_hits["flank_type"] == "last"]

            # Check each combination of first and last hits
            for _, first in first_hits.iterrows():
                for _, last in last_hits.iterrows():
                    if first["target_name"] == last["target_name"]:
                        # Ensure coordinates are positive
                        first_start = max(0, first["target_start"])
                        first_end = max(0, first["target_end"])
                        last_start = max(0, last["target_start"])
                        last_end = max(0, last["target_end"])

                        # Determine orientation and calculate metrics
                        if first["is_reverse"] == last["is_reverse"]:
                            if first_start < last_start:
                                orientation = (
                                    "forward" if not first["is_reverse"] else "reverse"
                                )
                                insert_size = last_start - first_end
                                total_length = last_end - first_start
                            else:
                                orientation = (
                                    "reverse" if not first["is_reverse"] else "forward"
                                )
                                insert_size = first_start - last_end
                                total_length = first_end - last_start

                            # Only consider if insert size is within acceptable range
                            if insert_size <= effective_max_gap and insert_size >= 0:
                                # Score based on insert size (smaller is better) and total length
                                score = -insert_size + (total_length / 1000)
                                match = {
                                    "genome_file": genome_file,
                                    "query_name": query,
                                    "target_name": first["target_name"],
                                    "first_start": first_start,
                                    "first_end": first_end,
                                    "last_start": last_start,
                                    "last_end": last_end,
                                    "insert_size": insert_size,
                                    "total_length": total_length,
                                    "orientation": orientation,
                                }
                                if score > best_score:
                                    best_score = score
                                    best_match = match

        # Add the best match for this query if one was found
        if best_match is not None:
            all_matches.append(best_match)

    # Convert to DataFrame
    matches_df = pd.DataFrame(all_matches)

    if matches_df.empty:
        return matches_df

    # Group by genomic location
    # Two hits are considered to be to the same region if they overlap significantly
    def group_overlapping_hits(df):
        # Sort by start position
        df = df.sort_values("first_start")
        groups = []
        current_group = [df.iloc[0]]

        for _, row in df.iloc[1:].iterrows():
            last_hit = current_group[-1]

            # Check if this hit overlaps significantly with the current group
            same_region = (
                row["target_name"] == last_hit["target_name"]
                and row["orientation"] == last_hit["orientation"]
                and abs(row["first_start"] - last_hit["first_start"]) < flank_size
                and abs(row["last_end"] - last_hit["last_end"]) < flank_size
            )

            if same_region:
                current_group.append(row)
            else:
                groups.append(current_group)
                current_group = [row]

        groups.append(current_group)
        return groups

    # Process each target sequence separately
    final_hits = []
    for target in matches_df["target_name"].unique():
        target_hits = matches_df[matches_df["target_name"] == target]
        groups = group_overlapping_hits(target_hits)

        for group in groups:
            # Take coordinates from the first hit in the group
            base_hit = group[0]
            queries = [hit["query_name"] for hit in group]

            final_hits.append(
                {
                    "genome_file": genome_file,
                    "query_names": queries,
                    "target_name": base_hit["target_name"],
                    "first_start": base_hit["first_start"],
                    "first_end": base_hit["first_end"],
                    "last_start": base_hit["last_start"],
                    "last_end": base_hit["last_end"],
                    "insert_size": base_hit["insert_size"],
                    "total_length": base_hit["total_length"],
                    "orientation": base_hit["orientation"],
                    "num_queries": len(queries),
                }
            )

    return pd.DataFrame(final_hits)


def process_genome(
    ship_fasta: str,
    genome_fasta: str,
    mmi_file: str,
    preset: str,
    flank_size: int = 5000,
    max_gap: Optional[int] = None,
):
    """Process ships against a single genome using minimap2"""
    try:
        sam_file = tempfile.NamedTemporaryFile(delete=False).name
        subprocess.run(
            [
                "minimap2",
                "-t",
                "1",
                "-ax",
                preset,
                "--end-bonus",
                "100",
                "-O",
                "5,56",
                "-E",
                "4,1",
                "-z",
                "400",
                "-B",
                "2",
                mmi_file,
                ship_fasta,
                genome_fasta,
                "-o",
                sam_file,
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        # Process alignments
        alignment_result = process_minimap2_alignments(sam_file, flank_size=flank_size)

        if not alignment_result.empty:
            paired_results = analyze_flank_pairs(
                alignment_result, genome_fasta, max_gap=max_gap
            )
            return paired_results

    finally:
        if os.path.exists(sam_file):
            os.remove(sam_file)

    return None


# TODO: handle creation of ship fasta
# TODO: parallelize genome search in batches
def process_genome_search(
    ship_dict: dict,
    genome_fasta: str,
    flank_size: int = 5000,
    max_gap: Optional[int] = None,
    seq_div: str = "0.1",
    threads: int = 20,
):
    """
    Process a genome against a set of query sequences

    Parameters:
    -----------
    ship_dict : dict
        Dictionary containing ship sequences
    genome_fasta : str
        FASTA file containing a genome sequence
    seq_div : str
        Amount of sequence divergence to tolerate (default: 0.1)
    flank_size : int
        Flank size around alignments (default: 5000)
    max_gap : int | None
        Maximum allowed gap between flanking hits. If None, uses query length for each query
    threads : int
        Number of threads to use (default: 20)
    """

    if seq_div == "0.1":
        preset = "asm5"
    elif seq_div == "1":
        preset = "asm10"
    elif seq_div == "5":
        preset = "asm20"
    else:
        preset = "asm5"

    mmi_file = tempfile.NamedTemporaryFile(delete=False).name

    minimap2_index = subprocess.run(
        [
            "minimap2",
            "-t",
            str(threads),
            "-d",
            mmi_file,
            genome_fasta,
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if minimap2_index.returncode != 0:
        logger.error(f"Minimap2 index failed for {genome_fasta}")
        return None

    paired_results = process_genome(
        ship_dict, genome_fasta, mmi_file, preset, flank_size, max_gap
    )

    return paired_results

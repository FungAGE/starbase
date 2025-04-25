import dash_mantine_components as dmc
from dash import dcc, html

import os
import tempfile
import subprocess
import json
import pandas as pd

from Bio import SearchIO

from src.utils.seq_utils import (
    clean_shipID,
)
from src.database.blastdb import blast_db_exists, create_dbs
from src.components.error_boundary import create_error_alert
import logging

logger = logging.getLogger(__name__)


def run_blast(db_list, query_type, query_fasta, tmp_blast, input_eval=0.01, threads=2):
    try:
        # Add input size check
        max_input_size = 10 * 1024 * 1024  # 10MB
        if os.path.getsize(query_fasta) > max_input_size:
            logger.error(
                f"Input FASTA file too large: {os.path.getsize(query_fasta)} bytes"
            )
            return None

        # db_list should now be a direct path string
        if not isinstance(db_list, (str, bytes, os.PathLike)):
            logger.error(f"db_list must be a path-like object, got {type(db_list)}")
            raise ValueError("Invalid database path")

        logger.debug(f"Using database path: {db_list}")

        if not blast_db_exists(db_list):
            logger.info("BLAST database not found. Creating new database...")
            create_dbs()

            if not blast_db_exists(db_list):
                logger.error("Failed to create BLAST database")
                raise ValueError("Failed to create BLAST database")

        blast_cmd = [
            "blastn" if query_type == "nucl" else "blastp",
            "-query",
            str(query_fasta),
            "-db",
            str(db_list),
            "-out",
            str(tmp_blast),
            "-evalue",
            str(input_eval),
            "-num_threads",
            str(threads),
            "-max_target_seqs",
            "10",
            "-max_hsps",
            "1",
            "-outfmt",
            "5",  # Changed to default BLAST output format
        ]

        subprocess.run(blast_cmd, check=True, timeout=300)

        return tmp_blast

    except subprocess.TimeoutExpired:
        logger.error("BLAST search timed out after 5 minutes")
        return None
    except Exception as e:
        logger.error(f"Error during BLAST search: {e}")
        return None


def hmmsearch(
    db_list=None,
    query_type=None,
    input_gene="tyr",
    input_eval=None,
    query_fasta=None,
    threads=None,
):
    """Run HMMER search."""
    try:
        logger.debug(
            f"Starting hmmsearch with params: query_type={query_type}, input_gene={input_gene}"
        )

        tmp_hmmer = tempfile.NamedTemporaryFile(suffix=".hmmer.txt").name
        logger.debug(f"Created temporary output file: {tmp_hmmer}")

        hmmer_db = db_list["gene"][input_gene]["hmm"][query_type]
        logger.debug(f"Using HMMER database: {hmmer_db}")

        if not os.path.exists(hmmer_db):
            logger.error(f"HMMER database not found: {hmmer_db}")
            raise ValueError(f"HMMER database {hmmer_db} not found.")

        if os.path.getsize(hmmer_db) == 0:
            logger.error(f"HMMER database is empty: {hmmer_db}")
            raise ValueError(f"HMMER database {hmmer_db} is empty.")

        # Calculate optimal thread count
        if threads is None:
            threads = 1
        else:
            threads = min(int(threads), 4)
        logger.debug(f"Using {threads} threads")

        hmmer_cmd = f"hmmsearch -o {tmp_hmmer} --cpu {threads} -E {input_eval} {hmmer_db} {query_fasta}"
        logger.debug(f"Running HMMER command: {hmmer_cmd}")

        result = subprocess.run(hmmer_cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"HMMER search failed with return code {result.returncode}")
            logger.error(f"STDERR: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, hmmer_cmd)

        logger.debug("HMMER search completed successfully")
        return tmp_hmmer

    except Exception as e:
        logger.error(f"Error in hmmsearch: {str(e)}")
        logger.exception("Full traceback:")
        raise


def run_hmmer(
    db_list=None,
    query_type=None,
    input_gene="tyr",
    input_eval=0.01,
    query_fasta=None,
    threads=2,
):
    """Run HMMER search on input sequence.

    Args:
        db_list: Dictionary containing database paths
        query_type: Type of sequence ('nucl' or 'prot')
        input_gene: Gene to search for (default: 'tyr')
        input_eval: E-value threshold
        query_fasta: Path to query FASTA file
        threads: Number of CPU threads to use

    Returns:
        For nucleotide input (query_type=="nucl"):
            Tuple[Dict, str]: (HMMER results dictionary, path to protein sequence file)
        For protein input (query_type=="prot"):
            Tuple[Dict, None]: (HMMER results dictionary, None)
    """
    try:
        logger.debug(
            f"Starting HMMER search with params: query_type={query_type}, input_gene={input_gene}"
        )
        logger.debug(f"Input FASTA file: {query_fasta}")

        protein_filename = query_fasta

        # Use diamond to extract protein sequence if input is nucleotide
        if query_type == "nucl":
            logger.debug("Input is nucleotide sequence, running Diamond first")
            protein_filename = run_diamond(
                db_list=db_list,
                query_type=query_type,
                input_gene=input_gene,
                input_eval=input_eval,
                query_fasta=query_fasta,
                threads=threads,
            )
            logger.debug(f"Diamond search complete. Output file: {protein_filename}")

        logger.debug(f"Running HMMER search on protein sequence: {protein_filename}")
        tmp_hmmer = hmmsearch(
            db_list=db_list,
            query_type="prot",  # Always use protein HMM
            input_gene=input_gene,
            input_eval=input_eval,
            query_fasta=protein_filename,
            threads=threads,
        )
        logger.debug(f"HMMER search complete. Output file: {tmp_hmmer}")

        logger.debug("Parsing HMMER results")
        tmp_hmmer_parsed = parse_hmmer(tmp_hmmer)
        logger.debug(f"Parsed HMMER results saved to: {tmp_hmmer_parsed}")

        hmmer_results = pd.read_csv(tmp_hmmer_parsed, sep="\t")
        logger.debug(
            f"HMMER results loaded into DataFrame with shape: {hmmer_results.shape}"
        )
        logger.debug(f"HMMER results columns: {hmmer_results.columns}")

        # Return protein_filename only for nucleotide input
        if query_type == "nucl":
            return hmmer_results.to_dict("records"), protein_filename
        else:
            return hmmer_results.to_dict("records"), None

    except Exception as e:
        logger.error(f"Error in HMMER search: {str(e)}")
        logger.exception("Full traceback:")
        return None, None


# Parse the HMMER results
def parse_hmmer(hmmer_output_file):
    """Parse HMMER output file."""
    try:
        logger.debug(f"Starting to parse HMMER output file: {hmmer_output_file}")

        parsed_file = tempfile.NamedTemporaryFile(suffix=".hmmer.parsed.txt").name
        logger.debug(f"Created temporary parsed output file: {parsed_file}")

        with open(parsed_file, "w") as tsv_file:
            tsv_file.write(
                "query_id\thit_IDs\taln_length\tquery_start\tquery_end\tgaps\tquery_seq\tsubject_seq\tevalue\tbitscore\n"
            )

            logger.debug("Parsing HMMER records")
            record_count = 0
            hit_count = 0

            for record in SearchIO.parse(hmmer_output_file, "hmmer3-text"):
                record_count += 1
                for hit in record.hits:
                    hit_count += 1
                    for hsp in hit.hsps:
                        query_seq = str(hsp.query.seq)
                        subject_seq = clean_shipID(str(hsp.hit.seq))
                        aln_length = hsp.aln_span
                        query_start = hsp.query_start
                        query_end = hsp.query_end
                        gaps = str("N/A")
                        bitscore = hsp.bitscore
                        evalue = hsp.evalue

                        tsv_file.write(
                            f"{hit.id}\t{record.id}\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t{query_seq}\t{subject_seq}\t{evalue}\t{bitscore}\n"
                        )
        return parsed_file

    except Exception as e:
        logger.error(f"Error parsing HMMER output: {str(e)}")
        logger.exception("Full traceback:")
        raise


# parse blast xml output to tsv
def parse_blast_xml(xml):
    parsed_file = tempfile.NamedTemporaryFile(suffix=".blast.parsed.txt").name
    with open(parsed_file, "w") as tsv_file:
        # Add quotes around field names to ensure proper parsing
        tsv_file.write(
            '"query_id"\t"hit_IDs"\t"aln_length"\t"query_start"\t"query_end"\t"gaps"\t"query_seq"\t"subject_seq"\t"evalue"\t"bitscore"\t"pident"\n'
        )
        record = SearchIO.read(xml, "blast-xml")
        for hit in record:
            for hsp in hit:
                try:
                    query_seq = str(hsp.query) if hsp.query else "N/A"
                    subject_seq = clean_shipID(str(hsp.hit)) if hsp.hit else "N/A"
                    aln_length = (
                        hsp.aln_span if hasattr(hsp, "aln_span") else len(query_seq)
                    )
                    query_start = hsp.query_start if hasattr(hsp, "query_start") else 0
                    query_end = (
                        hsp.query_end if hasattr(hsp, "query_end") else aln_length
                    )
                    gaps = str(hsp.gap_num) if hasattr(hsp, "gap_num") else "0"
                    bitscore = hsp.bitscore if hasattr(hsp, "bitscore") else 0.0
                    evalue = hsp.evalue if hasattr(hsp, "evalue") else 0.0

                    # Calculate percent identity
                    if hasattr(hsp, "ident_num") and hsp.aln_span > 0:
                        pident = (hsp.ident_num / hsp.aln_span) * 100
                    else:
                        identical_count = sum(
                            1 for q, h in zip(str(hsp.query), str(hsp.hit)) if q == h
                        )
                        pident = (
                            (identical_count / aln_length) * 100
                            if aln_length > 0
                            else 0
                        )

                    # Quote string values and format numbers
                    line = f'"{record.id}"\t"{hit.id}"\t{aln_length}\t{query_start}\t{query_end}\t{gaps}\t"{query_seq}"\t"{subject_seq}"\t{evalue}\t{bitscore}\t{pident}\n'
                    tsv_file.write(line)

                except Exception as e:
                    logger.error(f"Error processing HSP: {e}")
                    continue
    return parsed_file


def circos_prep(blast_output, links_output, layout_output):
    best_blast = blast_output.nsmallest(10, "evalue")
    links = [
        {
            "source": {
                "id": row["qseqid"],
                "start": int(row["qstart"]),
                "end": int(row["qend"]),
            },
            "target": {
                "id": row["sseqid"],
                "start": int(row["sstart"]),
                "end": int(row["send"]),
            },
            # Use pident as the single numeric value for simplicity
            "value": row["pident"],
        }
        for index, row in best_blast.iterrows()
    ]
    # Save to JSON file
    with open(links_output, "w") as json_file:
        json.dump(links, json_file, indent=4)

    # Determine unique sequence IDs and their lengths
    sequence_lengths = {}
    for index, row in best_blast.iterrows():
        if row["qseqid"] not in sequence_lengths:
            sequence_lengths[row["qseqid"]] = int(row["qend"])
        else:
            sequence_lengths[row["qseqid"]] = max(
                sequence_lengths[row["qseqid"]], int(row["qend"])
            )

        if row["sseqid"] not in sequence_lengths:
            sequence_lengths[row["sseqid"]] = int(row["send"])
        else:
            sequence_lengths[row["sseqid"]] = max(
                sequence_lengths[row["sseqid"]], int(row["send"])
            )

    # Define colors for each chromosome block
    colors = [
        "#FF5733",
        "#33FF57",
        "#3357FF",
        "#FF33A6",
        "#33FFF1",
        "#F1FF33",
        "#B633FF",
        "#33FFB6",
        "#FFA633",
        "#33A6FF",
    ]

    # Create Circos layout with colors
    layout = [
        {
            "id": seq_id,
            "len": length,
            "color": colors[i % len(colors)],
            "label": {"display": seq_id, "color": "#000000", "size": 14},
        }
        for i, (seq_id, length) in enumerate(sequence_lengths.items())
    ]

    # Save layout to JSON file
    with open(layout_output, "w") as json_file:
        json.dump(layout, json_file, indent=4)


def blast_download_button():
    """Creates the download button for BLAST results."""
    return html.Div(
        [
            dmc.Space(h="xl"),
            dmc.Center(
                dmc.Button(
                    "Download BLAST Results",
                    id="blast-dl-button",
                    variant="gradient",
                    gradient={"from": "indigo", "to": "cyan"},
                    size="lg",
                    leftSection=[html.I(className="bi bi-download")],
                )
            ),
            dcc.Download(id="blast-dl"),
        ]
    )


def select_ship_family(hmmer_results):
    """Process HMMER results and return family information."""
    if isinstance(hmmer_results, dict):
        # Handle dictionary input
        return (
            hmmer_results.get("hit_IDs"),
            hmmer_results.get("aln_length"),
            hmmer_results.get("evalue"),
        )

    # Handle DataFrame input
    try:
        hmmer_results["evalue"] = pd.to_numeric(
            hmmer_results["evalue"], errors="coerce"
        )
        hmmer_results.dropna(subset=["evalue"], inplace=True)

        if len(hmmer_results) == 0:
            logger.warning("HMMER results DataFrame is empty after dropping NaNs.")
            return None, None, None

        hmmer_results_sorted = hmmer_results.sort_values(by=["query_id", "evalue"])
        best_matches = hmmer_results_sorted.drop_duplicates(
            subset="query_id", keep="first"
        )

        if len(best_matches) > 0:
            best_match = best_matches.iloc[0]
            superfamily = best_match["hit_IDs"]
            aln_length = best_match["aln_length"]
            try:
                evalue_num = float(best_match["evalue"])
                evalue = "{:.2e}".format(evalue_num)
            except (ValueError, TypeError, KeyError) as e:
                logger.error(f"Error formatting e-value: {e}")
                evalue = str(best_match.get("evalue", "N/A"))

            return superfamily, aln_length, evalue

        logger.warning("No valid rows found in hmmer_results DataFrame.")
        return None, None, None

    except Exception as e:
        logger.error(f"Error processing HMMER results: {e}")
        return None, None, None


# TODO: add qseq_translated to the output
def run_diamond(
    db_list=None,
    query_type=None,
    input_gene=None,
    input_eval=0.01,
    query_fasta=None,
    threads=2,
):
    """Run DIAMOND search against protein database.

    Args:
        db_list: Dictionary containing database paths
        query_type: Type of query sequence ('nucl' or 'prot')
        input_gene: Gene type to search against (default: 'tyr')
        input_eval: E-value threshold
        query_fasta: Path to query FASTA file
        threads: Number of threads to use

    Returns:
        List of dictionaries containing DIAMOND results
    """
    try:
        diamond_out = tempfile.NamedTemporaryFile(suffix=".tsv", delete=False).name

        try:
            diamond_db = db_list["gene"][input_gene]["prot"]
        except KeyError:
            raise ValueError(f"Database path not found for gene type: {input_gene}")

        if not os.path.exists(diamond_db) or os.path.getsize(diamond_db) == 0:
            raise ValueError(f"DIAMOND database {diamond_db} not found or is empty")

        blast_type = "blastx" if query_type == "nucl" else "blastp"
        evalue = input_eval if input_eval else 0.001

        diamond_cmd = [
            "diamond",
            blast_type,
            "--db",
            diamond_db,
            "-q",
            query_fasta,
            "--outfmt",
            "6",
            "-e",
            str(evalue),
            "--strand",
            "both",
            "-p",
            str(threads),
            "-k",
            "1",
            "--skip-missing-seqids",
            "-o",
            diamond_out,
        ]

        try:
            subprocess.run(diamond_cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"DIAMOND error: {e.stderr}")
            raise

        if os.path.getsize(diamond_out) > 0:
            column_names = [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ]
            diamond_results = pd.read_csv(diamond_out, sep="\t", names=column_names)
            return diamond_results.to_dict("records")
        else:
            logger.warning("No DIAMOND hits found")
            return []

    except Exception as e:
        logger.error(f"Error running DIAMOND: {str(e)}")
        raise

    finally:
        if "diamond_out" in locals() and os.path.exists(diamond_out):
            os.unlink(diamond_out)


def make_captain_alert(family, aln_length, evalue, search_type):
    try:
        # Validate inputs
        if family is None or (isinstance(family, str) and not family.strip()):
            logger.error(f"Invalid family parameter: {family}")
            return create_error_alert("Invalid family name")

        # Convert pandas Series to scalar if needed
        if hasattr(family, "iloc"):
            family = family.iloc[0]
        if hasattr(aln_length, "iloc"):
            aln_length = aln_length.iloc[0]
        if hasattr(evalue, "iloc"):
            evalue = evalue.iloc[0]

        if not aln_length:
            logger.error(f"Invalid alignment length: {aln_length}")
            return create_error_alert("Invalid alignment length")

        try:
            evalue_num = float(evalue)
            if evalue_num == 0.0:
                formatted_evalue = "0.00e+00"
            else:
                formatted_evalue = "{:.2e}".format(evalue_num)
        except (ValueError, TypeError) as e:
            logger.error(f"Error formatting e-value {evalue}: {e}")
            return create_error_alert("Invalid e-value format")

        if search_type not in ["blast", "hmmsearch"]:
            logger.error(f"Invalid search type: {search_type}")
            return create_error_alert("Invalid search type")

        return dmc.Alert(
            title=f"Family Classification via {search_type.upper()} Search",
            children=[
                f"Your sequence is likely in Starship family: {family}",
                dmc.Space(h=5),
                dmc.Text(
                    f"BLAST Search: Alignment length = {aln_length}, E-value = {formatted_evalue}",
                    size="sm",
                    c="dimmed",
                ),
            ],
            color="blue",
            variant="light",
            withCloseButton=False,
            style={
                "width": "100%",
                "margin": "0 auto",
                "@media (min-width: 768px)": {"width": "50%"},
            },
        )

    except Exception as e:
        logger.error(f"Error in make_captain_alert: {e}")
        return create_error_alert(str(e))


def process_captain_results(captain_results_dict=None, evalue=None):
    """Process captain results and return appropriate alert component."""
    if captain_results_dict is None or not captain_results_dict:
        logger.warning("No captain results provided")
        return html.Div(create_no_matches_alert())

    try:
        # Handle case where captain_results_dict is a list
        if isinstance(captain_results_dict, list):
            if not captain_results_dict:  # Empty list
                return html.Div(create_no_matches_alert())
            result = captain_results_dict[0]  # Use first result
        else:
            # Handle single dictionary case
            result = captain_results_dict

        # Extract relevant information directly from the result
        superfamily = result.get("hit_IDs")
        family_aln_length = result.get("aln_length")
        family_evalue = result.get("evalue", evalue)  # Use provided evalue as fallback

        if superfamily:
            return make_captain_alert(
                superfamily, family_aln_length, family_evalue, search_type="hmmsearch"
            )

        return html.Div(create_no_matches_alert())

    except Exception as e:
        logger.error(f"Error processing captain results: {str(e)}")
        return html.Div(create_error_alert("Failed to process captain results"))


def create_no_matches_alert():
    return dmc.Alert(
        title="No Matches Found",
        children=[
            dmc.Text("Your sequence did not match any Starships in our database."),
            dmc.Space(h=10),
            dmc.Text("Suggestions:", size="sm", fw=500),
            dmc.Space(h=5),
            dmc.List(
                [
                    dmc.ListItem(
                        dmc.Text(
                            "Try searching the full database (including un-curated sequences)",
                            size="sm",
                        )
                    ),
                    dmc.ListItem(
                        dmc.Text(
                            "Try searching with a different region of your sequence",
                            size="sm",
                        )
                    ),
                    dmc.ListItem(
                        dmc.Text(
                            "Consider using a less stringent E-value threshold",
                            size="sm",
                        )
                    ),
                ],
                withPadding=True,
                spacing="xs",
                size="sm",
                type="ordered",
            ),
        ],
        color="yellow",
        variant="light",
        withCloseButton=False,
        style={
            "width": "100%",
            "margin": "0 auto",
            "@media (min-width: 768px)": {"width": "50%"},
        },
    )


def get_blast_db(query_type):
    """Get the appropriate BLAST database configuration based on query type."""
    from src.config.settings import BLAST_DB_PATHS

    if query_type == "nucl":
        return BLAST_DB_PATHS["ship"]["nucl"]
    elif query_type == "prot":
        return BLAST_DB_PATHS["gene"]["tyr"]["prot"]
    else:
        raise ValueError(f"Invalid query type: {query_type}")

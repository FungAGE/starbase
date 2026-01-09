import subprocess
import os
import glob
import fcntl
import pandas as pd
from Bio import SeqIO
from sourmash import MinHash, SourmashSignature, save_signatures, load_file_as_signatures
import screed

from src.config.cache import cache
from src.config.settings import BLAST_DB_PATHS
from src.utils.seq_utils import create_ncbi_style_header, write_fasta

from src.config.logging import get_logger

logger = get_logger(__name__)


class FileLock:
    """Simple filesystem lock using fcntl (Linux/Unix)."""

    def __init__(self, lock_path: str):
        self.lock_path = lock_path
        self._fh = None

    def __enter__(self):
        os.makedirs(os.path.dirname(self.lock_path), exist_ok=True)
        self._fh = open(self.lock_path, "w")
        fcntl.flock(self._fh.fileno(), fcntl.LOCK_EX)
        return self

    def __exit__(self, exc_type, exc, tb):
        try:
            if self._fh is not None:
                fcntl.flock(self._fh.fileno(), fcntl.LOCK_UN)
                self._fh.close()
        finally:
            self._fh = None

def create_blast_database(fasta_path, dbtype):
    """
    Create a BLAST database from the provided FASTA file.

    - Uses the FASTA basename (without extension) as the DB prefix
    - Proactively removes any stale BLAST index/taxonomy files to avoid taxonomy LMDB conflicts
    """
    # Use basename without extension as the BLAST DB prefix
    out_prefix, _ = os.path.splitext(fasta_path)                                                                                                                                                                                

    # Best-effort cleanup of old index/taxonomy files (prevents taxonomy-related crashes)
    stale_suffixes = [
        ".nhr",
        ".nin",
        ".nsq",
        ".ndb",
        ".not",
        ".ntf",
        ".nto",
        ".njs",
    ]
    for suffix in stale_suffixes:
        stale_path = f"{out_prefix}{suffix}"
        try:
            if os.path.exists(stale_path):
                os.remove(stale_path)
        except Exception:
            # Ignore cleanup failures
            pass

    cmd = [
        "makeblastdb",
        "-in",
        fasta_path,
        "-input_type",
        "fasta",
        "-dbtype",
        dbtype,
        "-out",
        out_prefix,
    ]

    logger.debug(f"Creating BLAST database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("BLAST database created successfully.")
    except subprocess.CalledProcessError as e:
        stderr = e.stderr or ""
        stdout = e.stdout or ""
        logger.error(
            "Failed to create BLAST database (exit %s).\nSTDOUT:\n%s\nSTDERR:\n%s",
            e.returncode,
            stdout.strip(),
            stderr.strip(),
        )
        raise


def blast_db_exists(blastdb):
    """
    Check if the BLAST database exists by looking for required file extensions.
    """
    # Determine the expected BLAST DB prefix (basename without extension)
    db_prefix, _ = os.path.splitext(blastdb)

    # Required core files for BLAST DB existence check
    required_suffixes = [".nhr", ".nin", ".nsq"]

    for suffix in required_suffixes:
        expected_path = f"{db_prefix}{suffix}"
        if not os.path.exists(expected_path):
            return False

    return True


def create_sourmash_signatures(fasta_path, seq_type="nucl"):
    """
    Create sourmash signatures from a FASTA file and save to disk.
    
    Args:
        fasta_path: Path to input FASTA file
        seq_type: 'nucl' for nucleotide or 'prot' for protein
        
    Returns:
        Path to the created signature file (.sig)
    """
    sig_path = fasta_path + ".sig"
    
    # Set parameters based on sequence type
    if seq_type == "nucl":
        k = 21
        scaled = 1000
        is_protein = False
    else:
        k = 7
        scaled = 100
        is_protein = True
    
    logger.debug(f"Creating sourmash signatures for {fasta_path} (type={seq_type})")
    
    try:
        signatures = []
        seq_count = 0
        
        # Read FASTA and create signatures
        for record in screed.open(fasta_path):
            mh = MinHash(n=0, ksize=k, scaled=scaled, is_protein=is_protein)
            
            # Add sequence data
            if is_protein:
                mh.add_protein(record.sequence)
            else:
                mh.add_sequence(record.sequence, force=True)
            
            # Create signature with sequence name
            sig = SourmashSignature(mh, name=record.name)
            signatures.append(sig)
            seq_count += 1
        
        # Save all signatures to file
        with open(sig_path, 'w') as f:
            save_signatures(signatures, f)
        
        logger.info(f"Created {seq_count} sourmash signatures: {sig_path}")
        return sig_path
        
    except Exception as e:
        logger.error(f"Failed to create sourmash signatures: {e}")
        raise


def sourmash_sig_exists(fasta_path):
    """
    Check if sourmash signature file exists for a given FASTA file.
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        bool: True if .sig file exists
    """
    sig_path = fasta_path + ".sig"
    return os.path.exists(sig_path)


def load_sourmash_signatures(fasta_path):
    """
    Load pre-computed sourmash signatures from disk.
    
    Args:
        fasta_path: Path to FASTA file (signature file is fasta_path + '.sig')
        
    Returns:
        List of (seq_name, signature) tuples or None if file doesn't exist
    """
    sig_path = fasta_path + ".sig"
    
    if not os.path.exists(sig_path):
        logger.warning(f"Signature file not found: {sig_path}")
        return None
    
    try:
        signatures = []
        for sig in load_file_as_signatures(sig_path):
            signatures.append((sig.name, sig))
        
        logger.debug(f"Loaded {len(signatures)} signatures from {sig_path}")
        return signatures
        
    except Exception as e:
        logger.error(f"Failed to load signatures from {sig_path}: {e}")
        return None


def create_diamond_database(fasta_path, threads=2):
    # Create output path with .dmnd extension
    diamond_db = fasta_path + ".dmnd"

    cmd = [
        "diamond",
        "makedb",
        "--in",
        fasta_path,
        "--db",
        diamond_db,
        "--threads",
        str(threads),
    ]

    logger.debug(f"Creating Diamond database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug("Diamond database created successfully.")
        return diamond_db  # Return the path to the created Diamond database
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create Diamond database: {e.stderr}")
        raise


def create_dbs():
    from src.database.sql_manager import (
        fetch_captains,
        fetch_ships,
        fetch_meta_data,
    )

    datasets = [
        ("all", False),
        ("curated", True),
    ]

    # Serialize rebuilds across processes using a file lock
    lock_file = os.path.join(
        os.path.dirname(BLAST_DB_PATHS["ship"]["all"]["nucl"]),
        ".rebuild.lock",
    )

    with FileLock(lock_file):
        for dataset_name, curated_filter in datasets:
            logger.info(f"Creating {dataset_name} ship BLAST database...")

            ship_fasta_path = BLAST_DB_PATHS["ship"][dataset_name]["nucl"]
            ship_fasta_dir = os.path.dirname(ship_fasta_path)
            os.makedirs(ship_fasta_dir, exist_ok=True)

            ship_sequences = fetch_ships(
                dereplicate=True, with_sequence=True, curated=curated_filter
            )

            if ship_sequences.empty:
                logger.warning(f"No sequences found for {dataset_name} dataset")
                continue

            ship_metadata = fetch_meta_data(
                accession_tags=ship_sequences["accession_tag"].tolist()
            )

            ship_sequences_dict = {}

            for _, row in ship_sequences.iterrows():
                # Convert Series to scalar values to avoid unhashable type errors
                accession_tag = str(row["accession_tag"]) if pd.notna(row["accession_tag"]) else None
                sequence = str(row["sequence"]) if pd.notna(row["sequence"]) else None

                if not accession_tag or not sequence:
                    continue

                metadata_rows = ship_metadata[
                    ship_metadata["accession_tag"] == accession_tag
                ]
                if not metadata_rows.empty:
                    # Ensure all values are scalars, not pandas objects
                    metadata = {}
                    for col in metadata_rows.columns:
                        val = metadata_rows.iloc[0][col]
                        # Convert to Python types
                        if pd.isna(val):
                            metadata[col] = None
                        else:
                            metadata[col] = val.item() if hasattr(val, 'item') else val

                    header = create_ncbi_style_header(metadata)
                    if header and sequence:
                        # Ensure header is a string to avoid unhashable type errors
                        header = str(header) if header else None
                        if header:
                            ship_sequences_dict[header] = sequence

            # rebuild BLAST database everytime?
            logger.info(f"Rebuilding {dataset_name} BLAST database...")
            write_fasta(ship_sequences_dict, ship_fasta_path)
            create_blast_database(ship_fasta_path, "nucl")
            logger.info(f"{dataset_name} BLAST database rebuilt successfully")
            
            # Create sourmash signatures
            logger.info(f"Creating sourmash signatures for {dataset_name} dataset...")
            try:
                create_sourmash_signatures(ship_fasta_path, seq_type="nucl")
                logger.info(f"{dataset_name} sourmash signatures created successfully")
            except Exception as e:
                logger.error(f"Failed to create sourmash signatures for {dataset_name}: {e}")

    # Create captain database
    captain_fasta_path = BLAST_DB_PATHS["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(captain_fasta_dir, exist_ok=True)

    captain_sequences = cache.get("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_captains(dereplicate=True, with_sequence=True)

    captain_sequences_dict = {}
    for index, row in captain_sequences.iterrows():
        captain_id_val = row.get("captain_id_col", row.get("captainID"))
        sequence_val = row["sequence"]

        if pd.notna(captain_id_val) and pd.notna(sequence_val):
            accession = str(captain_id_val)
            sequence = str(sequence_val)
            captain_sequences_dict[accession] = sequence

    # Captain DBs: skip if already exist
    if blast_db_exists(captain_fasta_path):
        logger.info("Captain BLAST DB already exists; skipping rebuild")
    else:
        write_fasta(captain_sequences_dict, captain_fasta_path)
        create_blast_database(captain_fasta_path, "prot")
        create_diamond_database(captain_fasta_path)
        
        # Create sourmash signatures for captain sequences
        logger.info("Creating sourmash signatures for captain sequences...")
        try:
            create_sourmash_signatures(captain_fasta_path, seq_type="prot")
            logger.info("Captain sourmash signatures created successfully")
        except Exception as e:
            logger.error(f"Failed to create captain sourmash signatures: {e}")
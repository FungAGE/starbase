import subprocess
import os
import glob
import fcntl

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
    # Use basename without extension as the BLAST DB prefix (avoid using the .fa suffix)
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
                accession_tag=ship_sequences["accession_tag"].tolist()
            )

            ship_sequences_dict = {}

            for _, row in ship_sequences.iterrows():
                accession_tag = row["accession_tag"]
                sequence = row["sequence"]

                metadata_rows = ship_metadata[
                    ship_metadata["accession_tag"] == accession_tag
                ]
                if not metadata_rows.empty:
                    metadata = metadata_rows.iloc[0]

                    header = create_ncbi_style_header(metadata)
                    if header and sequence:
                        ship_sequences_dict[header] = sequence

            # If DB already exists, skip rebuild
            if blast_db_exists(ship_fasta_path):
                logger.info(
                    f"BLAST DB already exists for {dataset_name}; skipping rebuild"
                )
            else:
                logger.info(
                    f"Writing {len(ship_sequences_dict)} sequences to {dataset_name} FASTA file"
                )
                write_fasta(ship_sequences_dict, ship_fasta_path)
                create_blast_database(ship_fasta_path, "nucl")

    # Create captain database
    captain_fasta_path = BLAST_DB_PATHS["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(captain_fasta_dir, exist_ok=True)

    captain_sequences = cache.get("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_captains(dereplicate=True, with_sequence=True)

    captain_sequences_dict = {}
    for index, row in captain_sequences.iterrows():
        accession = row["captainID"]
        sequence = row["sequence"]
        captain_sequences_dict[accession] = sequence

    # Captain DBs: skip if already exist
    if blast_db_exists(captain_fasta_path):
        logger.info("Captain BLAST DB already exists; skipping rebuild")
    else:
        write_fasta(captain_sequences_dict, captain_fasta_path)
        create_blast_database(captain_fasta_path, "prot")
        create_diamond_database(captain_fasta_path)

import subprocess
import os
import glob
import logging

from src.config.cache import cache
from src.config.settings import BLAST_DB_PATHS
from src.utils.seq_utils import write_fasta
logger = logging.getLogger(__name__)


def create_blast_database(fasta_path, dbtype):
    cmd = [
        "makeblastdb",
        "-in",
        fasta_path,
        "-input_type",
        "fasta",
        "-dbtype",
        dbtype,
        "-out",
        fasta_path,
    ]

    logger.info(f"Creating BLAST database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logger.info("BLAST database created successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create BLAST database: {e}")
        raise


def blast_db_exists(blastdb):
    """
    Check if the BLAST database exists by looking for required file extensions.
    """
    db_directory = os.path.dirname(blastdb)
    extensions = ["*.ndb", "*.nhr", "*.nin", "*.not", "*.nsq", "*.ntf", "*.nto"]

    for ext in extensions:
        if glob.glob(os.path.join(db_directory, ext)):
            return True

    return False

def create_diamond_database(fasta_path, threads=2):
    cmd = [
        "/usr/bin/diamond prepdb",
        "--db",
        fasta_path,
        "--threads",
        str(threads),
    ]

    logger.info(f"Creating Diamond database with command: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        logger.info("Diamond database created successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to create Diamond database: {e}")
        raise


def create_dbs():
    from src.database.sql_manager import fetch_captains, fetch_ships

    # Create ship database
    ship_fasta_path = BLAST_DB_PATHS["ship"]
    ship_fasta_dir = os.path.dirname(ship_fasta_path)
    os.makedirs(ship_fasta_dir, exist_ok=True)

    ship_sequences = cache.get("all_ships")
    if ship_sequences is None:
        ship_sequences = fetch_ships(with_sequence=True)

    ship_sequences_dict = {
        row["accession_tag"]: row["sequence"] 
        for _, row in ship_sequences.iterrows()
    }

    write_fasta(ship_sequences_dict, ship_fasta_path)
    create_blast_database(ship_fasta_path, "nucl")

    # Create captain database
    captain_fasta_path = BLAST_DB_PATHS["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(captain_fasta_dir, exist_ok=True)

    captain_sequences = cache.get("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_captains(with_sequence=True)

    captain_sequences_dict = {
        row["captainID"]: row["sequence"]
        for _, row in captain_sequences.iterrows()
    }

    write_fasta(captain_sequences_dict, captain_fasta_path)
    create_blast_database(captain_fasta_path, "prot")
    create_diamond_database(captain_fasta_path, "prot")

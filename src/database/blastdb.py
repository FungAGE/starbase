import subprocess
import os
import glob
import logging

from src.config.cache import cache
from src.config.database import DB_PATHS

logger = logging.getLogger(__name__)

DB_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "database", "db")

os.makedirs(DB_DIR, exist_ok=True)
db_list = {
    "ship": {"nucl": f"{DB_DIR}/ships/fna/blastdb/ships.fa"},
    "gene": {
        "tyr": {
            "nucl": f"{DB_DIR}/captain/tyr/fna/blastdb/captains.fna",
            "prot": f"{DB_DIR}/captain/tyr/faa/blastdb/captains.faa",
            "hmm": {
                "nucl": f"{DB_DIR}/captain/tyr/fna/hmm/combined.hmm",
                "prot": f"{DB_DIR}/captain/tyr/faa/hmm/combined.hmm",
            },
        },
    },
}


def write_fasta(sequences, fasta_path):
    with open(fasta_path, "w") as fasta_file:
        for name, sequence in sequences:
            fasta_file.write(f">{name}\n{sequence}\n")


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


def create_dbs():
    
    from src.database.sql_manager import fetch_all_captains, fetch_all_ships

    ship_fasta_path = db_list["ship"]
    ship_fasta_dir = os.path.dirname(ship_fasta_path)
    os.makedirs(ship_fasta_dir, exist_ok=True)

    ship_sequences_list = []
    ship_sequences = cache.get("all_ships")
    if ship_sequences is None:
        ship_sequences = fetch_all_ships()

    for index, row in ship_sequences.iterrows():
        name = row["accession_tag"]
        sequence = row["sequence"]
        ship_sequences_list.append((name, sequence))

    write_fasta(ship_sequences_list, ship_fasta_path)
    create_blast_database(ship_fasta_path, "nucl")

    captain_fasta_path = db_list["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(
        captain_fasta_dir, exist_ok=True
    )

    captain_sequences_list = []
    captain_sequences = cache.get("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_all_captains()

    for index, row in captain_sequences.iterrows():
        name = row["captainID"]
        sequence = row["sequence"]
        captain_sequences_list.append((name, sequence))

    write_fasta(captain_sequences_list, captain_fasta_path)
    create_blast_database(captain_fasta_path, "prot")

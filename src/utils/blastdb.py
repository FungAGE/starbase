import subprocess
import os
import glob
import logging

logger = logging.getLogger(__name__)


db_list = {
    "ship": {"nucl": "src/data/ships/fna/blastdb/ships.fa"},
    "gene": {
        "tyr": {
            "nucl": "src/data/captain/tyr/fna/blastdb/captains.fna",
            "prot": "src/data/captain/tyr/faa/blastdb/captains.faa",
            "hmm": {
                "nucl": "src/data/captain/tyr/fna/hmm/combined.hmm",
                "prot": "src/data/captain/tyr/faa/hmm/combined.hmm",
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

    # Check for the existence of files with the specified extensions
    for ext in extensions:
        if glob.glob(os.path.join(db_directory, ext)):
            return True

    return False


def create_dbs():
    from src.components.cache_manager import load_from_cache
    from src.components.sql_queries import fetch_all_captains, fetch_all_ships

    # Create BLAST database for ships
    ship_fasta_path = db_list["ship"]["nucl"]
    ship_fasta_dir = os.path.dirname(ship_fasta_path)
    os.makedirs(ship_fasta_dir, exist_ok=True)  # Create directory if it doesn't exist

    ship_sequences_list = []
    ship_sequences = load_from_cache("all_ships")
    if ship_sequences is None:
        ship_sequences = fetch_all_ships()

    for index, row in ship_sequences.iterrows():
        name = row["accession_tag"]
        sequence = row["sequence"]
        ship_sequences_list.append((name, sequence))

    # Fix the parameter passed to write_fasta
    write_fasta(ship_sequences_list, ship_fasta_path)
    create_blast_database(ship_fasta_path, "nucl")

    # Create BLAST database for captains
    captain_fasta_path = db_list["gene"]["tyr"]["prot"]
    captain_fasta_dir = os.path.dirname(captain_fasta_path)
    os.makedirs(
        captain_fasta_dir, exist_ok=True
    )  # Create directory if it doesn't exist

    captain_sequences_list = []
    captain_sequences = load_from_cache("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_all_captains()

    for index, row in captain_sequences.iterrows():
        name = row["captainID"]
        sequence = row["sequence"]
        captain_sequences_list.append((name, sequence))

    write_fasta(captain_sequences_list, captain_fasta_path)
    create_blast_database(captain_fasta_path, "prot")

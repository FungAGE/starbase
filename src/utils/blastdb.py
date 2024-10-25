import subprocess
import os
from src.components.cache_manager import load_from_cache
from src.components.sql_queries import fetch_all_captains, fetch_all_ships

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
        # "nlr": {
        #     "nucl": "src/data/cargo/nlr/fna/blastdb/nlr.fa",
        #     "prot": "src/data/cargo/nlr/faa/blastdb/nlr.mycoDB.faa",
        # },
        # "fre": {
        #     "nucl": "src/data/cargo/fre/fna/blastdb/fre.fa",
        #     "prot": "src/data/cargo/fre/faa/blastdb/fre.mycoDB.faa",
        # },
        # "plp": {
        #     "nucl": "src/data/cargo/plp/fna/blastdb/plp.fa",
        #     "prot": "src/data/cargo/plp/faa/blastdb/plp.mycoDB.faa",
        # },
        # "duf3723": {
        #     "nucl": "src/data/cargo/duf3723/fna/blastdb/duf3723.fa",
        #     "prot": "src/data/cargo/duf3723/faa/blastdb/duf3723.mycoDB.faa",
        # },
    },
}


def write_fasta(sequences, fasta_path):
    with open(fasta_path, "w") as fasta_file:
        for name, sequence in sequences:
            fasta_file.write(f">{name}\n{sequence}\n")


def create_blast_database(fasta_path, dbtype):
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            fasta_path,
            "-input_type",
            "fasta",
            "-dbtype",
            dbtype,
            # "-parse_seqids",
            "-out",
            fasta_path,
        ],
        check=True,
    )


def create_dbs():
    ship_fasta_path = db_list["ship"]["nucl"]
    if not os.path.exists(ship_fasta_path):
        os.makedirs(ship_fasta_path)

    ship_sequences_list = []
    ship_sequences = load_from_cache("all_ships")
    if ship_sequences is None:
        ship_sequences = fetch_all_ships()
    for index, row in ship_sequences.iterrows():
        name = row["accession_tag"]
        sequence = row["sequence"]
        ship_sequences_list.append((name, sequence))

    write_fasta(ship_sequences, ship_fasta_path)
    create_blast_database(ship_fasta_path, "nucl")

    captain_fasta_path = db_list["gene"]["tyr"]["prot"]
    if not os.path.exists(captain_fasta_path):
        os.makedirs(captain_fasta_path)

    captain_sequences_list = []
    captain_sequences = load_from_cache("all_captains")
    if captain_sequences is None:
        captain_sequences = fetch_all_captains()
    for index, row in captain_sequences.iterrows():
        name = row["captainID"]
        sequence = row["sequence"]
        captain_sequences_list.append((name, sequence))

    write_fasta(captain_sequences, captain_fasta_path)
    create_blast_database(captain_fasta_path, "prot")

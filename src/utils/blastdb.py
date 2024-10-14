from sqlalchemy import create_engine, MetaData, Table
import subprocess

db_list = {
    "ship": {"nucl": "src/data/db/ships/fna/blastdb/ships.fa"},
    "gene": {
        "tyr": {
            "nucl": "src/data/db/captain/tyr/fna/blastdb/captains.fna",
            "prot": "src/data/db/captain/tyr/faa/blastdb/captains.faa",
            "hmm": {
                "nucl": "src/data/db/captain/tyr/fna/hmm/combined.hmm",
                "prot": "src/data/db/captain/tyr/faa/hmm/combined.hmm",
            },
        },
        "nlr": {
            "nucl": "src/data/db/cargo/nlr/fna/blastdb/nlr.fa",
            "prot": "src/data/db/cargo/nlr/faa/blastdb/nlr.mycoDB.faa",
        },
        "fre": {
            "nucl": "src/data/db/cargo/fre/fna/blastdb/fre.fa",
            "prot": "src/data/db/cargo/fre/faa/blastdb/fre.mycoDB.faa",
        },
        "plp": {
            "nucl": "src/data/db/cargo/plp/fna/blastdb/plp.fa",
            "prot": "src/data/db/cargo/plp/faa/blastdb/plp.mycoDB.faa",
        },
        "duf3723": {
            "nucl": "src/data/db/cargo/duf3723/fna/blastdb/duf3723.fa",
            "prot": "src/data/db/cargo/duf3723/faa/blastdb/duf3723.mycoDB.faa",
        },
    },
}


def fetch_sequences(db_url, type):
    engine = create_engine(db_url)
    connection = engine.connect()
    metadata = MetaData()
    metadata.reflect(bind=engine)

    sequences_table = Table(type, metadata, autoload_with=engine)
    result = connection.execute(sequences_table.select())

    if type == "captains":
        rowID = "captainID"
        sequenceID = "sequence"
    elif type == "ships":
        rowID = "accession_tag"
        sequenceID = "ship_sequence"
    else:
        raise ValueError("Unsupported table type")

    sequences = []
    for row in result.mappings():  # Use .mappings() to return rows as dictionaries
        name = row[rowID]
        sequence = row[sequenceID]
        sequences.append((name, sequence))

    connection.close()
    return sequences


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


if __name__ == "__main__":
    db_url = "sqlite:///src/data/db/starbase.sqlite"
    ship_fasta_path = db_list["ship"]["nucl"]
    ship_sequences = fetch_sequences(db_url, "ships")
    write_fasta(ship_sequences, ship_fasta_path)
    create_blast_database(ship_fasta_path, "nucl")

    captain_fasta_path = db_list["gene"]["tyr"]["prot"]
    captain_sequences = fetch_sequences(db_url, "captains")
    write_fasta(captain_sequences, captain_fasta_path)
    create_blast_database(captain_fasta_path, "prot")

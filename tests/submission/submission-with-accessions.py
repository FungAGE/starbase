from src.utils.seq_utils import load_fasta_to_dict
from src.utils.classification_utils import assign_accession

from src.config.logging import get_logger
from src.database.sql_manager import fetch_ships

logger = get_logger(__name__)

fasta_file = "/home/adrian/Downloads/Hephaestus_all.fasta"
fasta_dict = load_fasta_to_dict(fasta_file)


existing_ships = fetch_ships(curated=False, dereplicate=True)

for header, sequence in fasta_dict.items():
    accession, needs_review = assign_accession(sequence, existing_ships)

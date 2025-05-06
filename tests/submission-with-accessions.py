from src.utils.seq_utils import load_fasta_to_dict
from src.utils.classification_utils import assign_accession
from src.config.database import StarbaseSession
import pandas as pd
from src.config.logging import get_logger

logger = get_logger(__name__)

fasta_file = "/home/adrian/Downloads/Hephaestus_all.fasta"
fasta_dict = load_fasta_to_dict(fasta_file)


def fetch_ships(accession_tags=None, curated=False, dereplicate=True):
    """
    Fetch ship data for specified accession tags.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all ships.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.

    Returns:
        pd.DataFrame: DataFrame containing ship data
    """
    session = StarbaseSession()

    query = """
    WITH valid_ships AS (
        SELECT DISTINCT 
            a.id as accession_id, 
            a.accession_tag,
            j.curated_status,
            j.elementBegin,
            j.elementEnd,
            j.contigID,
            t.species,
            t.genus,
            t.family,
            t.`order`,
            f.familyName,
            g.assembly_accession
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN taxonomy t ON j.taxid = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        WHERE 1=1
    """

    if accession_tags:
        query += " AND a.accession_tag IN ({})".format(
            ",".join(f"'{tag}'" for tag in accession_tags)
        )
    if curated:
        query += " AND j.curated_status = 'curated'"

    query += """
    )
    SELECT 
        v.*,
        s.sequence
    FROM valid_ships v
    INNER JOIN ships s ON s.accession = v.accession_id
    """

    try:
        df = pd.read_sql_query(query, session.bind)

        if dereplicate:
            df = df.drop_duplicates(subset="accession_tag")

        if df.empty:
            logger.warning("Fetched ships DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching ships data: {str(e)}")
        raise
    finally:
        session.close()


existing_ships = fetch_ships(curated=False, dereplicate=True)

for header, sequence in fasta_dict.items():
    accession, needs_review = assign_accession(sequence, existing_ships)

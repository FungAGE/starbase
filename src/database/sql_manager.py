import pandas as pd
from src.config.database import StarbaseSession
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)
from contextlib import contextmanager
import sqlalchemy.exc
from sqlalchemy import text
from src.config.cache import cache
from src.config.settings import PHYLOGENY_PATHS

from src.config.logging import get_logger

logger = get_logger(__name__)


# Define a common retry decorator for database operations
def db_retry_decorator(additional_retry_exceptions=()):
    """
    Create a retry decorator for database operations
    Args:
        additional_retry_exceptions: Tuple of additional exceptions to retry on
    """
    retry_exceptions = (sqlalchemy.exc.OperationalError,) + additional_retry_exceptions

    return retry(
        stop=stop_after_attempt(3),
        wait=wait_exponential(multiplier=1, min=4, max=10),
        retry=retry_if_exception_type(retry_exceptions),
        before_sleep=lambda retry_state: logger.warning(
            f"Retrying database operation after error: {retry_state.outcome.exception()}"
        ),
    )


# Context manager for database sessions with timeout
@contextmanager
def db_session_manager():
    """Context manager for database sessions with timeout"""
    session = None
    try:
        session = StarbaseSession()
        # SQLite doesn't support SET SESSION, so we'll skip the timeout setting
        if session.bind.dialect.name != "sqlite":
            session.execute(text("SET SESSION wait_timeout=30"))  # Only for MySQL
        yield session
    except Exception as e:
        logger.error(f"Database error: {str(e)}")
        if session:
            session.rollback()
        raise
    finally:
        if session:
            session.close()


@db_retry_decorator()
def fetch_meta_data(curated=False, accession_tags=None):
    """
    Fetch metadata from the database with caching.

    Args:
        curated (bool): If True, only return curated entries
        accession_tags (str or list): Single accession tag or list of accession tags

    Returns:
        pd.DataFrame: Metadata for the specified accession tags
    """
    session = StarbaseSession()

    meta_query = """
    SELECT j.curated_status, j.starshipID, j.ship_id,
           a.accession_tag, a.version_tag,
           CASE 
               WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
               THEN a.accession_tag || '.' || a.version_tag
               ELSE a.accession_tag
           END as accession_display,
           t.taxID, t.strain, t.`order`, t.family, t.name, 
           sf.elementLength, sf.upDR, sf.downDR, sf.contigID, sf.captainID, sf.elementBegin, sf.elementEnd, 
           f.familyName, f.type_element_reference, n.navis_name, h.haplotype_name,
           g.ome, g.version, g.genomeSource, g.citation, g.assembly_accession
    FROM joined_ships j
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN taxonomy t ON j.tax_id = t.id
    LEFT JOIN starship_features sf ON a.id = sf.accession_id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN navis_names n ON j.ship_navis_id = n.id
    LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    """

    params = []
    if curated:
        meta_query += " WHERE j.curated_status = 'curated'"

    if accession_tags:
        where_clause = " WHERE " if not curated else " AND "
        if isinstance(accession_tags, list):
            # Use ? for SQLite placeholders
            placeholders = ",".join(["?"] * len(accession_tags))
            meta_query += f"{where_clause}a.accession_tag IN ({placeholders})"
            params = tuple(accession_tags)
        else:
            # Use ? for SQLite placeholder
            meta_query += f"{where_clause}a.accession_tag = ?"
            params = (accession_tags,)

    try:
        if params:
            meta_df = pd.read_sql_query(meta_query, session.bind, params=params)
        else:
            meta_df = pd.read_sql_query(meta_query, session.bind)

        return meta_df
    except Exception as e:
        logger.error(f"Error fetching meta data: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_paper_data():
    """Fetch paper data from the database and cache the result."""
    session = StarbaseSession()

    paper_query = """
    SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, 
           p.shortCitation, f.familyName, f.type_element_reference
    FROM papers p
    LEFT JOIN family_names f ON p.shortCitation = f.type_element_reference
    """
    try:
        paper_df = pd.read_sql_query(paper_query, session.bind)
        if paper_df.empty:
            logger.warning("Fetched paper DataFrame is empty.")
        return paper_df
    except Exception as e:
        logger.error(f"Error fetching paper data: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    session = StarbaseSession()

    query = """
    SELECT a.accession_tag, a.version_tag, 
           CASE 
               WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
               THEN a.accession_tag || '.' || a.version_tag
               ELSE a.accession_tag
           END as accession_display,
           f.familyName, p.shortCitation, t.`order`, t.family, t.name 
    FROM joined_ships j
    LEFT JOIN taxonomy t ON j.tax_id = t.id
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN papers p ON f.type_element_reference = p.shortCitation
    -- Only show entries that have sequences
    INNER JOIN ships s ON s.accession_id = a.id
    WHERE 1=1
    """

    if curated:
        query += " AND j.curated_status = 'curated'"

    try:
        df = pd.read_sql_query(query.strip(), session.bind)

        if dereplicate:
            df = df.drop_duplicates(subset="accession_tag")

        if df.empty:
            logger.warning("Fetched Download DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching download data: {str(e)}")
        raise
    finally:
        session.close()


@db_retry_decorator()
def fetch_ships(
    accession_tags=None, curated=False, dereplicate=True, with_sequence=False
):
    """
    Fetch ship data for specified accession tags.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all ships.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing ship data
    """
    import re

    session = StarbaseSession()

    base_query = """
    WITH valid_ships AS (
        SELECT DISTINCT 
            a.id as accession_id, 
            a.accession_tag, a.version_tag,
            CASE 
                WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
                THEN a.accession_tag || '.' || a.version_tag
                ELSE a.accession_tag
            END as accession_display,
            j.curated_status,
            sf.elementBegin, sf.elementEnd, sf.contigID,
            t.name, t.family, t.`order`,
            f.familyName, n.navis_name, h.haplotype_name,
            g.assembly_accession, c.captainID"""

    if dereplicate:
        base_query += """,
            ROW_NUMBER() OVER (
                PARTITION BY a.accession_tag 
                ORDER BY CASE 
                    WHEN a.version_tag IS NULL OR a.version_tag = '' THEN 0 
                    ELSE CAST(a.version_tag AS INTEGER) 
                END DESC
            ) as rn"""

    base_query += """
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        LEFT JOIN starship_features sf ON a.id = sf.accession_id
        LEFT JOIN captains c ON j.captain_id = c.id
        WHERE 1=1
    """

    query = base_query

    use_accessions = []
    if accession_tags:
        for accession in accession_tags:
            if re.match("\..*", accession):
                use_accessions.append(re.sub(pattern="\..*", repl="", string=accession))
            else:
                use_accessions.append(accession)

        query += " AND a.accession_tag IN ({})".format(
            ",".join(f"'{tag}'" for tag in accession_tags)
        )
    if curated:
        query += " AND j.curated_status = 'curated'"

    if with_sequence:
        query += """
        )
        SELECT 
            v.accession_id,
            v.accession_tag,
            v.version_tag,
            v.accession_display,
            v.curated_status,
            v.elementBegin,
            v.elementEnd,
            v.contigID,
            v.name,
            v.family,
            v.`order`,
            v.familyName,
            v.navis_name,
            v.haplotype_name,
            v.assembly_accession,
            s.sequence,
            s.md5,
            s.rev_comp_md5,
            v.captainID
        FROM valid_ships v
        LEFT JOIN ships s ON s.accession_id = v.accession_id
        WHERE s.sequence IS NOT NULL"""

        if dereplicate:
            query += " AND v.rn = 1"

        query += """
        """
    else:
        query += """
        )
        SELECT 
            v.accession_id,
            v.accession_tag,
            v.version_tag,
            v.accession_display,
            v.curated_status,
            v.elementBegin,
            v.elementEnd,
            v.contigID,
            v.name,
            v.family,
            v.`order`,
            v.familyName,
            v.navis_name,
            v.haplotype_name,
            v.assembly_accession,
            v.captainID
        FROM valid_ships v"""

        if dereplicate:
            query += " WHERE v.rn = 1"

        query += """
        """

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched ships DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching ships data: {str(e)}")
        raise
    finally:
        session.close()


@db_retry_decorator()
def fetch_ship_table(curated=False):
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    session = StarbaseSession()

    query = """
    SELECT DISTINCT 
        a.accession_tag, a.version_tag,
        CASE 
            WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
            THEN a.accession_tag || '.' || a.version_tag
            ELSE a.accession_tag
        END as accession_display,
        f.familyName,
        t.name
    FROM joined_ships js
    INNER JOIN accessions a ON js.ship_id = a.id
    LEFT JOIN taxonomy t ON js.tax_id = t.id
    LEFT JOIN family_names f ON js.ship_family_id = f.id
    -- Filter for ships that have sequence data
    INNER JOIN ships s ON s.accession_id = a.id AND s.sequence IS NOT NULL
    -- Filter for ships that have GFF annotation data
    INNER JOIN gff g ON g.ship_id = a.id
    """

    if curated:
        query += " WHERE js.curated_status = 'curated'"

    query += " ORDER BY f.familyName ASC"

    try:
        df = pd.read_sql_query(query, session.bind)
        return df
    except Exception as e:
        logger.error(f"Error fetching ship table data: {str(e)}")
        raise

@db_retry_decorator()
def fetch_accession_ship(accession_tag):
    """Fetch sequence and GFF data for a specific ship."""
    session = StarbaseSession()

    sequence_query = """
    SELECT s.sequence
    FROM ships s
    LEFT JOIN accessions a ON s.accession_id = a.id
    WHERE a.accession_tag = :accession_tag
    """

    gff_query = """
    SELECT g.*
    FROM gff g
    LEFT JOIN accessions a ON g.ship_id = a.id
    WHERE a.accession_tag = :accession_tag
    """

    try:
        sequence = pd.read_sql_query(
            sequence_query, session.bind, params={"accession_tag": accession_tag}
        )
        gff_df = pd.read_sql_query(
            gff_query, session.bind, params={"accession_tag": accession_tag}
        )

        # Get the sequence string
        sequence_str = sequence.iloc[0]["sequence"] if not sequence.empty else None

        # Ensure GFF data is a DataFrame
        if gff_df.empty:
            logger.warning(f"No GFF data found for accession: {accession_tag}")
            gff_df = None

        logger.debug(f"GFF data type: {type(gff_df)}")
        logger.debug(f"GFF columns: {gff_df.columns if gff_df is not None else 'None'}")

        return {"sequence": sequence_str, "gff": gff_df}
    except Exception as e:
        logger.error(f"Error fetching ship data for {accession_tag}: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_captains(
    accession_tags=None, curated=False, dereplicate=True, with_sequence=False
):
    """
    Fetch captain data for specified accession tags.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all captains.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing captain data
    """
    session = StarbaseSession()

    query = """
    WITH valid_captains AS (
        SELECT DISTINCT 
            a.id, 
            a.accession_tag,
            a.version_tag,
            CASE 
                WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
                THEN a.accession_tag || '.' || a.version_tag
                ELSE a.accession_tag
            END as accession_display,
            j.curated_status,
            j.starshipID,
            sf.captainID,
            c."sequence",
            n.navis_name,
            h.haplotype_name,
            c.captainID
        FROM joined_ships j
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        LEFT JOIN captains c ON j.captain_id = c.id
        LEFT JOIN starship_features sf ON a.id = sf.accession_id
        WHERE 1=1
    """

    if accession_tags:
        query += " AND a.accession_tag IN ({})".format(
            ",".join(f"'{tag}'" for tag in accession_tags)
        )
    if curated:
        query += " AND j.curated_status = 'curated'"

    if with_sequence:
        query += """
        )
        SELECT 
            v.id,
            v.accession_tag,
            v.version_tag,
            v.accession_display,
            v.curated_status,
            v.starshipID,
            v.captainID,
            v.sequence,
            v.navis_name,
            v.haplotype_name,
            v.captainID
        FROM valid_captains v
        WHERE v.sequence IS NOT NULL
        """
    else:
        query += """
        )
        SELECT 
            v.id,
            v.accession_tag,
            v.version_tag,
            v.accession_display,
            v.curated_status,
            v.starshipID,
            v.captainID,
            v.navis_name,
            v.haplotype_name,
            v.captainID
        FROM valid_captains v
        """

    try:
        df = pd.read_sql_query(query, session.bind)

        if dereplicate:
            df = df.drop_duplicates(subset="accession_tag")

        if df.empty:
            logger.warning("Fetched captains DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching captains data: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_captain_tree():
    fallback_tree_path = PHYLOGENY_PATHS["tree"]

    with open(fallback_tree_path, "r") as f:
        return f.read()

@db_retry_decorator()
def fetch_sf_data():
    sf_data = pd.read_csv(PHYLOGENY_PATHS["clades"], sep="\t")

    # Add debug logging
    logger.debug(f"Loaded sf_data columns: {sf_data.columns.tolist()}")
    logger.debug(f"Loaded sf_data head: \n{sf_data.head()}")

    return sf_data


@db_retry_decorator()
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
    try:
        # use metadata from previous query
        meta_df = fetch_meta_data(curated=False)
        total_count = len(meta_df)
        curated_count = len(meta_df[meta_df["curated_status"] == "curated"])
        uncurated_count = total_count - curated_count
        species_count = len(meta_df["name"].dropna().unique())
        family_count = len(meta_df["familyName"].dropna().unique())

        stats = {
            "total_starships": total_count,
            "curated_starships": curated_count,
            "uncurated_starships": uncurated_count,
            "species_count": species_count,
            "family_count": family_count,
        }
        return stats
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        raise
    finally:
        session.close()

import logging
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

logger = logging.getLogger(__name__)


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


@cache.memoize()
def fetch_meta_data(curated=False, accession_tag=None):
    """
    Fetch metadata from the database with caching.

    Args:
        curated (bool): If True, only return curated entries
        accession_tag (str or list): Single accession tag or list of accession tags

    Returns:
        pd.DataFrame: Metadata for the specified accession tags
    """
    session = StarbaseSession()

    meta_query = """
    SELECT j.ship_family_id, j.curated_status, t.taxID, j.starshipID,
           j.ome, j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, j.contigID, 
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.name, 
           g.version, g.genomeSource, g.citation, a.accession_tag, j.strain, j.starship_navis, j.starship_haplotype, g.assembly_accession
    FROM joined_ships j
    INNER JOIN taxonomy t ON j.taxid = t.id
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    """

    params = []
    if curated:
        meta_query += " WHERE j.curated_status = 'curated'"

    if accession_tag:
        where_clause = " WHERE " if not curated else " AND "
        if isinstance(accession_tag, list):
            # Use ? for SQLite placeholders
            placeholders = ",".join(["?"] * len(accession_tag))
            meta_query += f"{where_clause}a.accession_tag IN ({placeholders})"
            params = tuple(accession_tag)
        else:
            # Use ? for SQLite placeholder
            meta_query += f"{where_clause}a.accession_tag = ?"
            params = (accession_tag,)

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


@cache.memoize()
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


@cache.memoize()
def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    session = StarbaseSession()

    query = """
    SELECT a.accession_tag, f.familyName, p.shortCitation, t.`order`, t.family, t.name 
    FROM joined_ships j
    INNER JOIN taxonomy t ON j.taxid = t.id
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN papers p ON f.type_element_reference = p.shortCitation
    -- Only show entries that have sequences
    INNER JOIN ships s ON s.accession = a.id
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


@cache.memoize()
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
            t.name,
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
    LEFT JOIN ships s ON s.accession = v.accession_id
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


@cache.memoize()
@db_retry_decorator()
def fetch_ship_table(curated=False):
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    session = StarbaseSession()

    query = """
    SELECT DISTINCT 
        a.accession_tag,
        f.familyName,
        t.name
    FROM joined_ships js
    LEFT JOIN accessions a ON js.ship_id = a.id
    LEFT JOIN taxonomy t ON js.taxid = t.id
    LEFT JOIN family_names f ON js.ship_family_id = f.id
    -- Filter for ships that have sequence data
    LEFT JOIN ships s ON s.accession = a.id AND s.sequence IS NOT NULL
    LEFT JOIN gff g ON g.ship_id = a.id
    WHERE js.orphan IS NULL
    """

    if curated:
        query += " AND js.curated_status = 'curated'"

    query += " ORDER BY f.familyName ASC"

    try:
        df = pd.read_sql_query(query, session.bind)
        return df
    except Exception as e:
        logger.error(f"Error fetching ship table data: {str(e)}")
        raise


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
def fetch_accession_ship(accession_tag):
    """Fetch sequence and GFF data for a specific ship."""
    session = StarbaseSession()

    sequence_query = """
    SELECT s.sequence
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
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
def fetch_captains(accession_tags=None, curated=False, dereplicate=True):
    """
    Fetch captain data for specified accession tags.

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
            c.*
        FROM captains c
        INNER JOIN joined_ships j ON c.id = j.captainID_new
        INNER JOIN accessions a ON j.ship_id = a.id
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
    FROM valid_captains v
    LEFT JOIN ships s ON s.accession = v.accession_id
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


@cache.memoize()
def fetch_captain_tree():
    fallback_tree_path = PHYLOGENY_PATHS["tree"]

    with open(fallback_tree_path, "r") as f:
        return f.read()


@cache.memoize()
def fetch_sf_data():
    sf_data = pd.read_csv(PHYLOGENY_PATHS["clades"], sep="\t")

    # Add debug logging
    logger.debug(f"Loaded sf_data columns: {sf_data.columns.tolist()}")
    logger.debug(f"Loaded sf_data head: \n{sf_data.head()}")

    return sf_data


@cache.memoize()
@db_retry_decorator()
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
    try:
        # Get curated and uncurated counts
        curated_count = (
            session.execute("""
            SELECT COUNT(DISTINCT a.accession_tag) 
            FROM accessions a
            LEFT JOIN joined_ships j ON j.ship_id = a.id
            WHERE j.curated_status = 'curated'
        """).scalar()
            or 0
        )

        uncurated_count = (
            session.execute("""
            SELECT COUNT(DISTINCT a.accession_tag) 
            FROM accessions a
            LEFT JOIN joined_ships j ON j.ship_id = a.id
            WHERE j.curated_status != 'curated' OR j.curated_status IS NULL
        """).scalar()
            or 0
        )

        stats = {
            "curated_starships": curated_count,
            "uncurated_starships": uncurated_count,
            "species_count": session.execute(
                "SELECT COUNT(DISTINCT name) FROM taxonomy"
            ).scalar()
            or 0,
            "family_count": session.execute(
                "SELECT COUNT(DISTINCT newFamilyID) FROM family_names WHERE newFamilyID IS NOT NULL"
            ).scalar()
            or 0,
        }
        return stats
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        raise
    finally:
        session.close()

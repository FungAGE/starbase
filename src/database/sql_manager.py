from typing import Any
import logging
import pandas as pd
from src.config.database import StarbaseSession
from src.utils.plot_utils import create_sunburst_plot
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
from contextlib import contextmanager
import sqlalchemy.exc
from sqlalchemy import text
from src.config.cache import cache
logger = logging.getLogger(__name__)

@cache.memoize(timeout=86400)
def fetch_meta_data(curated=False):
    """Fetch metadata from the database with caching."""
    session = StarbaseSession()
    
    meta_query = """
    SELECT j.ship_family_id, j.curated_status, t.taxID, j.starshipID,
           j.ome, j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, j.contigID, 
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.genus, t.species, 
           g.version, g.genomeSource, g.citation, a.accession_tag, g.strain, j.starship_navis, j.starship_haplotype, g.assembly_accession
    FROM joined_ships j
    INNER JOIN taxonomy t ON j.taxid = t.id
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    """

    if curated:
        meta_query += " AND j.curated_status = 'curated'"

    try:
        meta_df = pd.read_sql_query(meta_query, session.bind)
        return meta_df
    except Exception as e:
        logger.error(f"Error fetching meta data: {str(e)}")
        raise
    finally:
        session.close()

@cache.memoize(timeout=86400)
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

@cache.memoize(timeout=86400)
def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    session = StarbaseSession()

    query = """
    SELECT a.accession_tag, f.familyName, p.shortCitation, t.`order`, t.family, t.species 
    FROM joined_ships j
    INNER JOIN taxonomy t ON j.taxid = t.id
    INNER JOIN accessions a ON j.ship_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN papers p ON f.type_element_reference = p.shortCitation
    """
    
    if curated:
        query += " AND j.curated_status = 'curated'"
    
    try:
        df = pd.read_sql_query(query.strip(), session.bind)
        
        if dereplicate:
            df = df.drop_duplicates(subset='accession_tag')
            logger.info(f"Dereplicated to {len(df)} unique accession tags.")
        
        if df.empty:
            logger.warning("Fetched Download DataFrame is empty.")            
        return df
    except Exception as e:
        logger.error(f"Error fetching download data: {str(e)}")
        raise
    finally:
        session.close()

@cache.memoize(timeout=86400)
def fetch_all_ships(curated=True):
    session = StarbaseSession()

    query = """
    SELECT s.*, a.accession_tag
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
    LEFT JOIN joined_ships j ON j.ship_id = a.id
    WHERE 1=1
    """
    
    if curated:
        query += " AND j.curated_status = 'curated'"

    session = StarbaseSession()
    try:
        df = pd.read_sql_query(query, session.bind)
        if df.empty:
            logger.warning("Fetched all_ships DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching all_ships data: {str(e)}")
        raise
    finally:
        session.close()
    

@cache.memoize(timeout=86400)
def fetch_ship_table(curated=False):
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    session = StarbaseSession()
    
    query = """
    SELECT DISTINCT 
        a.accession_tag,
        f.familyName,
        t.species
    FROM joined_ships js
    LEFT JOIN accessions a ON js.ship_id = a.id
    LEFT JOIN taxonomy t ON js.taxid = t.id
    LEFT JOIN family_names f ON js.ship_family_id = f.id
    -- Filter for ships that have sequence data
    LEFT JOIN ships s ON s.accession = a.id AND s.sequence IS NOT NULL
    -- Filter for ships that have GFF data
    LEFT JOIN gff g ON g.ship_id = a.id
    WHERE js.orphan IS NULL
    """
    
    if curated:
        query += " AND js.curated_status = 'curated'"

    query += " ORDER BY f.familyName ASC"

@contextmanager
def db_session_manager():
    """Context manager for database sessions with timeout"""
    session = None
    try:
        session = StarbaseSession()
        # SQLite doesn't support SET SESSION, so we'll skip the timeout setting
        if session.bind.dialect.name != 'sqlite':
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
        )
    )

@db_retry_decorator()
def fetch_ship_table(curated=True):
    with db_session_manager() as session:
        try:
            query = """
            SELECT DISTINCT 
                a.accession_tag,
                f.familyName,
                t.species
            FROM joined_ships js
            INNER JOIN accessions a ON js.ship_id = a.id
            INNER JOIN family_names f ON js.ship_family_id = f.id
            INNER JOIN ships s ON s.accession = a.id
            INNER JOIN gff g ON g.ship_id = a.id
            INNER JOIN taxonomy t ON js.taxid = t.id
            """
            
            if curated:
                query += " AND js.curated_status = 'curated'"

            query += " ORDER BY f.familyName ASC LIMIT 1000"

            df = pd.read_sql_query(query, session.bind)
            return df
        except Exception as e:
            logger.error(f"Error fetching ship_table data: {str(e)}")
            return None

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
        sequence = pd.read_sql_query(sequence_query, session.bind, params={"accession_tag": accession_tag})
        gff_df = pd.read_sql_query(gff_query, session.bind, params={"accession_tag": accession_tag})
        
        # Get the sequence string
        sequence_str = sequence.iloc[0]["sequence"] if not sequence.empty else None
        
        # Ensure GFF data is a DataFrame
        if gff_df.empty:
            logger.warning(f"No GFF data found for accession: {accession_tag}")
            gff_df = None
            
        logger.debug(f"GFF data type: {type(gff_df)}")
        logger.debug(f"GFF columns: {gff_df.columns if gff_df is not None else 'None'}")
        
        return {
            "sequence": sequence_str,
            "gff": gff_df
        }
    except Exception as e:
        logger.error(f"Error fetching ship data for {accession_tag}: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_all_captains():
    session = StarbaseSession()

    query = f"""
    SELECT c.*
    FROM captains c
    """

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched captain DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching captain data: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_captain_tree():
    session = StarbaseSession()

    tree_query = """SELECT string FROM trees WHERE id=1"""

    try:
        tree_string = session.execute(tree_query).fetchone()

        if tree_string.empty:
            logger.warning("Fetched tree string is empty.")
        return tree_string
    except Exception as e:
        logger.error(f"Error fetching tree string: {str(e)}")
        raise
    finally:
        session.close()

@db_retry_decorator()
def fetch_sf_data():
    session = StarbaseSession()

    query = """
    SELECT sf.*
    FROM superfam-clades sf
    """

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched sf DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching sf data: {str(e)}")
        raise
    finally:
        session.close()

@cache.memoize(timeout=86400)
@db_retry_decorator()
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
    try:
        # Get curated and uncurated counts
        curated_count = session.execute("""
            SELECT COUNT(DISTINCT a.accession_tag) 
            FROM accessions a
            LEFT JOIN joined_ships j ON j.ship_id = a.id
            WHERE j.curated_status = 'curated'
        """).scalar() or 0
        
        uncurated_count = session.execute("""
            SELECT COUNT(DISTINCT a.accession_tag) 
            FROM accessions a
            LEFT JOIN joined_ships j ON j.ship_id = a.id
            WHERE j.curated_status != 'curated' OR j.curated_status IS NULL
        """).scalar() or 0
        
        stats = {
            "curated_starships": curated_count,
            "uncurated_starships": uncurated_count,
            "species_count": session.execute(
                "SELECT COUNT(DISTINCT species) FROM taxonomy"
            ).scalar() or 0,
            "family_count": session.execute(
                "SELECT COUNT(DISTINCT newFamilyID) FROM family_names WHERE newFamilyID IS NOT NULL"
            ).scalar() or 0
        }
        return stats
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        raise
    finally:
        session.close()
import logging
import pandas as pd
from src.components.cache_manager import save_to_cache, load_from_cache, cache_exists
from src.components.sql_engine import starbase_session_factory
from src.utils.plot_utils import create_sunburst_plot

logger = logging.getLogger(__name__)


def generate_cache_key(base_key, unique_identifier=None):
    """Generate a cache key by combining a base key with a unique identifier."""
    if unique_identifier:
        return f"{base_key}_{unique_identifier}"
    else:
        return base_key


def fetch_meta_data(curated=False):
    """Fetch metadata from the database and cache the result."""
    cache_key = generate_cache_key("meta_data")

    if cache_exists(cache_key):
        logger.info(f"Meta data found in cache for key '{cache_key}'")
        return load_from_cache(cache_key)

    logger.info(f"Fetching meta data from database for key '{cache_key}'")

    meta_query = """
    SELECT j.ship_family_id, j.curated_status, t.taxID, j.starshipID,
           j.ome, j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, j.contigID, 
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.species, 
           g.version, g.genomeSource, g.citation, a.accession_tag, g.strain
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN genomes g ON j.genome_id = g.id
    JOIN accessions a ON j.ship_id = a.id
    WHERE j.orphan IS NULL
    """

    session = starbase_session_factory()

    if curated:
        meta_query += " AND j.curated_status = 'curated'"

    try:
        meta_df = pd.read_sql_query(meta_query, session.bind)
        logger.info(
            f"Meta data successfully fetched from database, caching it under key '{cache_key}'"
        )
        save_to_cache(meta_df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching meta data: {str(e)}")
        raise

    return meta_df


def fetch_paper_data():
    """Fetch paper data from the database and cache the result."""
    cache_key = generate_cache_key("paper_data")

    if cache_exists(cache_key):
        logger.info(f"Paper data found in cache for key '{cache_key}'")
        return load_from_cache(cache_key)

    logger.info(f"Fetching paper data from database for key '{cache_key}'")
    paper_query = """
    SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, 
           p.shortCitation, f.familyName, f.type_element_reference
    FROM papers p
    LEFT JOIN family_names f ON p.shortCitation = f.type_element_reference
    """

    session = starbase_session_factory()

    try:
        # Use session.bind to execute the query through the starbase_engine
        paper_df = pd.read_sql_query(paper_query, session.bind)

        if paper_df.empty:
            logger.warning("Fetched paper DataFrame is empty.")
        else:
            logger.info(
                f"Paper data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(paper_df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching paper data: {str(e)}")
        return None
    finally:
        # Close the session instance instead of the session factory
        session.close()

    return paper_df


def cache_sunburst_plot(family, df):
    """Create sunburst plots for wiki page and cache the result."""
    cache_key = generate_cache_key("sunburst", family)

    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    sunburst = create_sunburst_plot(df=df, type="tax", title_switch=False)
    save_to_cache(sunburst, cache_key)

    return sunburst


def fetch_download_data():
    cache_key = generate_cache_key("download_data")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT a.accession_tag, f.familyName, t.`order`, t.family, t.species 
    FROM accessions a
    LEFT JOIN joined_ships j ON a.id = j.ship_id
    LEFT JOIN taxonomy t ON j.taxid = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    WHERE j.orphan IS NULL
    """

    session = starbase_session_factory()

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched Download DataFrame is empty.")
        else:
            logger.info(
                f"Download data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching download data: {str(e)}")
        return None
    finally:
        session.close()
    return df


def fetch_all_ships():
    cache_key = generate_cache_key("all_ships")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT s.*, a.accession_tag
    FROM ships s
    JOIN accessions a ON s.accession = a.id
    """
    session = starbase_session_factory()
    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched all_ships DataFrame is empty.")
        else:
            logger.info(
                f"all_ships data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching all_ships data: {str(e)}")
        return None
    finally:
        session.close()
    return df


def fetch_accession_gff(accession):
    cache_key = generate_cache_key("accession_gff", accession)

    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT g.*
    FROM gff g
    LEFT JOIN accessions a ON g.ship_id = a.id
    LEFT JOIN ships s on s.accession = a.id
    LEFT JOIN joined_ships j ON j.ship_id = a.id
    LEFT JOIN taxonomy t ON j.taxid = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    WHERE g.accession = ? AND j.orphan IS NULL
    """

    session = starbase_session_factory()
    try:
        df = pd.read_sql_query(query, session.bind, params=(accession,))

        if df.empty:
            logger.warning("Fetched gff DataFrame is empty.")
        else:
            logger.info(
                f"gff data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching gff data: {str(e)}")
        return None
    finally:
        session.close()
    return df


def fetch_ship_table(meta_df=None):
    """Fetch and filter ship table data based on metadata DataFrame."""
    cache_key = generate_cache_key("ship_table")

    if cache_exists(cache_key):
        logger.info(f"Ship table found in cache for key '{cache_key}'")
        return load_from_cache(cache_key)

    query = """
    SELECT DISTINCT a.accession_tag, f.familyName, t.species
    FROM gff g
    JOIN accessions a ON g.ship_id = a.id
    JOIN joined_ships js ON a.id = js.ship_id 
    JOIN taxonomy t ON js.taxid = t.id
    JOIN family_names f ON js.ship_family_id = f.id
    JOIN ships s ON s.accession = a.id
    WHERE s.sequence IS NOT NULL
    AND g.ship_id IS NOT NULL
    AND js.orphan IS NULL
    """

    session = starbase_session_factory()
    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched ship_table DataFrame is empty.")
        else:
            df = df.sort_values(by="familyName", ascending=True)

            logger.info(
                f"ship_table data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching ship_table data: {str(e)}")
        return None
    finally:
        session.close()
    return df


def fetch_all_captains():
    cache_key = generate_cache_key("all_captains")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = f"""
    SELECT c.*
    FROM captains c
    """

    session = starbase_session_factory()

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched captain DataFrame is empty.")
        else:
            logger.info(
                f"captain data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching captain data: {str(e)}")
        return None
    finally:
        session.close()
    return df


def fetch_captain_tree():
    cache_key = generate_cache_key("captain_tree")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    tree_query = """SELECT string FROM trees WHERE id=1"""

    session = starbase_session_factory()

    try:
        tree_string = session.execute(tree_query).fetchone()

        if tree_string.empty:
            logger.warning("Fetched tree string is empty.")
        else:
            logger.info(
                f"tree string successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(tree_string, cache_key)
    except Exception as e:
        logger.error(f"Error fetching tree string: {str(e)}")
        return None
    finally:
        session.close()
    return tree_string


def fetch_sf_data():
    cache_key = generate_cache_key("sf_data")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT sf.*
    FROM superfam-clades sf
    """
    session = starbase_session_factory()
    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched sf DataFrame is empty.")
        else:
            logger.info(
                f"sf data successfully fetched from database, caching it under key '{cache_key}'"
            )
            save_to_cache(df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching sf data: {str(e)}")
        return None
    finally:
        session.close()
    return df

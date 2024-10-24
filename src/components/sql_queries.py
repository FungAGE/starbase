import logging
import pandas as pd
from src.components.cache_manager import save_to_cache, load_from_cache, cache_exists
from src.components.mariadb import engine
from src.utils.plot_utils import create_sunburst_plot

logger = logging.getLogger(__name__)


def generate_cache_key(base_key, unique_identifier=None):
    """Generate a cache key by combining a base key with a unique identifier."""
    if unique_identifier:
        return f"{base_key}_{unique_identifier}"
    else:
        return base_key


def fetch_meta_data(curated=True):
    """Fetch metadata from the database and cache the result."""
    cache_key = generate_cache_key("meta_data")

    if cache_exists(cache_key):
        logger.info(f"Meta data found in cache for key '{cache_key}'")
        return load_from_cache(cache_key)

    logger.info(f"Fetching meta data from database for key '{cache_key}'")

    meta_query = """
    SELECT j.ship_family_id, j.curated_status, j.taxid, j.ship_id, j.genome_id, 
           j.ome, j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, j.contigID, 
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.species, 
           g.version, g.genomeSource, g.citation, a.accession_tag
    FROM joined_ships j
    LEFT JOIN taxonomy t ON j.taxid = t.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN accessions a ON j.ship_id = a.id
    WHERE j.orphan IS NULL
    """

    if curated:
        meta_query += " AND j.curated_status = 'curated'"

    try:
        meta_df = pd.read_sql_query(meta_query, engine)
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

    try:
        paper_df = pd.read_sql_query(paper_query, engine)
        if paper_df.empty:
            logger.warning("Fetched paper DataFrame is empty.")
        else:
            logger.info(
                f"Paper data successfully fetched from database, caching it under key '{cache_key}'"
            )
        save_to_cache(paper_df, cache_key)
    except Exception as e:
        logger.error(f"Error fetching paper data: {str(e)}")
        return None  # Ensure it returns None on error

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
    df = pd.read_sql_query(query, engine)

    save_to_cache(df, cache_key)

    return df


def fetch_all_ships():
    cache_key = generate_cache_key("all_ships")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT s.*, a.accession_tag
    FROM ships s
    LEFT JOIN accessions a ON s.accession = a.id
    """
    df = pd.read_sql_query(query, engine)

    save_to_cache(df, cache_key)
    return df


def fetch_accession_gff(accession):
    cache_key = generate_cache_key("accession_gff", accession)

    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT g.*
    FROM gff g
    WHERE g.accession = %s
    """

    df = pd.read_sql_query(query, engine, params=(accession,))

    save_to_cache(df, cache_key)

    return df


def fetch_ship_table(meta_df=None):
    """Fetch and filter ship table data based on metadata DataFrame."""
    cache_key = generate_cache_key("ship_table")

    if cache_exists(cache_key):
        logger.info(f"Ship table found in cache for key '{cache_key}'")
        return load_from_cache(cache_key)

    if meta_df is None:
        meta_df = fetch_meta_data(curated=True)

    filtered_df = meta_df[
        (meta_df["curated_status"] == "curated")
        & (meta_df["familyName"].notna())
        & (meta_df["accession_tag"].notna())
    ]

    ship_table_df = filtered_df[
        ["accession_tag", "familyName", "species"]
    ].drop_duplicates()
    ship_table_df = ship_table_df.sort_values(by="familyName", ascending=True)

    save_to_cache(ship_table_df, cache_key)
    logger.info(f"Ship table successfully created and cached under key '{cache_key}'")

    return ship_table_df


def fetch_all_captains():
    cache_key = generate_cache_key("all_captains")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = f"""
    SELECT c.*
    FROM captains c
    """


def fetch_captain_tree():
    cache_key = generate_cache_key("captain_tree")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    tree_query = """SELECT string FROM trees WHERE id=1"""
    with engine.connect() as connection:
        tree_string = connection.execute(tree_query).fetchone()
    save_to_cache(tree_string, cache_key)
    return tree_string


def fetch_sf_data():
    cache_key = generate_cache_key("sf_data")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    sf_query = """
    SELECT sf.*
    FROM superfam-clades sf
    """
    df = pd.read_sql_query(sf_query)
    save_to_cache(df, cache_key)
    return df

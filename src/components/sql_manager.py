import os
import pickle
from typing import Any
import logging
import pandas as pd
from src.components.cache import cache
from src.components.sql_engine import starbase_session_factory
from src.utils.plot_utils import create_sunburst_plot

logger = logging.getLogger(__name__)

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
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.genus, t.species, 
           g.version, g.genomeSource, g.citation, a.accession_tag, g.strain, j.starship_navis, j.starship_haplotype
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


def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    cache_key = generate_cache_key(f"download_data_curated_{curated}_derep_{dereplicate}")
    
    if cache_exists(cache_key):
        return load_from_cache(cache_key)
    
    query = """
    SELECT a.accession_tag, f.familyName, t.`order`, t.family, t.species 
    FROM joined_ships j
    JOIN taxonomy t ON j.taxid = t.id
    JOIN family_names f ON j.ship_family_id = f.id
    JOIN genomes g ON j.genome_id = g.id
    JOIN accessions a ON j.ship_id = a.id
    WHERE j.orphan IS NULL
    """
    
    if curated:
        query += " AND j.curated_status = 'curated'"
    
    session = starbase_session_factory()
    try:
        df = pd.read_sql_query(query.strip(), session.bind)
        
        if dereplicate:
            df = df.drop_duplicates(subset='accession_tag')
            logger.info(f"Dereplicated to {len(df)} unique accession tags.")
        
        if df.empty:
            logger.warning("Fetched Download DataFrame is empty.")
        else:
            logger.debug(f"Attempting to cache download data with key: {cache_key}")
            logger.debug(f"Cache status before save: exists={cache_exists(cache_key)}")
            save_to_cache(df, cache_key)
            logger.debug(f"Cache status after save: exists={cache_exists(cache_key)}")
            
        return df
    except Exception as e:
        logger.error(f"Error fetching download data: {str(e)}")
        return None
    finally:
        session.close()

def fetch_all_ships(curated=True):
    cache_key = generate_cache_key("all_ships", f"curated_{curated}")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    query = """
    SELECT s.*, a.accession_tag
    FROM ships s
    JOIN accessions a ON s.accession = a.id
    JOIN joined_ships j ON j.ship_id = a.id
    WHERE 1=1
    """
    
    if curated:
        query += " AND j.curated_status = 'curated'"

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

def get_database_stats():
    """Get statistics about the Starship database."""
    cache_key = generate_cache_key("stats")
    if cache_exists(cache_key):
        return load_from_cache(cache_key)

    session = starbase_session_factory()
    try:
        # Get curated and uncurated counts
        curated_count = session.execute("""
            SELECT COUNT(*) 
            FROM accessions a
            JOIN joined_ships j ON j.ship_id = a.id
            WHERE j.curated_status = 'curated'
        """).scalar() or 0
        
        uncurated_count = session.execute("""
            SELECT COUNT(*) 
            FROM accessions a
            JOIN joined_ships j ON j.ship_id = a.id
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
        save_to_cache(stats, cache_key)
        return stats
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        return {
            "curated_starships": 0,
            "uncurated_starships": 0,
            "species_count": 0,
            "family_count": 0
        }
    finally:
        session.close()

#########################
# Cache functions
#########################
CACHE_DIR = ".cache"  # Directory where all cached objects will be saved

def generate_cache_key(base_key, unique_identifier=None):
    """Generate a cache key by combining a base key with a unique identifier."""
    if unique_identifier:
        return f"{base_key}_{unique_identifier}"
    else:
        return base_key


if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR, exist_ok=True)
    logger.info(f"Cache directory created at {CACHE_DIR}")


def get_cache_filepath(cache_key: str) -> str:
    """Generate the full path for a cache file based on a cache key."""
    filepath = os.path.join(CACHE_DIR, f"{cache_key}.pkl")
    logger.debug(f"Cache file path for key '{cache_key}': {filepath}")
    return filepath


def save_to_cache(obj: Any, cache_key: str):
    """Save an object to a cache file."""
    filepath = get_cache_filepath(cache_key)
    try:
        with open(filepath, "wb") as f:
            pickle.dump(obj, f)
        logger.info(f"Cached object under key '{cache_key}' at {filepath}")
    except Exception as e:
        logger.error(f"Failed to save cache for key '{cache_key}': {str(e)}")


def load_from_cache(cache_key: str) -> Any:
    """Load an object from the cache."""
    filepath = get_cache_filepath(cache_key)
    if os.path.exists(filepath):
        try:
            with open(filepath, "rb") as f:
                logger.info(f"Loaded cached object for key '{cache_key}'")
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Failed to load cache for key '{cache_key}': {str(e)}")
    else:
        logger.warning(f"No cache file found for key '{cache_key}'")
    return None


def cache_exists(cache_key: str) -> bool:
    """Check if a cached file exists for a given cache key."""
    exists = os.path.exists(get_cache_filepath(cache_key))
    if exists:
        logger.debug(f"Cache exists for key '{cache_key}'")
    else:
        logger.debug(f"Cache does not exist for key '{cache_key}'")
    return exists


def invalidate_cache(cache_key: str):
    """Delete the cached file for a specific cache key."""
    filepath = get_cache_filepath(cache_key)
    if os.path.exists(filepath):
        os.remove(filepath)
        logger.info(
            f"Cache invalidated for key '{cache_key}' (file deleted: {filepath})"
        )
    else:
        logger.warning(f"No cache file found to invalidate for key '{cache_key}'")

def initialize_cache() -> None:
    """Initialize cache with commonly used data during app startup."""
    try:
        # Cache metadata
        meta_data = fetch_meta_data()
        cache.set("meta_data", meta_data)
        
        # Cache paper data
        paper_data = fetch_paper_data()
        cache.set("paper_data", paper_data)
        
        # Pre-calculate and cache filtered options
        df = pd.DataFrame(meta_data)
        cache.set("taxonomy_options", sorted(df["genus"].dropna().unique()))
        cache.set("family_options", sorted(df["familyName"].dropna().unique()))
        cache.set("navis_options", sorted(df["starship_navis"].dropna().unique()))
        cache.set("haplotype_options", sorted(df["starship_haplotype"].dropna().unique()))
        
        logger.info("Cache initialized successfully")
    except Exception as e:
        logger.error(f"Error initializing cache: {str(e)}", exc_info=True)
        raise


#########################
# Precompute cache
#########################

families = [
    "Phoenix",
    "Hephaestus",
    "Tardis",
    "Serenity",
    "Prometheus",
    "Enterprise",
    "Galactica",
    "Moya",
    "Arwing",
    "Voyager",
    "Family-11",
]


def precompute_all():
    """Precompute and cache all necessary data and figures."""
    cache_status = {}

    precompute_tasks = {
        "meta_data": lambda: fetch_meta_data(curated=True),
        "paper_data": fetch_paper_data,
        "download_data_curated_true_derep_false": lambda: fetch_download_data(curated=True, dereplicate=False),
        "download_data_curated_true_derep_true": lambda: fetch_download_data(curated=True, dereplicate=True),
        "download_data_curated_false_derep_false": lambda: fetch_download_data(curated=False, dereplicate=False),
        "download_data_curated_false_derep_true": lambda: fetch_download_data(curated=False, dereplicate=True),
        "all_ships": fetch_all_ships,
        "ship_table": fetch_ship_table,
        "all_captains": fetch_all_captains,
        "captain_tree": fetch_captain_tree,
        "sf_data": fetch_sf_data,
    }
    
    try:
        # Handle meta_data first as other operations depend on it
        meta_data = precompute_tasks["meta_data"]()
        df = pd.DataFrame(meta_data)
        
        # Cache both in memory and file
        cache.set("meta_data", meta_data)
        save_to_cache(meta_data, generate_cache_key("meta_data"))
        
        # Cache search options
        search_options = {
            "taxonomy_options": sorted(df["genus"].dropna().unique()),
            "family_options": sorted(df["familyName"].dropna().unique()),
            "navis_options": sorted(df["starship_navis"].dropna().unique()),
            "haplotype_options": sorted(df["starship_haplotype"].dropna().unique())
        }
        for key, value in search_options.items():
            cache.set(key, value)
            save_to_cache(value, generate_cache_key(key))
        
        cache_status["meta_data"] = True
        
        # Create sunburst plots (depends on meta_data)
        logger.info("Creating sunburst figures...")
        for family in families:
            try:
                sunburst_figure = cache_sunburst_plot(family, meta_data)
                cache_status[f"sunburst_{family}"] = True
            except Exception as e:
                logger.error(f"Failed to create sunburst for {family}: {str(e)}")
                cache_status[f"sunburst_{family}"] = False
        
        # Execute remaining precompute tasks
        for key, func in precompute_tasks.items():
            if key != "meta_data":  # Skip meta_data as it's already done
                try:
                    logger.info(f"Precomputing {key}...")
                    result = func()
                    if result is not None:
                        # Save both to memory cache and file cache
                        cache.set(key, result)
                        save_to_cache(result, generate_cache_key(key))
                        cache_status[key] = True
                        logger.info(f"Successfully cached {key}")
                    else:
                        cache_status[key] = False
                        logger.error(f"Failed to precompute {key}: returned None")
                except Exception as e:
                    logger.error(f"Failed to precompute {key}: {str(e)}")
                    cache_status[key] = False
                    
    except Exception as e:
        logger.error(f"Critical error during precomputation: {str(e)}")
        cache_status["meta_data"] = False
    
    # Log final status
    success_rate = sum(1 for v in cache_status.values() if v) / len(cache_status)
    logger.info(f"Precomputation complete. Success rate: {success_rate:.1%}")
    logger.info(f"Cache status: {cache_status}")
    
    return cache_status

def refresh_cache():
    """Refresh all cached data. Can be called after database updates."""
    try:
        logger.info("Starting cache refresh...")
        
        # Invalidate all existing cache
        for family in families:
            invalidate_cache(generate_cache_key("sunburst", family))
        
        # Invalidate all precomputed data
        for key in [
            "meta_data",
            "paper_data",
            "download_data_curated_true_derep_false",
            "download_data_curated_true_derep_true",
            "download_data_curated_false_derep_false",
            "download_data_curated_false_derep_true",
            "all_ships",
            "ship_table",
            "all_captains",
            "captain_tree",
            "sf_data",
            "taxonomy_options",
            "family_options",
            "navis_options",
            "haplotype_options"
        ]:
            invalidate_cache(generate_cache_key(key))
            cache.delete(key)  # Also clear from memory cache
        
        # Recompute everything
        cache_status = precompute_all()
        
        logger.info(f"Cache refresh complete. Status: {cache_status}")
        return {"success": True, "status": cache_status}
    except Exception as e:
        logger.error(f"Error refreshing cache: {str(e)}")
        return {"success": False, "error": str(e)}
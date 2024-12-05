import os
import pickle
from typing import Any
import logging
import pandas as pd
from src.components.cache import cache
from src.components.sql_engine import starbase_session_factory
from src.utils.plot_utils import create_sunburst_plot

logger = logging.getLogger(__name__)

@cache.memoize(timeout=300)
def fetch_meta_data(curated=False):
    """Fetch metadata from the database with caching."""
    session = starbase_session_factory()
    
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

@cache.memoize(timeout=300)
def fetch_paper_data():
    """Fetch paper data from the database and cache the result."""
    session = starbase_session_factory()

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

@cache.memoize(timeout=300)
def cache_sunburst_plot(family, df):
    """Create sunburst plots for wiki page and cache the result."""
    return create_sunburst_plot(df=df, type="tax", title_switch=False)


@cache.memoize(timeout=300)
def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    session = starbase_session_factory()

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

@cache.memoize(timeout=300)
def fetch_all_ships(curated=True):
    session = starbase_session_factory()

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
        return df
    except Exception as e:
        logger.error(f"Error fetching all_ships data: {str(e)}")
        raise
    finally:
        session.close()
    
@cache.memoize(timeout=300)
def fetch_accession_gff(accession):
    session = starbase_session_factory()

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

    try:
        df = pd.read_sql_query(query, session.bind, params=(accession,))

        if df.empty:
            logger.warning("Fetched gff DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching gff data: {str(e)}")
        raise
    finally:
        session.close()

@cache.memoize(timeout=300)
def fetch_ship_table(meta_df=None):
    """Fetch and filter ship table data based on metadata DataFrame."""
    session = starbase_session_factory()
    
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

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched ship_table DataFrame is empty.")
        df = df.sort_values(by="familyName", ascending=True)
        return df
    except Exception as e:
        logger.error(f"Error fetching ship_table data: {str(e)}")
        raise
    finally:
        session.close()


@cache.memoize(timeout=300)
def fetch_all_captains():
    session = starbase_session_factory()

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

@cache.memoize(timeout=300)
def fetch_captain_tree():
    session = starbase_session_factory()

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

@cache.memoize(timeout=300)
def fetch_sf_data():
    session = starbase_session_factory()

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

@cache.memoize(timeout=300)
def get_database_stats():
    """Get statistics about the Starship database."""
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
        return stats
    except Exception as e:
        logger.error(f"Error fetching database stats: {str(e)}")
        raise
    finally:
        session.close()

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
        "stats": get_database_stats,
    }
    
    try:
        # Handle meta_data first as other operations depend on it
        meta_data = precompute_tasks["meta_data"]()
        df = pd.DataFrame(meta_data)
        
        # Cache both in memory and file
        cache.set("meta_data", meta_data)
        
        # Cache search options
        search_options = {
            "taxonomy_options": sorted(df["genus"].dropna().unique()),
            "family_options": sorted(df["familyName"].dropna().unique()),
            "navis_options": sorted(df["starship_navis"].dropna().unique()),
            "haplotype_options": sorted(df["starship_haplotype"].dropna().unique())
        }
        for key, value in search_options.items():
            cache.set(key, value)
        
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
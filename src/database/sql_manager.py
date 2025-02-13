import os
import pickle
from typing import Any
import logging
import pandas as pd
from flask_caching import Cache
from src.config.cache import cache
from src.config.database import StarbaseSession
from src.utils.plot_utils import create_sunburst_plot

logger = logging.getLogger(__name__)

@cache.memoize()
def fetch_meta_data(curated=False):
    """Fetch metadata from the database with caching."""
    session = StarbaseSession()
    
    meta_query = """
    SELECT j.ship_family_id, j.curated_status, t.taxID, j.starshipID,
           j.ome, j.size, j.upDR, j.downDR, f.familyName, f.type_element_reference, j.contigID, 
           j.elementBegin, j.elementEnd, t.`order`, t.family, t.genus, t.species, 
           g.version, g.genomeSource, g.citation, a.accession_tag, g.strain, j.starship_navis, j.starship_haplotype, g.assembly_accession
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
def cache_sunburst_plot(family, df):
    """Create sunburst plots for wiki page and cache the result."""
    return create_sunburst_plot(df=df, type="tax", title_switch=False)


@cache.memoize()
def fetch_download_data(curated=True, dereplicate=False):
    """Fetch download data from the database and cache the result."""
    session = StarbaseSession()

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

@cache.memoize()
def fetch_all_ships(curated=True):
    session = StarbaseSession()

    query = """
    SELECT s.*, a.accession_tag
    FROM ships s
    JOIN accessions a ON s.accession = a.id
    JOIN joined_ships j ON j.ship_id = a.id
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
    
@cache.memoize()
def fetch_ship_table():
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    session = StarbaseSession()
    
    query = """
    SELECT DISTINCT 
        a.accession_tag,
        f.familyName,
        t.species
    FROM joined_ships js
    JOIN accessions a ON js.ship_id = a.id
    JOIN taxonomy t ON js.taxid = t.id
    JOIN family_names f ON js.ship_family_id = f.id
    -- Filter for ships that have sequence data
    JOIN ships s ON s.accession = a.id AND s.sequence IS NOT NULL
    -- Filter for ships that have GFF data
    JOIN gff g ON g.ship_id = a.id
    WHERE js.orphan IS NULL
    ORDER BY f.familyName ASC
    """

    try:
        df = pd.read_sql_query(query, session.bind)
        if df.empty:
            logger.warning("Fetched ship_table DataFrame is empty.")
        return df
    except Exception as e:
        logger.error(f"Error fetching ship_table data: {str(e)}")
        raise
    finally:
        session.close()

@cache.memoize()
def fetch_accession_ship(accession_tag):
    """Fetch sequence and GFF data for a specific ship."""
    session = StarbaseSession()
    
    sequence_query = """
    SELECT s.sequence
    FROM ships s
    JOIN accessions a ON s.accession = a.id
    WHERE a.accession_tag = :accession_tag
    """
    
    gff_query = """
    SELECT g.*
    FROM gff g
    JOIN accessions a ON g.ship_id = a.id
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

@cache.memoize()
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

@cache.memoize()
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

@cache.memoize()
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

@cache.memoize()
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
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

def precompute_all():
    """Precompute and cache all necessary data and figures."""
    from src.config.cache import cache
    from src.config.precompute import precompute_tasks
    
    cache_status = {}
    
    try:
        for key, func in precompute_tasks.items():
            try:
                logger.info(f"Precomputing {key}...")
                result = func()
                if result is not None:
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
    
    # Log final status
    success_rate = sum(1 for v in cache_status.values() if v) / len(cache_status)
    logger.info(f"Precomputation complete. Success rate: {success_rate:.1%}")
    logger.info(f"Cache status: {cache_status}")
    
    return cache_status
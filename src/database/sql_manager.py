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
from src.config.cache import smart_cache
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

# TODO: caching can be much more efficient, if we only cache the full dataset, and then apply curation/dereplication filters to the cached dataset afterwards
@db_retry_decorator()
@smart_cache(timeout=3600)
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
    SELECT j.curated_status, j.starshipID,
           a.accession_tag, a.version_tag,
           j.ship_id, j.id as joined_ship_id,
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
    INNER JOIN accessions a ON j.accession_id = a.id
    LEFT JOIN taxonomy t ON j.tax_id = t.id
    LEFT JOIN starship_features sf ON a.id = sf.accession_id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN navis_names n ON j.ship_navis_id = n.id
    LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    """

    # Build WHERE clause conditions
    where_conditions = []
    params = []
    
    if curated:
        where_conditions.append("j.curated_status = 'curated'")

    if accession_tags:
        if isinstance(accession_tags, list):
            placeholders = ",".join(["?"] * len(accession_tags))
            where_conditions.append(f"a.accession_tag IN ({placeholders})")
            params = list(accession_tags)
        else:
            where_conditions.append("a.accession_tag = ?")
            params = [accession_tags]
    
    if where_conditions:
        meta_query += f"\n    WHERE {' AND '.join(where_conditions)}"
    
    params = tuple(params) if params else []

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
@smart_cache(timeout=None)
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
@smart_cache(timeout=7200)
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
    INNER JOIN accessions a ON j.accession_id = a.id
    LEFT JOIN family_names f ON j.ship_family_id = f.id
    LEFT JOIN genomes g ON j.genome_id = g.id
    LEFT JOIN papers p ON f.type_element_reference = p.shortCitation
    -- Only show entries that have sequences
    INNER JOIN ships s ON s.id = j.ship_id
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
@smart_cache(timeout=3600)
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
            j.ship_id,
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
        INNER JOIN accessions a ON j.accession_id = a.id
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
            v.ship_id,
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
        LEFT JOIN ships s ON s.id = v.ship_id
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
def fetch_ship_table(curated=True, with_sequence=False, with_gff_entries=False):
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
    LEFT JOIN accessions a ON js.accession_id = a.id
    LEFT JOIN taxonomy t ON js.tax_id = t.id
    LEFT JOIN family_names f ON js.ship_family_id = f.id
    WHERE 1=1
    """

    if with_sequence:
        query += " AND js.ship_id IS NOT NULL"


    if with_gff_entries:
        with_gff_entries_query = """
        SELECT DISTINCT js.ship_id
        FROM joined_ships js
        LEFT JOIN gff g ON g.ship_id = js.ship_id
        WHERE g.source IS NOT NULL
        """
        query += f" AND js.ship_id IN ({with_gff_entries_query})"
    
    if curated:
        query += " AND js.curated_status = 'curated'"

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
    FROM joined_ships j
    LEFT JOIN ships s ON s.id = j.ship_id
    LEFT JOIN accessions a ON a.id = j.accession_id
    WHERE a.accession_tag = :accession_tag AND s.sequence IS NOT NULL
    """

    gff_query = """
    SELECT g.source, g.type, g.start, g.end, g.phase, g.strand, g.score, g.attributes
    FROM joined_ships j
    LEFT JOIN gff g ON g.ship_id = j.ship_id
    LEFT JOIN accessions a ON a.id = j.accession_id
    WHERE a.accession_tag = :accession_tag AND g.source IS NOT NULL
    """

    try:
        sequence_df = pd.read_sql_query(
            sequence_query, session.bind, params={"accession_tag": accession_tag}
        )
        if sequence_df.empty:
            logger.warning(f"No sequence data found for accession: {accession_tag}")
            sequence_df = None
        gff_df = pd.read_sql_query(
            gff_query, session.bind, params={"accession_tag": accession_tag}
        )
        if gff_df.empty:
            logger.warning(f"No GFF data found for accession: {accession_tag}")
            gff_df = None

        return {"sequence": sequence_df, "gff": gff_df}
    except Exception as e:
        logger.error(f"Error fetching sequence data for {accession_tag}: {str(e)}")
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
        INNER JOIN accessions a ON j.accession_id = a.id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
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
@smart_cache(timeout=0)
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
    try:
        # use metadata from previous query
        meta_df = fetch_meta_data(curated=False)
        filtered_df = meta_df[meta_df["accession_tag"].notna()]
        total_count = len(filtered_df["accession_tag"].unique())
        curated_count = len(filtered_df[filtered_df["curated_status"] == "curated"]["accession_tag"].unique())
        uncurated_count = total_count - curated_count
        species_count = len(filtered_df["name"].unique())
        families = filtered_df["familyName"].dropna().loc[~filtered_df["familyName"].isin(["NA", "None", None, "NULL"])].unique()
        logger.info(f"families: {families}")
        family_count = len(families)

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


def add_quality_tag(joined_ship_id, tag_type, tag_value=None, created_by="auto"):
    """
    Add a quality tag to a ship.
    
    Args:
        joined_ship_id (int): ID of the joined_ships record
        tag_type (str): Type of tag (e.g., "incomplete", "fragmented", "verified")
        tag_value (str, optional): Optional value for the tag
        created_by (str): Who created this tag (default: "auto")
    
    Returns:
        int: ID of the created tag, or existing tag ID if duplicate
    
    Example:
        >>> add_quality_tag(123, "incomplete", created_by="curator_name")
        >>> add_quality_tag(123, "nested", "inside_SS-1.1", created_by="auto")
    """
    from src.database.models.schema import ShipQualityTags
    from datetime import datetime
    
    session = StarbaseSession()
    try:
        # Check if tag already exists (unique constraint on joined_ship_id + tag_type)
        existing_tag = session.query(ShipQualityTags).filter_by(
            joined_ship_id=joined_ship_id,
            tag_type=tag_type
        ).first()
        
        if existing_tag:
            # Update the tag value if provided
            if tag_value is not None:
                existing_tag.tag_value = tag_value
                existing_tag.created_by = created_by
                session.commit()
            logger.info(f"Updated existing tag {tag_type} for ship {joined_ship_id}")
            return existing_tag.id
        
        # Create new tag
        new_tag = ShipQualityTags(
            joined_ship_id=joined_ship_id,
            tag_type=tag_type,
            tag_value=tag_value,
            created_at=datetime.now(),
            created_by=created_by
        )
        session.add(new_tag)
        session.commit()
        logger.info(f"Added tag {tag_type} to ship {joined_ship_id}")
        return new_tag.id
    except Exception as e:
        session.rollback()
        logger.error(f"Error adding quality tag: {str(e)}")
        raise
    finally:
        session.close()


def remove_quality_tag(joined_ship_id, tag_type):
    """
    Remove a quality tag from a ship.
    
    Args:
        joined_ship_id (int): ID of the joined_ships record
        tag_type (str): Type of tag to remove
    
    Returns:
        bool: True if tag was removed, False if not found
    """
    from src.database.models.schema import ShipQualityTags
    
    session = StarbaseSession()
    try:
        tag = session.query(ShipQualityTags).filter_by(
            joined_ship_id=joined_ship_id,
            tag_type=tag_type
        ).first()
        
        if tag:
            session.delete(tag)
            session.commit()
            logger.info(f"Removed tag {tag_type} from ship {joined_ship_id}")
            return True
        else:
            logger.warning(f"Tag {tag_type} not found for ship {joined_ship_id}")
            return False
    except Exception as e:
        session.rollback()
        logger.error(f"Error removing quality tag: {str(e)}")
        raise
    finally:
        session.close()


def get_quality_tags(joined_ship_id):
    """
    Get all quality tags for a ship.
    
    Args:
        joined_ship_id (int): ID of the joined_ships record
    
    Returns:
        list: List of dicts with tag information
    """
    from src.database.models.schema import ShipQualityTags
    
    session = StarbaseSession()
    try:
        tags = session.query(ShipQualityTags).filter_by(
            joined_ship_id=joined_ship_id
        ).all()
        
        return [
            {
                "tag_type": tag.tag_type,
                "tag_value": tag.tag_value,
                "created_at": tag.created_at,
                "created_by": tag.created_by
            }
            for tag in tags
        ]
    except Exception as e:
        logger.error(f"Error fetching quality tags: {str(e)}")
        raise
    finally:
        session.close()

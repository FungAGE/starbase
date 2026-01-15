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


@db_retry_decorator()
def fetch_meta_data(curated=False, accession_tags=None):
    """
    Fetch metadata from the database with efficient caching.

    Caches the full dataset and applies filters in memory for maximum efficiency.

    Args:
        curated (bool): If True, only return curated entries
        accession_tags (str or list): Single accession tag or list of accession tags

    Returns:
        pd.DataFrame: Metadata for the specified filters
    """
    from src.config.cache import cache

    # Always cache the full dataset with a fixed key
    cache_key = "fetch_meta_data:full_dataset"

    # Try to get cached full dataset
    full_df = cache.get(cache_key)
    if full_df is not None:
        if isinstance(full_df, dict) and "pandas_df" in full_df:
            full_df = pd.DataFrame.from_dict(full_df["pandas_df"])
    else:
        # Fetch and cache the full dataset
        session = StarbaseSession()
        try:
            meta_query = """
            SELECT j.curated_status, j.starshipID,
                   a.accession_tag, a.version_tag,
                   j.ship_id, j.id as joined_ship_id,
                   CASE
                       WHEN a.version_tag IS NOT NULL AND a.version_tag != ''
                       THEN a.accession_tag || '.' || a.version_tag
                       ELSE a.accession_tag
                   END as accession_display,
                   sa.ship_accession_tag,
                   CASE
                       WHEN sa.version_tag IS NOT NULL AND sa.version_tag != ''
                       THEN sa.ship_accession_tag || '.' || sa.version_tag
                       ELSE sa.ship_accession_tag
                   END as ship_accession_display,
                   t.taxID, t.strain, t.`order`, t.family, t.name,
                   sf.elementLength, sf.upDR, sf.downDR, sf.contigID, sf.captainID, sf.elementBegin, sf.elementEnd,
                   f.familyName, f.type_element_reference, n.navis_name, h.haplotype_name,
                   g.ome, g.version, g.genomeSource, g.citation, g.assembly_accession,
                   s.md5, s.rev_comp_md5, s.sequence_length
            FROM joined_ships j
            LEFT JOIN accessions a ON j.accession_id = a.id
            LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
            LEFT JOIN taxonomy t ON j.tax_id = t.id
            LEFT JOIN starship_features sf ON a.id = sf.accession_id
            LEFT JOIN family_names f ON j.ship_family_id = f.id
            LEFT JOIN navis_names n ON j.ship_navis_id = n.id
            LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
            LEFT JOIN genomes g ON j.genome_id = g.id
            LEFT JOIN ships s ON s.id = j.ship_id
            WHERE j.accession_id IS NOT NULL
            """

            full_df = pd.read_sql_query(meta_query, session.bind)
            # Cache as dictionary for DataFrame serialization
            cache.set(cache_key, {"pandas_df": full_df.to_dict()}, timeout=None)

        except Exception as e:
            logger.error(f"Error fetching meta data: {str(e)}")
            raise
        finally:
            session.close()

    # Apply filters in memory
    filtered_df = full_df.copy()

    if curated:
        filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]

    if accession_tags:
        if isinstance(accession_tags, list):
            filtered_df = filtered_df[filtered_df["accession_tag"].isin(accession_tags)]
        else:
            filtered_df = filtered_df[filtered_df["accession_tag"] == accession_tags]

    return filtered_df


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


def dereplicate_sequences(df):
    """Deduplicate sequences based on MD5 hashes (forward and reverse complement)."""
    seen_sequences = set()
    indices_to_keep = []

    for idx, row in df.iterrows():
        md5_val = row.get("md5", "")
        rev_comp_md5_val = row.get("rev_comp_md5", "")

        if not md5_val or not rev_comp_md5_val:
            indices_to_keep.append(idx)
            continue
        if md5_val not in seen_sequences and rev_comp_md5_val not in seen_sequences:
            indices_to_keep.append(idx)
            seen_sequences.add(md5_val)
            seen_sequences.add(rev_comp_md5_val)

    filtered_df = df.loc[indices_to_keep]
    return filtered_df


# TODO: figure out a way to handle caching with queries related to this query
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
            j.ship_id,
            CASE 
                WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
                THEN a.accession_tag || '.' || a.version_tag
                ELSE a.accession_tag
            END as accession_display,
            sa.ship_accession_tag,
            CASE
                WHEN sa.version_tag IS NOT NULL AND sa.version_tag != ''
                THEN sa.ship_accession_tag || '.' || sa.version_tag
                ELSE sa.ship_accession_tag
            END as ship_accession_display,
            j.curated_status,
            j.starshipID,
            sf.elementBegin, sf.elementEnd, sf.contigID,
            t.name, t.family, t.`order`,
            f.familyName, n.navis_name, h.haplotype_name,
            g.assembly_accession, c.captainID"""

    base_query += """
        FROM joined_ships j
        INNER JOIN accessions a ON j.accession_id = a.id
        LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
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
            v.ship_accession_tag,
            v.ship_accession_display,
            v.curated_status,
            v.starshipID,
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
            v.ship_accession_tag,
            v.ship_accession_display,
            v.curated_status,
            v.starshipID,
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

        query += """
        """

    try:
        df = pd.read_sql_query(query, session.bind)

        if df.empty:
            logger.warning("Fetched ships DataFrame is empty.")
            return df

        # Apply MD5-based deduplication if sequences are available and deduplication is requested
        if (
            with_sequence
            and dereplicate
            and "md5" in df.columns
            and "rev_comp_md5" in df.columns
        ):
            df = dereplicate_sequences(df)

        return df
    except Exception as e:
        logger.error(f"Error fetching ships data: {str(e)}")
        raise
    finally:
        session.close()


@smart_cache(timeout=None)
@db_retry_decorator()
def fetch_ship_table(curated=True, with_sequence=False, with_gff_entries=False):
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    from src.config.cache import cache

    cache_key = "fetch_ship_table:full_dataset"

    full_df = cache.get(cache_key)
    if full_df is not None:
        if isinstance(full_df, dict) and "pandas_df" in full_df:
            full_df = pd.DataFrame.from_dict(full_df["pandas_df"])
    else:
        session = StarbaseSession()

        try:
            query = """
            SELECT DISTINCT
                js.ship_id,
                js.source,
                js.curated_status,
                a.accession_tag, a.version_tag,
                CASE
                    WHEN a.version_tag IS NOT NULL AND a.version_tag != ''
                    THEN a.accession_tag || '.' || a.version_tag
                    ELSE a.accession_tag
                END as accession_display,
                sa.ship_accession_tag,
                CASE
                    WHEN sa.version_tag IS NOT NULL AND sa.version_tag != ''
                    THEN sa.ship_accession_tag || '.' || sa.version_tag
                    ELSE sa.ship_accession_tag
                END as ship_accession_display,
                f.familyName,
                t.name
            FROM joined_ships js
            LEFT JOIN accessions a ON js.accession_id = a.id
            LEFT JOIN ship_accessions sa ON sa.ship_id = js.ship_id
            LEFT JOIN taxonomy t ON js.tax_id = t.id
            LEFT JOIN family_names f ON js.ship_family_id = f.id
            WHERE 1=1
            """

            full_df = pd.read_sql_query(query, session.bind)
            cache.set(cache_key, {"pandas_df": full_df.to_dict()}, timeout=None)

        except Exception as e:
            logger.error(f"Error fetching ship table data: {str(e)}")
            raise
        finally:
            session.close()

    filtered_df = full_df.copy()

    if with_sequence:
        filtered_df = filtered_df[filtered_df["ship_id"].notna()]

    if with_gff_entries:
        # generate a list of ship_ids that have GFF entries, using a separate query
        session = StarbaseSession()
        try:
            gff_query = """
            SELECT DISTINCT g.ship_id
            FROM gff g
            WHERE g.ship_id IS NOT NULL AND g.ship_id != ''
            """
            gff_df = pd.read_sql_query(gff_query, session.bind)
            gff_ship_ids = gff_df["ship_id"].dropna().tolist()
            filtered_df = filtered_df[filtered_df["ship_id"].isin(gff_ship_ids)]
        except Exception as e:
            logger.error(f"Error fetching GFF data for filtering: {str(e)}")
            # If GFF query fails, return empty dataframe to ensure no entries without GFF data are shown
            filtered_df = filtered_df.iloc[0:0]
        finally:
            session.close()

    if curated:
        filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]

    filtered_df = filtered_df.sort_values(by="familyName")

    return filtered_df


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
            sa.ship_accession_tag,
            CASE
                WHEN sa.version_tag IS NOT NULL AND sa.version_tag != ''
                THEN sa.ship_accession_tag || '.' || sa.version_tag
                ELSE sa.ship_accession_tag
            END as ship_accession_display,
            j.curated_status,
            j.starshipID,
            c.captainID as captain_id,
            c."sequence",
            n.navis_name,
            h.haplotype_name,
            c.captainID
        FROM joined_ships j
        INNER JOIN accessions a ON j.accession_id = a.id
        LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
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
            v.ship_accession_tag,
            v.ship_accession_display,
            v.curated_status,
            v.starshipID,
            v.captain_id,
            v.sequence,
            v.navis_name,
            v.haplotype_name,
            v.captain_id as captain_id_col
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
            v.ship_accession_tag,
            v.ship_accession_display,
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
def get_database_version():
    """Get the current database semantic version from the database_versions table."""
    session = StarbaseSession()
    try:
        result = session.execute(
            text("""
            SELECT semantic_version FROM database_versions
            ORDER BY created_at DESC LIMIT 1
        """)
        ).fetchone()

        return result[0] if result else "unknown"
    except Exception as e:
        logger.error(f"Error fetching database version: {str(e)}")
        return "unknown"
    finally:
        session.close()


@db_retry_decorator()
def set_database_version(semantic_version, description="", created_by="manual"):
    """Manually set a new semantic version for the database."""
    session = StarbaseSession()
    try:
        session.execute(
            text("""
            INSERT INTO database_versions (semantic_version, description, created_by)
            VALUES (:version, :desc, :creator)
        """),
            {"version": semantic_version, "desc": description, "creator": created_by},
        )
        session.commit()
        logger.info(f"Database version manually set to {semantic_version}")
        return True
    except Exception as e:
        session.rollback()
        logger.error(f"Error setting database version: {str(e)}")
        raise
    finally:
        session.close()


@db_retry_decorator()
def get_alembic_schema_version():
    """Get the current Alembic schema version (for schema tracking)."""
    try:
        from alembic.migration import MigrationContext

        with StarbaseSession() as session:
            conn = session.connection()
            context = MigrationContext.configure(conn)
            current_rev = context.get_current_revision()

        return current_rev if current_rev else "unknown"
    except Exception as e:
        logger.error(f"Error fetching Alembic schema version: {str(e)}")
        return "unknown"


@db_retry_decorator()
@smart_cache(timeout=None)
def get_database_stats():
    """Get statistics about the Starship database."""
    session = StarbaseSession()
    try:
        # new sql query for stats
        # stats metdata
        stats_metadata_query = """
        SELECT j.curated_status, j.starshipID,
                a.accession_tag, a.version_tag,
                j.ship_id, j.id as joined_ship_id,
                CASE
                    WHEN a.version_tag IS NOT NULL AND a.version_tag != ''
                    THEN a.accession_tag || '.' || a.version_tag
                    ELSE a.accession_tag
                END as accession_display,
                sa.ship_accession_tag,
                CASE
                    WHEN sa.version_tag IS NOT NULL AND sa.version_tag != ''
                    THEN sa.ship_accession_tag || '.' || sa.version_tag
                    ELSE sa.ship_accession_tag
                END as ship_accession_display,
                t.taxID, t.strain, t.`order`, t.family, t.name,
                sf.elementLength, sf.upDR, sf.downDR, sf.contigID, sf.captainID, sf.elementBegin, sf.elementEnd,
                f.familyName, f.type_element_reference, n.navis_name, h.haplotype_name,
                g.ome, g.version, g.genomeSource, g.citation, g.assembly_accession, s.md5, s.rev_comp_md5
        FROM joined_ships j
        LEFT JOIN accessions a ON j.accession_id = a.id
        LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN starship_features sf ON a.id = sf.accession_id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
        LEFT JOIN genomes g ON j.genome_id = g.id
        LEFT JOIN ships s ON s.id = j.ship_id
        """
        stats_df = pd.read_sql_query(stats_metadata_query, session.bind)

        # total numer of ships (regardless of duplicates or sequencing similarity)
        total_count = len(stats_df)
        # total number of unique sequences (by md5 or rev_comp_md5)
        unique_sequences_df = stats_df[["md5", "rev_comp_md5"]].drop_duplicates()
        unique_sequences_count = len(unique_sequences_df)

        curated_count = len(stats_df[stats_df["curated_status"] == "curated"])
        uncurated_count = total_count - curated_count
        species_count = len(stats_df["name"].unique())
        family_count = len(
            stats_df["familyName"]
            .dropna()
            .loc[~stats_df["familyName"].isin(["NA", "None", None, "NULL"])]
            .unique()
        )

        stats = {
            "total_starships": total_count,
            "unique_sequences": unique_sequences_count,
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
        existing_tag = (
            session.query(ShipQualityTags)
            .filter_by(joined_ship_id=joined_ship_id, tag_type=tag_type)
            .first()
        )

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
            created_by=created_by,
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
        tag = (
            session.query(ShipQualityTags)
            .filter_by(joined_ship_id=joined_ship_id, tag_type=tag_type)
            .first()
        )

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
        tags = (
            session.query(ShipQualityTags)
            .filter_by(joined_ship_id=joined_ship_id)
            .all()
        )

        return [
            {
                "tag_type": tag.tag_type,
                "tag_value": tag.tag_value,
                "created_at": tag.created_at,
                "created_by": tag.created_by,
            }
            for tag in tags
        ]
    except Exception as e:
        logger.error(f"Error fetching quality tags: {str(e)}")
        raise
    finally:
        session.close()

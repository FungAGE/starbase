import pandas as pd
from src.database.sql_engine import get_starbase_session
from sqlalchemy import text
from src.config.cache import smart_cache
from src.config.settings import PHYLOGENY_PATHS

from src.config.logging import get_logger

logger = get_logger(__name__)


def _get_accession_mode(accessions):
    # the set of accessions should all start with either "SSA" or all start with "SSB"
    if accessions:
        if all(accession.startswith("SSA") for accession in accessions):
            return "SSA"
        if all(accession.startswith("SSB") for accession in accessions):
            return "SSB"
    return None


def fetch_meta_data(curated=False, accessions=None):
    if isinstance(accessions, str):
        accessions = [accessions]
    """
    Fetch metadata from the database with efficient caching.

    Caches the full dataset and applies filters in memory for maximum efficiency.

    Args:
        curated (bool): If True, only return curated entries
        accessions (str or list): Single accession tag or list of accession tags
        
    Returns:
        pd.DataFrame: Metadata for the specified filters
    """
    from src.config.cache import cache

    # Always cache the full dataset with a fixed key (bump version when schema changes)
    cache_key = "fetch_meta_data:full_dataset:v2"

    accession_mode = _get_accession_mode(accessions)

    with get_starbase_session() as session:
        try:
            # Try to get cached full dataset
            full_df = cache.get(cache_key)
            if full_df is not None:
                if isinstance(full_df, dict) and "pandas_df" in full_df:
                    full_df = pd.DataFrame.from_dict(full_df["pandas_df"])
            else:
                meta_query = """
                SELECT curated_status, starshipID,
                    ship_id, joined_ship_id,
                    ship_accession_tag,
                    ship_version_tag,
                    ship_accession_display,
                    taxID, strain, "order", family, name,
                    elementLength, upDR, downDR, contigID, captainID, elementBegin, elementEnd,
                    familyName, type_element_reference, navis_name, navis_activity, haplotype_name, haplotype_activity,
                    ome, genome_version as version, genomeSource, citation, assembly_accession,
                    md5, rev_comp_md5, type_ship,
                    accession_tag, version_tag, accession_display
                FROM ships_with_metadata
                """

                full_df = pd.read_sql_query(meta_query, session.bind)
                # Cache as dictionary for DataFrame serialization
                cache.set(cache_key, {"pandas_df": full_df.to_dict()}, timeout=None)

        except Exception as e:
            logger.error(f"Error fetching meta data: {str(e)}")
            raise

        filtered_df = full_df.copy()

        if curated:
            filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]
        if accessions:
            formatted_values = [str(tag).strip("'\"") for tag in accessions]

            if accession_mode == "SSA":
                filtered_df = filtered_df[
                    filtered_df["accession_tag"].isin(formatted_values)
                ]
            elif accession_mode == "SSB":
                filtered_df = filtered_df[
                    filtered_df["ship_accession_tag"].isin(formatted_values)
                ]
            elif accession_mode is None:
                # Mixed or unknown: match either column
                mask = filtered_df["accession_tag"].isin(
                    formatted_values
                ) | filtered_df["ship_accession_tag"].isin(formatted_values)
                filtered_df = filtered_df[mask]
            else:
                raise ValueError(f"Invalid accession mode: {accession_mode}")

    return filtered_df


@smart_cache(timeout=None)
def fetch_paper_data():
    """Fetch paper data from the database and cache the result."""
    with get_starbase_session() as session:
        try:
            paper_query = """
            SELECT p.Title, p.Author, p.PublicationYear, p.DOI, p.Url, 
                   p.shortCitation, f.familyName, f.type_element_reference
            FROM papers p
            LEFT JOIN family_names f ON p.shortCitation = f.type_element_reference
            """
            paper_df = pd.read_sql_query(paper_query, session.bind)
            if paper_df.empty:
                logger.warning("Fetched paper DataFrame is empty.")
            return paper_df
        except Exception as e:
            logger.error(f"Error fetching paper data: {str(e)}")
            raise


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
def fetch_ships(
    accessions=None,
    curated=False,
    dereplicate=True,
    with_sequence=False,
):
    """
    Fetch ship data for specified accession tags.

    Args:
        accessions (list, optional): List of accession tags to fetch. If None, fetches all ships.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing ship data
    """
    accession_mode = _get_accession_mode(accessions)

    with get_starbase_session() as session:
        # ships_with_metadata is a real SQL view (see alembic/versions/i2j3k4l5m6n7_*)
        # mirroring the join previously duplicated here as an inline CTE. Column
        # list kept narrow (not SELECT *) so DISTINCT collapses the same rows it
        # used to -- the view carries extra columns (sequence, feature detail)
        # that would otherwise widen the distinct set.
        if with_sequence:
            select_cols = """
                sm.ship_id,
                sm.accession_tag,
                sm.version_tag,
                sm.accession_display,
                sm.ship_accession_tag,
                sm.ship_accession_display,
                sm.curated_status,
                sm.starshipID,
                sm.elementBegin,
                sm.elementEnd,
                sm.contigID,
                sm.name,
                sm.family,
                sm."order",
                sm.familyName,
                sm.navis_name,
                sm.haplotype_name,
                sm.assembly_accession,
                sm.type_ship,
                sm.sequence,
                sm.md5,
                sm.rev_comp_md5,
                sm.captainID"""
        else:
            select_cols = """
                sm.accession_tag,
                sm.version_tag,
                sm.accession_display,
                sm.ship_accession_tag,
                sm.ship_accession_display,
                sm.curated_status,
                sm.starshipID,
                sm.elementBegin,
                sm.elementEnd,
                sm.contigID,
                sm.name,
                sm.family,
                sm."order",
                sm.familyName,
                sm.navis_name,
                sm.haplotype_name,
                sm.assembly_accession,
                sm.type_ship,
                sm.captainID"""

        query = f"""
        SELECT DISTINCT
        {select_cols}
        FROM ships_with_metadata sm
        WHERE 1=1
        """

        if accessions:
            formatted_values = [str(tag).strip("'\"") for tag in accessions]
            quoted_sql = ", ".join(
                "'" + str(v).replace("'", "''") + "'" for v in formatted_values
            )
            # Match base accession with or without version suffix
            like_clauses = " OR ".join(
                "sm.accession_tag LIKE '" + str(v).replace("'", "''") + ".%'"
                for v in formatted_values
            )
            ship_like_clauses = " OR ".join(
                "sm.ship_accession_tag LIKE '" + str(v).replace("'", "''") + ".%'"
                for v in formatted_values
            )

            if accession_mode == "SSA":
                query += " AND (sm.accession_tag IN ({}) OR ({}))".format(
                    quoted_sql, like_clauses
                )
            elif accession_mode == "SSB":
                query += " AND (sm.ship_accession_tag IN ({}) OR ({}))".format(
                    quoted_sql, ship_like_clauses
                )
            elif accession_mode is None:
                # Mixed or unknown: match either column (e.g. download with mixed SSA/SSB or null SSA)
                query += " AND ((sm.accession_tag IN ({}) OR ({})) OR (sm.ship_accession_tag IN ({}) OR ({})))".format(
                    quoted_sql,
                    like_clauses,
                    quoted_sql,
                    ship_like_clauses,
                )
            else:
                raise ValueError(f"Invalid accession mode: {accession_mode}")
        if curated:
            query += " AND sm.curated_status = 'curated'"

        if with_sequence:
            query += " AND sm.sequence IS NOT NULL"

        try:
            df = pd.read_sql_query(query, session.bind)

            if (
                with_sequence
                and dereplicate
                and "md5" in df.columns
                and "rev_comp_md5" in df.columns
            ):
                # Apply MD5-based deduplication if sequences are available and deduplication is requested
                df = dereplicate_sequences(df)

        except Exception as e:
            logger.error(f"Error fetching ships data: {str(e)}")
            raise
        return df


@smart_cache(timeout=None)
def fetch_ship_table(curated=True, with_sequence=False, with_gff_entries=False):
    """Fetch ship metadata and filter for those with sequence and GFF data."""
    from src.config.cache import cache

    cache_key = "fetch_ship_table:full_dataset"

    full_df = cache.get(cache_key)
    with get_starbase_session() as session:
        if full_df is not None:
            if isinstance(full_df, dict) and "pandas_df" in full_df:
                full_df = pd.DataFrame.from_dict(full_df["pandas_df"])
        else:
            try:
                # ship_table_view mirrors this narrower join (see
                # alembic/versions/i2j3k4l5m6n7_*) — no features/genomes/
                # navis/haplotype/captains, to avoid row fan-out.
                query = "SELECT * FROM ship_table_view"

                full_df = pd.read_sql_query(query, session.bind)

                cache.set(cache_key, {"pandas_df": full_df.to_dict()}, timeout=None)

            except Exception as e:
                logger.error(f"Error fetching ship table data: {str(e)}")
                raise

        filtered_df = full_df.copy()

        if with_sequence:
            filtered_df = filtered_df[filtered_df["ship_id"].notna()]

        if with_gff_entries:
            # generate a list of ship_ids that have GFF entries, using a separate query
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

        if curated:
            filtered_df = filtered_df[filtered_df["curated_status"] == "curated"]

        filtered_df = filtered_df.sort_values(by="familyName")

    return filtered_df


def fetch_accession_ship(ship_accession_tag):
    """Fetch sequence and GFF data for a specific ship."""

    sequence_query = """
    SELECT s.sequence
    FROM joined_ships j
    LEFT JOIN ships s ON s.id = j.ship_id
    LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
    LEFT JOIN accessions a ON a.id = j.accession_id
    WHERE j.ship_id IS NOT NULL AND sa.ship_accession_tag = :ship_accession_tag AND s.sequence IS NOT NULL
    """

    gff_query = """
    SELECT g.source, g.type, g.start, g.end, g.phase, g.strand, g.score, g.attributes
    FROM joined_ships j
    LEFT JOIN gff g ON g.ship_id = j.ship_id
    LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
    LEFT JOIN accessions a ON a.id = j.accession_id
    WHERE j.ship_id IS NOT NULL AND sa.ship_accession_tag = :ship_accession_tag AND g.source IS NOT NULL
    """

    with get_starbase_session() as session:
        try:
            sequence_df = pd.read_sql_query(
                sequence_query,
                session.bind,
                params={"ship_accession_tag": ship_accession_tag},
            )
            if sequence_df.empty:
                logger.warning(
                    f"No sequence data found for accession: {ship_accession_tag}"
                )
                sequence_df = None
            gff_df = pd.read_sql_query(
                gff_query,
                session.bind,
                params={"ship_accession_tag": ship_accession_tag},
            )
            if gff_df.empty:
                logger.warning(f"No GFF data found for accession: {ship_accession_tag}")
                gff_df = None

            return {"sequence": sequence_df, "gff": gff_df}
        except Exception as e:
            logger.error(
                f"Error fetching sequence data for {ship_accession_tag}: {str(e)}"
            )
            raise


def fetch_captains(
    accessions=None, curated=False, dereplicate=True, with_sequence=False
):
    """
    Fetch captain data for specified accession tags.

    Args:
        accessions (list, optional): List of accession tags to fetch. If None, fetches all captains.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
            Rows are restricted to ``sequence IS NOT NULL`` (no captain sequence ⇒ not useful
            for mmseqs). When True, ``accession_display`` uses COALESCE over SSA/SSB tags,
            ``starshipID``, and ``captain_<id>`` because those join columns can still be NULL
            even when a captain sequence exists—those NULLs broke FASTA filenames / dict keys,
            not missing captains.

    Returns:
        pd.DataFrame: DataFrame containing captain data
    """

    # captains_with_metadata is a real SQL view (see
    # alembic/versions/i2j3k4l5m6n7_*) mirroring the join previously
    # duplicated here as an inline CTE.
    where_clauses = ["1=1"]
    if accessions:
        where_clauses.append(
            "cm.ship_accession_tag IN ({})".format(
                ",".join(f"'{tag}'" for tag in accessions)
            )
        )
    if curated:
        where_clauses.append("cm.curated_status = 'curated'")

    if with_sequence:
        where_clauses.append("cm.captain_sequence IS NOT NULL")
        query = """
        SELECT
            cm.ship_accession_tag,
            cm.version_tag,
            cm.ship_accession_display,
            cm.accession_tag,
            COALESCE(
                NULLIF(TRIM(COALESCE(CAST(cm.accession_display AS TEXT), '')), ''),
                NULLIF(TRIM(COALESCE(CAST(cm.accession_tag AS TEXT), '')), ''),
                NULLIF(TRIM(COALESCE(CAST(cm.ship_accession_display AS TEXT), '')), ''),
                NULLIF(TRIM(COALESCE(CAST(cm.ship_accession_tag AS TEXT), '')), ''),
                NULLIF(TRIM(COALESCE(CAST(cm.starshipID AS TEXT), '')), ''),
                'captain_' || CAST(cm.captain_row_id AS TEXT)
            ) AS accession_display,
            cm.curated_status,
            cm.starshipID,
            cm.captainID as captain_id,
            cm.captain_sequence as "sequence",
            cm.navis_name,
            cm.haplotype_name,
            cm.captainID as captain_id_col
        FROM captains_with_metadata cm
        WHERE """ + " AND ".join(where_clauses)
    else:
        query = """
        SELECT
            cm.accession_tag,
            cm.version_tag,
            cm.accession_display,
            cm.ship_accession_tag,
            cm.ship_accession_display,
            cm.curated_status,
            cm.starshipID,
            cm.captainID,
            cm.navis_name,
            cm.haplotype_name
        FROM captains_with_metadata cm
        WHERE """ + " AND ".join(where_clauses)

    with get_starbase_session() as session:
        try:
            df = pd.read_sql_query(query, session.bind)

            if dereplicate:
                df = df.drop_duplicates(subset="ship_accession_tag")

            if df.empty:
                logger.warning("Fetched captains DataFrame is empty.")
            return df
        except Exception as e:
            logger.error(f"Error fetching captains data: {str(e)}")
            raise


def fetch_captain_tree():
    fallback_tree_path = PHYLOGENY_PATHS["tree"]

    with open(fallback_tree_path, "r") as f:
        return f.read()


def fetch_sf_data():
    sf_data = pd.read_csv(PHYLOGENY_PATHS["clades"], sep="\t")

    # Add debug logging
    logger.debug(f"Loaded sf_data columns: {sf_data.columns.tolist()}")
    logger.debug(f"Loaded sf_data head: \n{sf_data.head()}")

    return sf_data


def get_database_version():
    """Get the current database semantic version from the database_versions table."""
    with get_starbase_session() as session:
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


def set_database_version(semantic_version, description="", created_by="manual"):
    """Manually set a new semantic version for the database."""
    with get_starbase_session() as session:
        try:
            session.execute(
                text("""
                INSERT INTO database_versions (semantic_version, description, created_by)
                VALUES (:version, :desc, :creator)
            """),
                {
                    "version": semantic_version,
                    "desc": description,
                    "creator": created_by,
                },
            )
            session.commit()
            logger.info(f"Database version manually set to {semantic_version}")
            return True
        except Exception as e:
            session.rollback()
            logger.error(f"Error setting database version: {str(e)}")
            raise


def get_alembic_schema_version():
    """
    Get the current Alembic schema version (for schema tracking).
    - Try to get single revision first
    - If multiple heads exist, get all heads
    - If no heads exist, return "unknown"
    """
    try:
        from alembic.migration import MigrationContext

        with get_starbase_session() as session:
            conn = session.connection()
            context = MigrationContext.configure(conn)

            try:
                current_rev = context.get_current_revision()
                return current_rev if current_rev else "unknown"
            except Exception:
                heads = context.get_current_heads()
                if heads:
                    return ", ".join(heads) if len(heads) > 1 else heads[0]
                return "unknown"
    except Exception as e:
        logger.error(f"Error fetching Alembic schema version: {str(e)}")
        return "unknown"


@smart_cache(timeout=None)
def get_database_stats():
    """Get statistics about the Starship database."""

    stats_metadata_query = """
    SELECT j.curated_status,
            j.ship_id,
            sa.ship_accession_tag,
            t.name,
            s.md5, s.rev_comp_md5
    FROM joined_ships j
    LEFT JOIN ship_accessions sa ON sa.ship_id = j.ship_id
    LEFT JOIN taxonomy t ON j.tax_id = t.id
    LEFT JOIN ships s ON s.id = j.ship_id
    WHERE j.ship_id IS NOT NULL
    """
    with get_starbase_session() as session:
        try:
            stats_df = pd.read_sql_query(stats_metadata_query, session.bind)

            # total numer of ships (regardless of duplicates or sequencing similarity)
            total_count = len(
                stats_df["ship_accession_tag"]
                .dropna()
                .loc[~stats_df["ship_accession_tag"].isin(["NA", "None", None, "NULL"])]
                .unique()
            )

            # total number of unique sequences (by md5 or rev_comp_md5)
            unique_sequences_df = stats_df[["md5", "rev_comp_md5"]].drop_duplicates()
            unique_sequences_count = len(unique_sequences_df)

            curated_count = len(stats_df[stats_df["curated_status"] == "curated"])
            uncurated_count = total_count - curated_count
            species_count = len(stats_df["name"].unique())

            # count from table "family_names"
            family_query = """
            SELECT DISTINCT familyName
            FROM family_names
            """
            family_df = pd.read_sql_query(family_query, session.bind)
            family_count = len(family_df["familyName"].unique())

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

    with get_starbase_session() as session:
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
                logger.info(
                    f"Updated existing tag {tag_type} for ship {joined_ship_id}"
                )
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

    with get_starbase_session() as session:
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


def set_ship_deleted(joined_ship_id, deleted=True):
    """
    Toggle the soft-delete flag on a joined_ships record.

    Soft-deleted ships are excluded from ships_with_metadata,
    captains_with_metadata, and ship_table_view automatically.

    Args:
        joined_ship_id (int): ID of the joined_ships record
        deleted (bool): True to soft-delete, False to restore

    Returns:
        bool: True if a row was updated, False if not found
    """
    from src.database.models.schema import JoinedShips

    with get_starbase_session() as session:
        try:
            joined_ship = (
                session.query(JoinedShips).filter_by(id=joined_ship_id).first()
            )
            if not joined_ship:
                logger.warning(f"joined_ships row {joined_ship_id} not found")
                return False

            joined_ship.is_deleted = deleted
            session.commit()
            logger.info(
                f"Set is_deleted={deleted} for joined_ships row {joined_ship_id}"
            )
            return True
        except Exception as e:
            session.rollback()
            logger.error(f"Error setting soft-delete flag: {str(e)}")
            raise


def get_quality_tags(joined_ship_id):
    """
    Get all quality tags for a ship.

    Args:
        joined_ship_id (int): ID of the joined_ships record

    Returns:
        list: List of dicts with tag information
    """
    from src.database.models.schema import ShipQualityTags

    with get_starbase_session() as session:
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

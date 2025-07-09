"""
Test helper utilities for running classification tests outside Flask app context.
Provides non-cached versions of database functions that can be used in standalone tests.
"""

import pandas as pd
import os
import sys
from contextlib import contextmanager

from src.config.database import StarbaseSession
from src.config.logging import get_logger
from tenacity import (
    retry,
    stop_after_attempt,
    wait_exponential,
    retry_if_exception_type,
)
import sqlalchemy.exc

# Add src to path if not already there
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(os.path.dirname(current_dir))
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

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


@contextmanager
def test_db_session():
    """Context manager for test database sessions"""
    session = StarbaseSession()
    try:
        yield session
    except Exception as e:
        logger.error(f"Database error: {str(e)}")
        session.rollback()
        raise
    finally:
        session.close()


def test_fetch_meta_data(curated=False, accession_tag=None):
    """
    Test-friendly version of fetch_meta_data without caching.

    Args:
        curated (bool): If True, only return curated entries
        accession_tag (str or list): Single accession tag or list of accession tags

    Returns:
        pd.DataFrame: Metadata for the specified accession tags
    """
    with test_db_session() as session:
        meta_query = """
        SELECT j.curated_status, j.starshipID, j.ship_id,
               a.accession_tag, a.version_tag,
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
        INNER JOIN accessions a ON j.ship_id = a.id
        LEFT JOIN taxonomy t ON j.tax_id = t.id
        LEFT JOIN starship_features sf ON a.id = sf.accession_id
        LEFT JOIN family_names f ON j.ship_family_id = f.id
        LEFT JOIN navis_names n ON j.ship_navis_id = n.id
        LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
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

            # Deduplicate by accession_tag to prevent multiple rows per sequence
            meta_df = meta_df.drop_duplicates(subset="accession_tag")

            return meta_df
        except Exception as e:
            logger.error(f"Error fetching meta data: {str(e)}")
            raise


@db_retry_decorator()
def test_fetch_ships(
    accession_tags=None, curated=False, dereplicate=True, with_sequence=False
):
    """
    Test-friendly version of fetch_ships without caching.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all ships.
        curated (bool, optional): If True, only fetch curated ships.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing ship data
    """
    with test_db_session() as session:
        query = """
        WITH valid_ships AS (
            SELECT DISTINCT 
                a.id as accession_id, 
                a.accession_tag, a.version_tag,
                CASE 
                    WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
                    THEN a.accession_tag || '.' || a.version_tag
                    ELSE a.accession_tag
                END as accession_display,
                j.curated_status,
                sf.elementBegin, sf.elementEnd, sf.contigID,
                t.name, t.family, t.`order`,
                f.familyName, n.navis_name, h.haplotype_name,
                g.assembly_accession
            FROM joined_ships j
            INNER JOIN accessions a ON j.ship_id = a.id
            LEFT JOIN taxonomy t ON j.tax_id = t.id
            LEFT JOIN family_names f ON j.ship_family_id = f.id
            LEFT JOIN navis_names n ON j.ship_navis_id = n.id
            LEFT JOIN haplotype_names h ON j.ship_haplotype_id = h.id
            LEFT JOIN genomes g ON j.genome_id = g.id
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
                s.md5
            FROM valid_ships v
            LEFT JOIN ships s ON s.accession_id = v.accession_id
            WHERE s.sequence IS NOT NULL
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
                v.assembly_accession
            FROM valid_ships v
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


@db_retry_decorator()
def test_fetch_captains(
    accession_tags=None, curated=False, dereplicate=True, with_sequence=False
):
    """
    Test-friendly version of fetch_captains without caching.

    Args:
        accession_tags (list, optional): List of accession tags to fetch. If None, fetches all captains.
        curated (bool, optional): If True, only fetch curated captains.
        dereplicate (bool, optional): If True, only return one entry per accession tag. Defaults to True.
        with_sequence (bool, optional): If True, fetch sequence data. Defaults to False.
    Returns:
        pd.DataFrame: DataFrame containing captain data
    """
    with test_db_session() as session:
        query = """
        WITH valid_captains AS (
            SELECT DISTINCT 
                a.id,
                a.accession_tag, a.version_tag,
                CASE 
                    WHEN a.version_tag IS NOT NULL AND a.version_tag != '' 
                    THEN a.accession_tag || '.' || a.version_tag
                    ELSE a.accession_tag
                END as accession_display,
                j.curated_status,
                j.starshipID,
                sf.captainID,
                c.sequence,
                n.navis_name
            FROM joined_ships j
            INNER JOIN accessions a ON j.ship_id = a.id
            LEFT JOIN starship_features sf ON a.id = sf.accession_id
            LEFT JOIN captains c ON sf.captainID = c.captainID
            LEFT JOIN navis_names n ON j.ship_navis_id = n.id
            WHERE sf.captainID IS NOT NULL
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
                v.navis_name
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
                v.navis_name
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


# Test-friendly version of the classification workflow
def test_run_classification_workflow(upload_data, meta_dict=None, stages=None):
    """
    Test-friendly version of run_classification_workflow that uses test database functions.

    Args:
        upload_data: Object with fasta_file, seq_type, and fetch parameters
        meta_dict: List of metadata dictionaries
        stages: List of stages to run. If None, runs all stages.
    Returns:
        dict: Classification workflow results
    """
    # Import here to avoid circular imports
    from src.utils.classification_utils import (
        check_exact_match,
        check_contained_match,
        check_similar_match,
        classify_family,
        classify_navis,
        classify_haplotype,
        WORKFLOW_STAGES,
    )

    if stages:
        WORKFLOW_STAGES = [stage for stage in WORKFLOW_STAGES if stage["id"] in stages]

    # Initialize workflow state
    workflow_state = {
        "complete": False,
        "error": None,
        "found_match": False,
        "match_stage": None,
        "match_result": None,
        "stages": {
            stage["id"]: {"progress": 0, "status": "pending"}
            for stage in WORKFLOW_STAGES
        },
        "task_id": "",
        "status": "initialized",
        "workflow_started": True,
        "current_stage": None,
        "current_stage_idx": 0,
        "start_time": 0.0,
        "class_dict": {},
    }

    try:
        # Fetch data using test functions
        logger.info("Fetching ships data for classification...")
        ships_df = test_fetch_ships(
            curated=getattr(upload_data.fetch_ship_params, "curated", False),
            with_sequence=getattr(upload_data.fetch_ship_params, "with_sequence", True),
            dereplicate=getattr(upload_data.fetch_ship_params, "dereplicate", True),
        )

        logger.info("Fetching captains data for classification...")
        captains_df = test_fetch_captains(
            curated=getattr(upload_data.fetch_captain_params, "curated", False),
            with_sequence=getattr(
                upload_data.fetch_captain_params, "with_sequence", True
            ),
            dereplicate=getattr(upload_data.fetch_captain_params, "dereplicate", True),
        )

        # Initialize similarities to None
        similarities = None

        # Run through workflow stages
        for i, stage in enumerate(WORKFLOW_STAGES):
            stage_id = stage["id"]

            # Update state
            workflow_state["current_stage"] = stage_id
            workflow_state["current_stage_idx"] = i
            workflow_state["stages"][stage_id]["progress"] = 10
            workflow_state["stages"][stage_id]["status"] = "running"

            logger.debug(f"Processing stage {i + 1}/{len(WORKFLOW_STAGES)}: {stage_id}")

            if stage_id == "exact":
                logger.debug("Running exact match check")
                result = check_exact_match(
                    fasta=upload_data.fasta_file, existing_ships=ships_df
                )

                if result:
                    logger.debug(f"Found exact match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "exact"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            elif stage_id == "contained":
                logger.debug("Running contained match check")
                result = check_contained_match(
                    fasta=upload_data.fasta_file,
                    existing_ships=ships_df,
                    min_coverage=0.95,
                    min_identity=0.95,
                )

                if result:
                    logger.debug(f"Found contained match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 30
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "contained"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            elif stage_id == "similar":
                logger.debug("Running similarity match check")
                result, similarities = check_similar_match(
                    fasta=upload_data.fasta_file,
                    existing_ships=ships_df,
                    threshold=0.95,
                )

                if result:
                    logger.debug(f"Found similar match: {result}")
                    workflow_state["stages"][stage_id]["progress"] = 50
                    workflow_state["stages"][stage_id]["status"] = "complete"
                    workflow_state["found_match"] = True
                    workflow_state["match_stage"] = "similar"
                    workflow_state["match_result"] = result
                    workflow_state["complete"] = True
                    return workflow_state

            elif stage_id == "family":
                logger.debug("Running family classification")
                family_dict, protein_file = classify_family(
                    fasta=upload_data.fasta_file,
                    seq_type=upload_data.seq_type,
                    meta_dict=meta_dict,
                    pident_thresh=90,
                    input_eval=0.001,
                    threads=1,
                )

                if family_dict:
                    family_name = family_dict["family"]
                    logger.debug(f"Found family classification: {family_name}")
                    workflow_state["stages"][stage_id]["progress"] = 70
                    workflow_state["stages"][stage_id]["status"] = "complete"

                    # Store family result but don't return yet - continue to other stages
                    if not workflow_state["found_match"]:
                        workflow_state["found_match"] = True
                        workflow_state["match_stage"] = "family"

                        # Simplify the result
                        if isinstance(family_dict, dict) and "family" in family_dict:
                            workflow_state["match_result"] = family_dict["family"]
                        else:
                            workflow_state["match_result"] = family_dict

                    # Store family result for testing purposes
                    workflow_state["class_dict"]["family"] = family_name
                else:
                    logger.debug("No family classification found")
                    workflow_state["stages"][stage_id]["progress"] = 100
                    workflow_state["stages"][stage_id]["status"] = "complete"

                    # Early stopping: If no family classification found, skip navis and haplotype
                    logger.debug(
                        "Skipping navis and haplotype classification due to family failure"
                    )
                    for remaining_stage in ["navis", "haplotype"]:
                        if remaining_stage in workflow_state["stages"]:
                            workflow_state["stages"][remaining_stage]["progress"] = 100
                            workflow_state["stages"][remaining_stage]["status"] = (
                                "skipped"
                            )

                    # Complete the workflow
                    workflow_state["complete"] = True
                    return workflow_state

            elif stage_id == "navis":
                logger.debug("Running navis classification")
                if captains_df.empty:
                    logger.warning("No captain data available for navis classification")
                    workflow_state["stages"][stage_id]["progress"] = 80
                    workflow_state["stages"][stage_id]["status"] = "skipped"

                    # Early stopping: If no captain data, skip haplotype
                    logger.debug(
                        "Skipping haplotype classification due to no captain data"
                    )
                    if "haplotype" in workflow_state["stages"]:
                        workflow_state["stages"]["haplotype"]["progress"] = 100
                        workflow_state["stages"]["haplotype"]["status"] = "skipped"

                    # Complete the workflow
                    workflow_state["complete"] = True
                    return workflow_state
                else:
                    result = classify_navis(
                        fasta=upload_data.fasta_file,
                        existing_captains=captains_df,
                        threads=1,
                    )

                    if result:
                        logger.debug(f"Found navis classification: {result}")
                        workflow_state["stages"][stage_id]["progress"] = 90
                        workflow_state["stages"][stage_id]["status"] = "complete"

                        # Store navis result but don't return yet - continue to other stages
                        if not workflow_state["found_match"]:
                            workflow_state["found_match"] = True
                            workflow_state["match_stage"] = "navis"
                            workflow_state["match_result"] = result

                        # Store navis result for testing purposes
                        workflow_state["class_dict"]["navis"] = result
                    else:
                        logger.debug("No navis classification found")
                        workflow_state["stages"][stage_id]["progress"] = 90
                        workflow_state["stages"][stage_id]["status"] = "complete"

                        # Early stopping: If no navis classification found, skip haplotype
                        logger.debug(
                            "Skipping haplotype classification due to navis failure"
                        )
                        if "haplotype" in workflow_state["stages"]:
                            workflow_state["stages"]["haplotype"]["progress"] = 100
                            workflow_state["stages"]["haplotype"]["status"] = "skipped"

                        # Complete the workflow
                        workflow_state["complete"] = True
                        return workflow_state

            elif stage_id == "haplotype":
                logger.debug("Running haplotype classification")
                if captains_df.empty or ships_df.empty:
                    logger.warning("Missing data for haplotype classification")
                    workflow_state["stages"][stage_id]["progress"] = 90
                    workflow_state["stages"][stage_id]["status"] = "skipped"
                else:
                    # Extract navis value
                    navis_value = None
                    try:
                        if (
                            not captains_df.empty
                            and "navis_name" in captains_df.columns
                        ):
                            navis_values = captains_df["navis_name"].dropna().unique()
                            if len(navis_values) > 0:
                                navis_value = navis_values[0]

                        if navis_value is None and "navis_name" in ships_df.columns:
                            navis_values = ships_df["navis_name"].dropna().unique()
                            if len(navis_values) > 0:
                                navis_value = navis_values[0]
                    except Exception as e:
                        logger.error(f"Error extracting navis value: {e}")

                    if navis_value is None:
                        navis_value = "UNK"

                    try:
                        result = classify_haplotype(
                            fasta=upload_data.fasta_file,
                            existing_ships=ships_df,
                            navis=navis_value,
                            similarities=similarities,
                        )

                        if result:
                            logger.debug(f"Found haplotype classification: {result}")
                            workflow_state["stages"][stage_id]["progress"] = 100
                            workflow_state["stages"][stage_id]["status"] = "complete"

                            # Store haplotype result but don't return yet - continue to other stages
                            if not workflow_state["found_match"]:
                                workflow_state["found_match"] = True
                                workflow_state["match_stage"] = "haplotype"
                                workflow_state["match_result"] = result

                            # Store haplotype result for testing purposes
                            if isinstance(result, dict) and "haplotype_name" in result:
                                workflow_state["class_dict"]["haplotype"] = result[
                                    "haplotype_name"
                                ]
                            else:
                                workflow_state["class_dict"]["haplotype"] = result
                        else:
                            logger.debug("No haplotype classification found")
                            workflow_state["stages"][stage_id]["progress"] = 100
                            workflow_state["stages"][stage_id]["status"] = "complete"
                    except Exception as e:
                        logger.error(f"Error in haplotype classification: {e}")
                        workflow_state["stages"][stage_id]["status"] = "error"
                        workflow_state["error"] = (
                            f"Haplotype classification error: {str(e)}"
                        )

            # Mark stage as complete
            workflow_state["stages"][stage_id]["progress"] = 100
            workflow_state["stages"][stage_id]["status"] = "complete"

        # Try BLAST fallback if no matches found
        if (
            not workflow_state.get("found_match", False)
            and hasattr(upload_data, "blast_df")
            and upload_data.blast_df is not None
        ):
            logger.debug("Trying BLAST fallback")
            try:
                blast_df = upload_data.blast_df
                if isinstance(blast_df, list):
                    blast_df = pd.DataFrame(blast_df)

                if not blast_df.empty:
                    blast_df = blast_df.sort_values(
                        ["evalue", "pident"], ascending=[True, False]
                    )
                    top_hit = blast_df.iloc[0]

                    top_pident = float(top_hit["pident"])
                    if top_pident >= 90:
                        hit_IDs = top_hit["hit_IDs"]
                        hit_IDs_list = (
                            [hit_IDs] if isinstance(hit_IDs, str) else hit_IDs
                        )

                        if meta_dict:
                            meta_df = pd.DataFrame(meta_dict)
                            meta_df_sub = meta_df[
                                meta_df["accession_tag"].isin(hit_IDs_list)
                            ]

                            if not meta_df_sub.empty:
                                top_family = meta_df_sub["familyName"].iloc[0]
                                workflow_state["found_match"] = True
                                workflow_state["match_stage"] = "blast_hit"
                                workflow_state["match_result"] = {
                                    "source": "blast_hit",
                                    "family": top_family,
                                    "closest_match": hit_IDs,
                                    "confidence": "High"
                                    if top_pident >= 90
                                    else "Medium",
                                }
            except Exception as e:
                logger.error(f"Error processing BLAST fallback: {e}")

        # Mark workflow as complete
        workflow_state["complete"] = True
        return workflow_state

    except Exception as e:
        error_message = str(e)
        logger.error(f"Error in classification workflow: {error_message}")
        workflow_state["error"] = error_message
        workflow_state["complete"] = True
        return workflow_state

import traceback
import pandas as pd

from src.utils.seq_utils import extract_accession
from src.database.sql_manager import (
    fetch_meta_data,
    get_quality_tags,
)

from src.config.logging import get_logger

logger = get_logger(__name__)


def safe_get_value(df, column, index=0, default="N/A", format_func=None):
    """
    Safely extract a value from a DataFrame, handling NA/Null values.

    Args:
        df: DataFrame to extract from
        column: Column name
        index: Row index (default 0)
        default: Default value if NA/Null (default "N/A")
        format_func: Optional function to format the value

    Returns:
        Formatted value or default
    """
    try:
        if column not in df.columns:
            return default

        value = df[column].iloc[index]

        # Check for various forms of null/empty values
        if (
            pd.isna(value)
            or value is None
            or value == ""
            or str(value).lower() in ["nan", "none", "null"]
        ):
            return default

        # Apply formatting function if provided
        if format_func and callable(format_func):
            try:
                return format_func(value)
            except (ValueError, TypeError):
                return default

        return str(value)

    except (IndexError, KeyError):
        return default


def safe_get_numeric(df, column, index=0, default="N/A"):
    """Helper for numeric values"""
    return safe_get_value(df, column, index, default, lambda x: str(int(float(x))))


def safe_get_position(df, begin_col, end_col, index=0, default="N/A"):
    """Helper for position ranges"""
    begin = safe_get_value(df, begin_col, index, None, lambda x: int(float(x)))
    end = safe_get_value(df, end_col, index, None, lambda x: int(float(x)))

    if begin != "N/A" and end != "N/A" and begin is not None and end is not None:
        return f"{begin} - {end}"
    return default


def create_ship_accession_modal_data(ship_accession_id):
    """Create structured data for accession modal instead of Dash components."""
    try:
        base_accession = extract_accession(ship_accession_id)

        modal_data = fetch_meta_data(accessions=[base_accession])

        if modal_data.empty:
            return {
                "title": f"Accession: {ship_accession_id}",
                "error": f"No data found for accession: {ship_accession_id}",
            }

        # Validate modal_data
        if not isinstance(modal_data, pd.DataFrame) or modal_data.empty:
            return {
                "title": f"Accession: {ship_accession_id}",
                "error": "Invalid modal data received",
            }

        # Check for required columns
        required_columns = ["accession_tag", "familyName"]
        missing_columns = [
            col for col in required_columns if col not in modal_data.columns
        ]
        if missing_columns:
            return {
                "title": f"Accession: {ship_accession_id}",
                "error": f"Missing required columns: {missing_columns}",
            }

        # Fetch quality tags separately using joined_ship_id
        joined_ship_id = safe_get_numeric(modal_data, "joined_ship_id")
        quality_tags = []
        accepted_quality_tags = [
            "missing_direct_repeats",
            "missing_tir",
            "missing_boundaries",
            "missing_genome_context",
            "unannotated",
            "missing_empty_site",
        ]
        if joined_ship_id:
            try:
                quality_tags_data = get_quality_tags(joined_ship_id)
                # Format tags as "tag_type:tag_value" or just "tag_type" if no value
                for tag in quality_tags_data:
                    if tag.get("tag_value"):
                        tag_value = f"{tag['tag_type']}:{tag['tag_value']}"
                    else:
                        tag_value = tag["tag_type"]
                    if tag_value in accepted_quality_tags:
                        quality_tags.append(tag_value)
            except Exception as e:
                logger.warning(
                    f"Error fetching quality tags for joined_ship_id {joined_ship_id}: {e}"
                )

        # Create structured data
        result = {
            "title": base_accession,
            "version_tag": safe_get_value(modal_data, "version_tag"),
            "familyName": safe_get_value(modal_data, "familyName"),
            "navis_name": safe_get_value(modal_data, "navis_name"),
            "haplotype_name": safe_get_value(modal_data, "haplotype_name"),
            "taxonomic_family": safe_get_value(modal_data, "family"),
            "order": safe_get_value(modal_data, "order"),
            "species_name": safe_get_value(modal_data, "name"),
            "tax_id": safe_get_numeric(modal_data, "taxID"),
            "assembly_accession": safe_get_value(
                modal_data, "assembly_accession", default=""
            ),
            "genome_source": safe_get_value(modal_data, "genomeSource", default=""),
            "contig_id": safe_get_value(modal_data, "contigID", default=""),
            "element_length": safe_get_numeric(modal_data, "elementLength", default=""),
            "element_position": safe_get_position(
                modal_data, "elementBegin", "elementEnd"
            ),
            "curated_status": safe_get_value(
                modal_data, "curated_status", default="unknown"
            ),
            "quality_tags": quality_tags,
        }

        return result

    except Exception as e:
        logger.error(f"Error in create_ship_accession_modal_data: {str(e)}")
        logger.error(traceback.format_exc())
        raise


def create_accession_modal_data(accession):
    """Create structured data for accession modal instead of Dash components."""
    try:
        base_accession = extract_accession(accession)

        modal_data = fetch_meta_data(accessions=[base_accession])

        if modal_data.empty:
            return {
                "title": f"Accession: {accession}",
                "error": f"No data found for accession: {accession}",
            }

        # Validate modal_data
        if not isinstance(modal_data, pd.DataFrame) or modal_data.empty:
            return {
                "title": f"Accession: {accession}",
                "error": "Invalid modal data received",
            }

        # Check for required columns
        required_columns = ["accession_tag", "familyName"]
        missing_columns = [
            col for col in required_columns if col not in modal_data.columns
        ]
        if missing_columns:
            return {
                "title": f"Accession: {accession}",
                "error": f"Missing required columns: {missing_columns}",
            }

        # HACK: applying a fix for extra rows in the starship_features table, only take the first begin/end coordinates for each ship_id/accession_id
        # ! this might cause some issues if coordinates are not updated for all rows for a ship_id/accession_id pair, updated only if begin/end coordinates are the same
        # TODO: split features table or move coordinate information to separate table or another existing table
        modal_data = modal_data.groupby("accession_tag").first().reset_index()

        # Fetch quality tags separately using joined_ship_id
        joined_ship_id = safe_get_numeric(modal_data, "joined_ship_id")
        quality_tags = []
        accepted_quality_tags = [
            "missing_direct_repeats",
            "missing_tir",
            "missing_boundaries",
            "missing_genome_context",
            "unannotated",
            "missing_empty_site",
        ]
        if joined_ship_id:
            try:
                quality_tags_data = get_quality_tags(joined_ship_id)
                # Format tags as "tag_type:tag_value" or just "tag_type" if no value
                for tag in quality_tags_data:
                    if tag.get("tag_value"):
                        tag_value = f"{tag['tag_type']}:{tag['tag_value']}"
                    else:
                        tag_value = tag["tag_type"]
                    if tag_value in accepted_quality_tags:
                        quality_tags.append(tag_value)
            except Exception as e:
                logger.warning(
                    f"Error fetching quality tags for joined_ship_id {joined_ship_id}: {e}"
                )

        # TODO: Create more comprehensive structured data
        # - some output will be the same across all ships within this accession
        # - some output we will have to aggregate across all ships within this accession
        result = {
            "title": f"Starship Accession: {accession}",
            "familyName": safe_get_value(modal_data, "familyName"),
            "genomes_present": str(len(modal_data)),
            "navis_name": safe_get_value(modal_data, "navis_name"),
            "haplotype_name": safe_get_value(modal_data, "haplotype_name"),
        }

        return result

    except Exception as e:
        logger.error(f"Error in create_accession_modal_data: {str(e)}")
        logger.error(traceback.format_exc())
        raise

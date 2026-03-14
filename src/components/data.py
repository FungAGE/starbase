import traceback
import pandas as pd

from src.utils.seq_utils import extract_accession, clean_contigIDs
from src.database.sql_manager import (
    fetch_meta_data,
    get_quality_tags,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

ACCEPTED_QUALITY_TAGS = [
    "missing_direct_repeats",
    "missing_tir",
    "missing_boundaries",
    "missing_genome_context",
    "unannotated",
    "missing_empty_site",
]


def _fetch_quality_tags_for_modal(joined_ship_id):
    """Fetch and format quality tags for modal display. Returns list of tag strings."""
    if not joined_ship_id:
        return []
    try:
        quality_tags_data = get_quality_tags(joined_ship_id)
        tags = []
        for tag in quality_tags_data:
            if tag.get("tag_value"):
                tag_value = f"{tag['tag_type']}:{tag['tag_value']}"
            else:
                tag_value = tag["tag_type"]
            if tag_value in ACCEPTED_QUALITY_TAGS:
                tags.append(tag_value)
        return tags
    except Exception as e:
        logger.warning(
            f"Error fetching quality tags for joined_ship_id {joined_ship_id}: {e}"
        )
        return []


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


def _get_value_if_active(modal_data, value_col, activity_col, index=0):
    """
    Return value only if activity is not 0. NULL/None activity = active (show).
    When activity=0, return None so the modal omits the field.
    """
    if activity_col in modal_data.columns:
        try:
            act = modal_data[activity_col].iloc[index]
            if not pd.isna(act) and act is not None and int(act) == 0:
                return None  # inactive - don't show
        except (IndexError, KeyError, ValueError, TypeError):
            pass
    return safe_get_value(modal_data, value_col, index, default=None)


def _modal_element_length(modal_data, index=0):
    """
    Element length for modals: always use computed value
    elementEnd - elementBegin + 1 when both coordinates exist (canonical).
    Only use starship_features.elementLength when coordinates are missing.
    Never use ships.sequence_length.
    """
    begin = safe_get_value(
        modal_data, "elementBegin", index, None, lambda x: int(float(x))
    )
    end = safe_get_value(modal_data, "elementEnd", index, None, lambda x: int(float(x)))
    if begin not in (None, "N/A") and end not in (None, "N/A"):
        return str(abs(int(end) - int(begin) + 1))
    length = safe_get_numeric(modal_data, "elementLength", index=index, default="")
    return length if length and length != "" else ""


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

        joined_ship_id = safe_get_numeric(modal_data, "joined_ship_id")
        quality_tags = _fetch_quality_tags_for_modal(joined_ship_id)

        # Build genomes array: one entry per genome this ship is present in
        genomes = []
        for i in range(len(modal_data)):
            genome = {
                "assembly_accession": safe_get_value(
                    modal_data, "assembly_accession", index=i, default=""
                ),
                "genome_source": safe_get_value(
                    modal_data, "genomeSource", index=i, default=""
                ),
                "contig_id": clean_contigIDs(
                    safe_get_value(modal_data, "contigID", index=i, default="")
                ) or "",
                "element_length": _modal_element_length(modal_data, index=i),
                "element_position": safe_get_position(
                    modal_data, "elementBegin", "elementEnd", index=i
                ),
            }
            genomes.append(genome)

        # Create structured data (omit navis/haplotype when activity=0)
        result = {
            "title": base_accession,
            "version_tag": safe_get_value(modal_data, "version_tag"),
            "familyName": safe_get_value(modal_data, "familyName"),
            "navis_name": _get_value_if_active(
                modal_data, "navis_name", "navis_activity"
            ),
            "haplotype_name": _get_value_if_active(
                modal_data, "haplotype_name", "haplotype_activity"
            ),
            "family": safe_get_value(modal_data, "family"),
            "order": safe_get_value(modal_data, "order"),
            "species_name": safe_get_value(modal_data, "name"),
            "tax_id": safe_get_numeric(modal_data, "taxID"),
            "genomes_present": str(len(modal_data)),
            "genomes": genomes,
            "curated_status": safe_get_value(
                modal_data, "curated_status", default="uncurated"
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
        title = f"Group Accession: {base_accession}"

        modal_data = fetch_meta_data(accessions=[base_accession])

        if modal_data.empty:
            return {
                "title": title,
                "error": f"No data found for accession: {accession}",
            }

        if not isinstance(modal_data, pd.DataFrame) or modal_data.empty:
            return {
                "title": title,
                "error": "Invalid modal data received",
            }

        required_columns = ["accession_tag", "familyName"]
        missing_columns = [
            col for col in required_columns if col not in modal_data.columns
        ]
        if missing_columns:
            return {
                "title": title,
                "error": f"Missing required columns: {missing_columns}",
            }

        # HACK: applying a fix for extra rows in the starship_features table, only take the first begin/end coordinates for each ship_id/accession_id
        # ! this might cause some issues if coordinates are not updated for all rows for a ship_id/accession_id pair, updated only if begin/end coordinates are the same
        # TODO: split features table or move coordinate information to separate table or another existing table
        modal_data = modal_data.groupby("accession_tag").first().reset_index()

        # TODO: Create more comprehensive structured data
        # - some output will be the same across all ships within this accession
        # - some output we will have to aggregate across all ships within this accession

        result = {
            "title": title,
            "familyName": safe_get_value(modal_data, "familyName"),
            "genomes_present": str(len(modal_data)),
            "navis_name": _get_value_if_active(
                modal_data, "navis_name", "navis_activity"
            ),
            "haplotype_name": _get_value_if_active(
                modal_data, "haplotype_name", "haplotype_activity"
            ),
        }

        return result

    except Exception as e:
        logger.error(f"Error in create_accession_modal_data: {str(e)}")
        logger.error(traceback.format_exc())
        raise

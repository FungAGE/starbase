import re
import traceback
from urllib.parse import quote

import pandas as pd
import dash_mantine_components as dmc
from dash import html

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
        accessions_to_try = list({str(ship_accession_id).strip(), base_accession})

        modal_data = fetch_meta_data(accessions=accessions_to_try)

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


# --- Native dmc.Modal rendering (replaces the vanilla-JS universal-modal renderers) ---

_ACCESSION_RE = re.compile(r"^[A-Za-z]{2,}_?\d{6,}(\.\d+)?$", re.IGNORECASE)
_POSITION_RE = re.compile(r"^(\d+)\s*-\s*(\d+)$")


def _has_value(v):
    return v is not None and v != "" and str(v).strip().upper() != "N/A"


def _is_valid_accession(value):
    return bool(value) and isinstance(value, str) and bool(_ACCESSION_RE.match(value.strip()))


def _genome_urls(genome):
    """Mirror the URL-construction logic from universal-modal.js's renderAccessionModal."""
    assembly = genome.get("assembly_accession")
    source = (genome.get("genome_source") or "").lower()
    contig_id = genome.get("contig_id")
    position = genome.get("element_position")

    genome_url = None
    sequence_viewer_url = None

    if _has_value(assembly) and _has_value(source) and _is_valid_accession(assembly):
        if source == "jgi":
            genome_url = f"https://mycocosm.jgi.doe.gov/{quote(assembly)}/{quote(assembly)}.home.html"
        elif source == "ncbi":
            genome_url = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{quote(assembly)}/"

        if (
            source in ("jgi", "ncbi")
            and _has_value(contig_id)
            and _is_valid_accession(contig_id)
            and _has_value(position)
        ):
            match = _POSITION_RE.match(position)
            if match:
                sequence_viewer_url = (
                    "https://www.ncbi.nlm.nih.gov/projects/sviewer/"
                    f"?id={quote(contig_id)}&from={match.group(1)}&to={match.group(2)}"
                )

    return genome_url, sequence_viewer_url


def _modal_row(label, value):
    return dmc.Group(
        [
            dmc.Text(label, fw=600, size="sm", c="dimmed"),
            value if hasattr(value, "to_plotly_json") else dmc.Text(str(value), size="sm"),
        ],
        justify="space-between",
        wrap="nowrap",
        gap="md",
    )


def _modal_section(title, children):
    return dmc.Paper(
        children=[
            dmc.Text(title, fw=700, size="lg", ta="center", mb="md"),
            dmc.Stack(children, gap="xs"),
        ],
        p="lg",
        radius="md",
        withBorder=True,
    )


def render_ship_accession_modal(data):
    """Build the dmc content for a ship-accession modal from create_ship_accession_modal_data()'s output."""
    if data.get("error"):
        return dmc.Alert(data["error"], title="Error", color="red", variant="light")

    genomes = data.get("genomes") or []
    starship_size_bp = None
    if len(genomes) == 1 and _has_value(genomes[0].get("element_length")):
        starship_size_bp = str(genomes[0]["element_length"])

    info_rows = []
    if starship_size_bp:
        info_rows.append(_modal_row("Size", f"{starship_size_bp} bp"))
    if _has_value(data.get("familyName")):
        info_rows.append(_modal_row("Starship Family", data["familyName"]))
    if _has_value(data.get("navis_name")):
        info_rows.append(_modal_row("Starship Navis", data["navis_name"]))
    if _has_value(data.get("haplotype_name")):
        info_rows.append(_modal_row("Haplotype", data["haplotype_name"]))
    genomes_present = data.get("genomes_present")
    if _has_value(genomes_present) and int(genomes_present) > 1:
        info_rows.append(
            _modal_row("Genomes Present", dmc.Badge(genomes_present, color="blue", variant="light"))
        )

    genome_sections = []
    for i, genome in enumerate(genomes):
        assembly = genome.get("assembly_accession")
        source = genome.get("genome_source")
        contig_id = genome.get("contig_id")
        position = genome.get("element_position")
        length = genome.get("element_length")

        has_content = any(_has_value(v) for v in (assembly, source, contig_id, position)) or (
            not starship_size_bp and _has_value(length)
        )
        if not has_content:
            continue

        genome_url, sequence_viewer_url = _genome_urls(genome)

        if _has_value(assembly) and _has_value(source):
            section_title = ["Genome: ", dmc.Anchor(assembly, href=genome_url, target="_blank")] if genome_url else f"Genome: {assembly}"
        elif _has_value(assembly):
            section_title = assembly
        elif _has_value(source):
            section_title = source
        else:
            section_title = "Genome" if len(genomes) == 1 else f"Genome {i + 1}"

        rows = []
        if _has_value(assembly):
            value = dmc.Anchor(assembly, href=genome_url, target="_blank") if genome_url else assembly
            rows.append(_modal_row("Assembly Accession", value))
        if _has_value(source):
            rows.append(_modal_row("Genome Source", source))
        if _has_value(contig_id):
            value = dmc.Anchor(contig_id, href=sequence_viewer_url, target="_blank") if sequence_viewer_url else contig_id
            rows.append(_modal_row("Contig ID", value))
        if _has_value(position):
            rows.append(_modal_row("Element Position", position))
        if not starship_size_bp and _has_value(length):
            rows.append(_modal_row("Size", f"{length} bp"))

        genome_sections.append(
            dmc.Paper(
                children=[
                    dmc.Text(section_title, fw=700, size="md", ta="center", mb="sm"),
                    dmc.Stack(rows, gap="xs"),
                ],
                p="md",
                radius="md",
                withBorder=True,
                mt="md",
            )
        )

    sections = []
    if info_rows or genome_sections:
        sections.append(
            _modal_section(
                "Starship Information",
                (info_rows if info_rows else []) + genome_sections,
            )
        )

    taxonomy_rows = []
    if _has_value(data.get("order")):
        taxonomy_rows.append(_modal_row("Order", data["order"]))
    if _has_value(data.get("family")):
        taxonomy_rows.append(_modal_row("Family", data["family"]))
    if _has_value(data.get("species_name")):
        taxonomy_rows.append(_modal_row("Species", data["species_name"]))
    if _has_value(data.get("tax_id")):
        taxonomy_rows.append(
            _modal_row(
                "NCBI Taxonomy ID",
                dmc.Anchor(
                    str(data["tax_id"]),
                    href=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={quote(str(data['tax_id']))}",
                    target="_blank",
                ),
            )
        )
    if taxonomy_rows:
        sections.append(_modal_section("Taxonomy", taxonomy_rows))

    quality_rows = []
    if _has_value(data.get("curated_status")):
        color = "green" if data["curated_status"] == "curated" else "yellow"
        quality_rows.append(
            _modal_row("Curation Status", dmc.Badge(data["curated_status"], color=color, variant="light"))
        )
    quality_tags = data.get("quality_tags") or []
    if quality_tags:
        quality_rows.append(
            _modal_row(
                "Quality Tags",
                dmc.Group(
                    [dmc.Badge(tag, color="yellow", variant="light") for tag in quality_tags],
                    gap="xs",
                    justify="flex-end",
                ),
            )
        )
    if quality_rows:
        sections.append(_modal_section("Data Quality", quality_rows))

    if not sections:
        return dmc.Text("No details available for this accession.", c="dimmed", ta="center")

    return dmc.Stack(sections, gap="md")


def render_group_accession_modal(data):
    """Build the dmc content for a group-accession (SSA) modal from create_accession_modal_data()'s output."""
    if data.get("error"):
        return dmc.Alert(data["error"], title="Error", color="red", variant="light")

    rows = []
    if _has_value(data.get("familyName")):
        rows.append(_modal_row("Starship Family", data["familyName"]))
    if _has_value(data.get("genomes_present")):
        rows.append(
            _modal_row("Genomes Present", dmc.Badge(data["genomes_present"], color="blue", variant="light"))
        )
    if _has_value(data.get("navis_name")):
        rows.append(_modal_row("Starship Navis", data["navis_name"]))
    if _has_value(data.get("haplotype_name")):
        rows.append(_modal_row("Haplotype", data["haplotype_name"]))

    if not rows:
        return dmc.Text("No details available for this accession.", c="dimmed", ta="center")

    return _modal_section("Group Information", rows)

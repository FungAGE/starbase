"""Editable location metadata grid for the submit page."""

from typing import Any, Dict, List, Optional

import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash import dcc

from src.utils.web_submission_adapter import infer_metadata_strand

LOCATION_GRID_ID = "submit-location-grid"

# ship_header: read-only, populated from FASTA
# genus, species, hostchr, shipstart, shipend: required before submit (may be filled via TSV)
# strain, assembly_accession, strand: optional (strand inferred from start/end when omitted)
EDITABLE_COLUMNS = {
    "genus",
    "species",
    "strain",
    "assembly_accession",
    "hostchr",
    "shipstart",
    "shipend",
    "strand",
}

COLUMN_DEFS = [
    {
        "field": "ship_header",
        "headerName": "Ship header",
        "editable": False,
        "pinned": "left",
        "minWidth": 140,
        "flex": 1,
    },
    {
        "field": "genus",
        "headerName": "Genus",
        "editable": True,
        "minWidth": 110,
    },
    {
        "field": "species",
        "headerName": "Species",
        "editable": True,
        "minWidth": 110,
    },
    {
        "field": "strain",
        "headerName": "Strain",
        "editable": True,
        "minWidth": 90,
    },
    {
        "field": "assembly_accession",
        "headerName": "Genome Accession",
        "editable": True,
        "minWidth": 140,
    },
    {
        "field": "hostchr",
        "headerName": "Host contig",
        "editable": True,
        "minWidth": 110,
    },
    {
        "field": "shipstart",
        "headerName": "Start",
        "editable": True,
        "minWidth": 90,
    },
    {
        "field": "shipend",
        "headerName": "End",
        "editable": True,
        "minWidth": 90,
    },
    {
        "field": "strand",
        "headerName": "Strand",
        "editable": True,
        "minWidth": 80,
        "cellEditor": "agSelectCellEditor",
        "cellEditorParams": {"values": ["+", "-"]},
    },
]


def empty_location_row(ship_header: str) -> Dict[str, Any]:
    return {
        "ship_header": ship_header,
        "genus": "",
        "species": "",
        "strain": "",
        "assembly_accession": "",
        "hostchr": "",
        "shipstart": "",
        "shipend": "",
        "strand": "+",
    }


def init_location_rows(
    seq_ids: List[str], existing: Optional[List[Dict[str, Any]]] = None
) -> List[Dict[str, Any]]:
    """Build grid rows for FASTA headers, preserving edits for matching ship headers."""
    existing_by_header = {}
    for row in existing or []:
        key = row.get("ship_header") or row.get("sequence")
        if key:
            existing_by_header[key] = row

    rows = []
    for seq_id in seq_ids:
        if seq_id in existing_by_header:
            row = dict(existing_by_header[seq_id])
            row["ship_header"] = seq_id
            rows.append(row)
        else:
            rows.append(empty_location_row(seq_id))
    return rows


def grid_height(row_count: int) -> str:
    if row_count <= 0:
        return "200px"
    height = min(600, max(200, 44 + 40 * row_count))
    return f"{height}px"


def create_location_grid(row_data: Optional[List[Dict[str, Any]]] = None) -> dag.AgGrid:
    rows = row_data or []
    dash_grid_options = {
        "suppressPropertyNamesCheck": True,
        "rowHeight": 40,
        "headerHeight": 44,
        "enableCellTextSelection": True,
        "ensureDomOrder": True,
    }

    col_defs = []
    for col_def in COLUMN_DEFS:
        updated = dict(col_def)
        if updated["field"] in EDITABLE_COLUMNS:
            updated["cellStyle"] = {
                "backgroundColor": "var(--mantine-color-gray-0)",
            }
        col_defs.append(updated)

    return dag.AgGrid(
        id=LOCATION_GRID_ID,
        columnDefs=col_defs,
        rowData=rows,
        defaultColDef={"resizable": True, "minWidth": 80},
        dashGridOptions=dash_grid_options,
        getRowId="params.data.ship_header",
        className="ag-theme-alpine",
        style={"width": "100%", "height": grid_height(len(rows))},
    )


def create_location_status(seq_ids: Optional[List[str]] = None):
    if not seq_ids:
        return dmc.Text(
            "Upload a FASTA file to fill in organism and location details below.",
            size="sm",
            c="dimmed",
        )

    label = f"{len(seq_ids)} ship{'s' if len(seq_ids) != 1 else ''} — edit the table or import a TSV/CSV."
    if len(seq_ids) > 1:
        label += " GFF rows match FASTA headers automatically."
    return dmc.Text(label, size="sm", c="dimmed")


def create_location_tsv_upload():
    return dmc.Group(
        [
            dcc.Upload(
                id="submit-location-tsv-upload",
                children=dmc.Button(
                    "Import TSV/CSV",
                    variant="light",
                    color="indigo",
                    size="sm",
                ),
                accept=".tsv,.csv,.txt",
                multiple=False,
            ),
            dmc.Accordion(
                variant="contained",
                chevronPosition="right",
                style={"flex": 1},
                children=[
                    dmc.AccordionItem(
                        value="format",
                        children=[
                            dmc.AccordionControl(
                                dmc.Text("File format", size="sm", c="dimmed"),
                            ),
                            dmc.AccordionPanel(
                                dmc.Stack(
                                    [
                                        dmc.Text(
                                            "Rows match by ship_header. Include any columns you want to fill; "
                                            "required before submit: genus, species, hostchr, shipstart, shipend.",
                                            size="sm",
                                            c="dimmed",
                                        ),
                                        dmc.Code(
                                            "ship_header\tgenus\tspecies\tstrain\tassembly_accession\thostchr\tshipstart\tshipend\n"
                                            "seq_1\tAlternaria\talternata\t\tGCA_000001305.1\tchr1\t1000\t5000",
                                            block=True,
                                            className="submit-tsv-example",
                                        ),
                                    ],
                                    gap="xs",
                                ),
                            ),
                        ],
                    ),
                ],
            ),
        ],
        gap="sm",
        align="flex-start",
        wrap="wrap",
        className="submit-tsv-toolbar",
    )


def _normalize_strand(strand: Any, shipstart: Any = None, shipend: Any = None) -> int:
    """Return strand_radio (1=forward, 2=reverse), inferring from coordinates when needed."""
    symbol = infer_metadata_strand(strand, shipstart, shipend)
    return 2 if symbol == "-" else 1


def _parse_coordinate(value: Any) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def rows_to_locations(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Convert grid row data to the locations list expected by build_submission_entries."""
    locations = []
    for row in rows:

        def _text(key: str):
            value = row.get(key)
            if value is None:
                return None
            if isinstance(value, str):
                value = value.strip()
            return value or None

        locations.append(
            {
                "genus": _text("genus"),
                "species": _text("species"),
                "strain": _text("strain"),
                "assembly_accession": _text("assembly_accession"),
                "hostchr": _text("hostchr"),
                "shipstart": _parse_coordinate(row.get("shipstart")),
                "shipend": _parse_coordinate(row.get("shipend")),
                "strand_radio": _normalize_strand(
                    row.get("strand"),
                    row.get("shipstart"),
                    row.get("shipend"),
                ),
            }
        )
    return locations

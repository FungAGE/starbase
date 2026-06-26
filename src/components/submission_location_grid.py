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
                "backgroundColor": "var(--mantine-color-yellow-0)",
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
            "Upload a FASTA file above to enter organism and location details for each sequence.",
            size="sm",
            c="dimmed",
        )

    children = [
        dmc.Text(
            "Starship submissions",
            size="sm",
            fw=600,
            c="dimmed",
            tt="uppercase",
        ),
        dmc.Text(
            f"{len(seq_ids)} sequence{'s' if len(seq_ids) != 1 else ''} — enter organism and location details for each row.",
            size="sm",
            c="dimmed",
        ),
    ]
    if len(seq_ids) > 1:
        children.append(
            dmc.Alert(
                "Annotations in your GFF file are matched to FASTA headers automatically.",
                color="var(--mantine-color-blue-6)",
                variant="light",
                title="Multiple sequences",
            )
        )
    return dmc.Stack(children, gap="xs")


def create_location_tsv_upload():
    return dmc.Stack(
        [
            dmc.Group(
                [
                    dmc.Text("Import metadata", size="sm", fw=500),
                    dcc.Upload(
                        id="submit-location-tsv-upload",
                        children=dmc.Button(
                            "Upload TSV/CSV",
                            variant="light",
                            color="indigo",
                            size="xs",
                        ),
                        accept=".tsv,.csv,.txt",
                        multiple=False,
                    ),
                ],
                gap="sm",
                align="center",
            ),
            dmc.Accordion(
                variant="contained",
                chevronPosition="right",
                children=[
                    dmc.AccordionItem(
                        value="format",
                        children=[
                            dmc.AccordionControl(
                                dmc.Text(
                                    "TSV/CSV format (rows matched by ship header)",
                                    size="sm",
                                    c="dimmed",
                                ),
                            ),
                            dmc.AccordionPanel(
                                dmc.Stack(
                                    [
                                        dmc.Text(
                                            "Column headers must match exactly. "
                                            "Only ship_header is required in the file; "
                                            "include any of the other columns below to fill the table. "
                                            "Rows are matched to FASTA headers; omitted ships keep existing table values.",
                                            size="sm",
                                            c="dimmed",
                                        ),
                                        dmc.Text(
                                            "Allowed columns: ship_header, genus, species, strain, "
                                            "assembly_accession, hostchr, shipstart, shipend, strand. "
                                            "Required before submit: genus, species, hostchr, shipstart, shipend. "
                                            "If strand is omitted, it is inferred from start/end order "
                                            "(start > end implies reverse).",
                                            size="sm",
                                            c="dimmed",
                                        ),
                                        dmc.Code(
                                            "ship_header\tgenus\tspecies\tstrain\tassembly_accession\thostchr\tshipstart\tshipend\n"
                                            "seq_1\tAlternaria\talternata\t\tGCA_000001305.1\tchr1\t1000\t5000",
                                            block=True,
                                        ),
                                    ],
                                    gap="xs",
                                ),
                            ),
                        ],
                    ),
                ],
            ),
            dmc.Text(
                "Select cells in the table and paste from a spreadsheet, or upload a TSV/CSV file. "
                "The ship header column is read-only.",
                size="xs",
                c="dimmed",
            ),
        ],
        gap="xs",
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

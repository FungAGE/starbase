import warnings
import logging

from dash import html
import pandas as pd
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag

from src.config.cache import cache
from src.database.sql_manager import fetch_paper_data

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)


def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


def table_loading_alert():
    return html.Div(
        dmc.Alert(
            title="No Data Available",
            children="Waiting for data to load...",
            color="blue",
            variant="light",
            icon=[DashIconify(icon="line-md:loading-loop")],
        ),
        style={"padding": "20px"},
    )


def table_no_results_alert():
    return html.Div(
        dmc.Alert(
            title="No Results Found",
            children="The query returned no results.",
            color="yellow",
            variant="light",
            icon=[DashIconify(icon="clarity:empty-line")],
        ),
        style={"padding": "20px"},
    )


def table_error(e):
    return html.Div(
        dmc.Alert(
            title="Error",
            children=f"Failed to create table: {str(e)}",
            color="red",
            variant="filled",
        ),
        style={"padding": "20px"},
    )


# Common AG Grid configuration
DEFAULT_GRID_OPTIONS = {
    "pagination": True,
    "paginationAutoPageSize": False,
    "domLayout": "autoHeight",
    "tooltipShowDelay": 0,
    "tooltipHideDelay": 1000,
    "enableCellTextSelection": True,
    "ensureDomOrder": True,
    "suppressRowClickSelection": True,
    "rowMultiSelectWithClick": False,
    "onFirstDataRendered": "function(params) { if(params.api) { params.api.sizeColumnsToFit(); }}",
    "rowHeight": 48,
    "headerHeight": 48,
    "suppressRowHoverHighlight": False,
    "suppressHorizontalScroll": False,
    "suppressPropertyNamesCheck": True,
    "suppressReactUi": False,
    "suppressLoadingOverlay": True,
    "suppressNoRowsOverlay": True,
}

DEFAULT_COL_DEF = {
    "resizable": True,
    "sortable": True,
    "filter": True,
    "minWidth": 100,
}


def create_ag_grid(df, id, columns=None, select_rows=False, pg_sz=10):
    """
    Creates an AG Grid component with consistent styling.

    Args:
        df (pd.DataFrame or list): Data to display
        id (str): Unique identifier for the grid
        columns (list): Column definitions
        select_rows (bool): Enable row selection
        pg_sz (int): Number of rows per page
    """
    # Convert DataFrame to row data
    if isinstance(df, pd.DataFrame):
        row_data = df.to_dict("records")
    else:
        row_data = df  # Assume it's already in records format

    # Default column definitions if none provided
    if not columns:
        grid_columns = [{"field": col} for col in df.columns]
    else:
        grid_columns = columns

    # Default column properties
    defaultColDef = {
        "resizable": True,
        "sortable": True,
        "filter": True,
        "minWidth": 100,
    }

    # Simplified getRowId function with proper JavaScript syntax
    get_row_id = "function getRowId(params) { return params.data.accession_tag ? params.data.accession_tag + '_' + Date.now() + '_' + Math.random().toString(36).substr(2, 9) : Date.now() + '_' + Math.random().toString(36).substr(2, 9); }"

    # Create grid component with updated configuration
    grid = dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=row_data,
        defaultColDef=defaultColDef,
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": pg_sz,
            "rowSelection": "multiple" if select_rows else None,
            "domLayout": "autoHeight",
            "tooltipShowDelay": 0,
            "tooltipHideDelay": 1000,
            "enableCellTextSelection": True,
            "ensureDomOrder": True,
            "suppressRowClickSelection": True,
            "rowMultiSelectWithClick": False,
            "onFirstDataRendered": "function(params) { if(params.api) { params.api.sizeColumnsToFit(); }}",
            "rowHeight": 48,
            "headerHeight": 48,
            "suppressRowHoverHighlight": False,
            "suppressHorizontalScroll": False,
            "suppressPropertyNamesCheck": True,
            "suppressReactUi": False,
            "suppressLoadingOverlay": True,
            "suppressNoRowsOverlay": True,
        },
        className="ag-theme-alpine",
        style={"width": "100%", "height": "100%"},
        getRowId=get_row_id,
        persistence=True,
        persistence_type="memory",
    )

    logger.info(f"Successfully created grid {id}")
    return grid


def make_ship_table(df, id, columns=None, select_rows=False, pg_sz=None):
    """
    Specific table constructor for ship data with accession tag handling.

    Args:
        df (pd.DataFrame): Ship data to display
        id (str): Unique identifier for the table
        columns (list): Column definitions
        select_rows (bool): Enable row selection
        pg_sz (int): Number of rows per page
    """
    # Handle empty or None DataFrame
    if df is None or (isinstance(df, pd.DataFrame) and df.empty):
        if columns:
            df = pd.DataFrame(columns=[col["field"] for col in columns])
        else:
            df = pd.DataFrame()

    # Create column definitions
    if columns:
        grid_columns = []
        for col in columns:
            col_def = {
                "field": col["field"],
                "headerName": col["name"]
                if "name" in col
                else col["field"].replace("_", " ").title(),
                "flex": 1,
            }

            # Add special styling for accession_tag
            if col["field"] == "accession_tag":
                col_def.update(
                    {
                        "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
                        "cellClass": "clickable-cell",
                    }
                )

            grid_columns.append(col_def)
    else:
        grid_columns = [{"field": col} for col in df.columns]

    grid_options = {
        **DEFAULT_GRID_OPTIONS,
        "rowSelection": "multiple" if select_rows else None,
        "paginationPageSize": pg_sz or 10,
    }

    return dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=df.to_dict("records") if isinstance(df, pd.DataFrame) else df,
        defaultColDef=DEFAULT_COL_DEF,
        dashGridOptions=grid_options,
        className="ag-theme-alpine",
        style={"width": "100%", "height": "100%"},
        persistence=True,
        persistence_type="memory",
    )


def make_pgv_table(df, id, columns=None, select_rows=False, pg_sz=None):
    """
    Specific table constructor for ship data using AG Grid.
    """
    # Handle empty or None DataFrame
    if df is None or (isinstance(df, pd.DataFrame) and df.empty):
        if columns:
            # Handle both "field" and "id" column formats
            df = pd.DataFrame(
                columns=[col.get("field") or col.get("id") for col in columns]
            )
        else:
            df = pd.DataFrame()

    # Convert column format to AG Grid format
    if columns:
        grid_columns = []
        for col in columns:
            col_def = {
                "field": col.get("id") or col.get("field"),
                "headerName": col.get("name")
                or col.get("headerName")
                or col.get("field", "").replace("_", " ").title(),
                "flex": 1,
            }

            # Add checkboxes for the first column if selection is enabled
            if select_rows and col == columns[0]:
                col_def.update(
                    {
                        "checkboxSelection": True,
                        "headerCheckboxSelection": True,
                    }
                )

            # Add special styling for accession_tag
            if col.get("id") == "accession_tag" or col.get("field") == "accession_tag":
                col_def.update(
                    {
                        "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
                        "cellClass": "clickable-cell",
                    }
                )

            grid_columns.append(col_def)
    else:
        grid_columns = [
            {"field": col, "headerName": col.replace("_", " ").title()}
            for col in df.columns
        ]

    grid_options = {
        **DEFAULT_GRID_OPTIONS,
        "rowSelection": "multiple" if select_rows else None,
        "suppressRowClickSelection": True,  # Don't select rows when clicking cells
        "rowMultiSelectWithClick": False,  # Don't allow multi-select with clicks
        "paginationPageSize": pg_sz or 10,
        "sortable": True,
        "sort": [{"colId": "familyName", "sort": "asc"}],
    }

    return dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=df.to_dict("records") if isinstance(df, pd.DataFrame) else df,
        defaultColDef=DEFAULT_COL_DEF,
        dashGridOptions=grid_options,
        className="ag-theme-alpine",
        style={"width": "100%", "height": "100%"},
        rowSelection="multiple"
        if select_rows
        else None,  # Add this property at component level
        persistence=True,
        persistence_type="memory",
    )


def make_paper_table():
    """Table for displaying paper data."""
    df = cache.get("paper_data")
    if df is None:
        df = fetch_paper_data()
    if df is None:
        df = pd.DataFrame(columns=["Title", "PublicationYear", "Author", "DOI"])
    elif not df.empty:
        df_summary = (
            df.groupby("Title")
            .agg(
                {
                    "familyName": lambda x: ", ".join(sorted(filter(None, x.unique()))),
                    "Author": "first",
                    "PublicationYear": "first",
                    "DOI": "first",
                    "Url": "first",
                }
            )
            .reset_index()
        )
        df = df_summary.sort_values(by="PublicationYear", ascending=False)

        df["DOI"] = df["DOI"].apply(
            lambda x: f"[{x}](https://doi.org/{x})" if pd.notnull(x) else ""
        )

    columns = [
        {
            "field": "Title",
            "headerName": "Title",
            "flex": 2,
            "tooltipField": "Title",
            "wrapText": True,
        },
        {
            "field": "PublicationYear",
            "headerName": "Publication Year",
            "flex": 1,
        },
        {
            "field": "Author",
            "headerName": "Authors",
            "flex": 1.5,
            "tooltipField": "Author",
            "wrapText": True,
        },
        {"field": "DOI", "headerName": "DOI", "flex": 1, "cellRenderer": "markdown"},
    ]

    if isinstance(df, pd.DataFrame):
        row_data = df.fillna("").to_dict("records")
    else:
        row_data = []

    return html.Div(
        dag.AgGrid(
            id="papers-table",
            columnDefs=columns,
            rowData=row_data,
            defaultColDef=DEFAULT_COL_DEF,
            dashGridOptions=DEFAULT_GRID_OPTIONS,
            className="ag-theme-alpine",
            style={"width": "100%", "height": "100%"},
            persistence=True,
            persistence_type="memory",
        ),
        style={"width": "100%"},
    )


def make_dl_table(df, id, table_columns):
    """Table for displaying download data."""
    # Ensure we have a valid data structure
    if df is None or (isinstance(df, list) and len(df) == 0):
        row_data = []
    elif isinstance(df, pd.DataFrame):
        row_data = df.to_dict("records")
    else:
        row_data = df  # Assume it's already in records format

    columns = [
        {
            "field": col["id"],
            "headerName": col["name"],
            "flex": 1,
            "checkboxSelection": col["id"]
            == "accession_tag",  # Add checkbox to first column
            "headerCheckboxSelection": col["id"]
            == "accession_tag",  # Add header checkbox
            "headerCheckboxSelectionFilteredOnly": col["id"]
            == "accession_tag",  # Only select filtered rows
            **(
                {"cellStyle": {"cursor": "pointer", "color": "#1976d2"}}
                if col["id"] == "accession_tag"
                else {}
            ),
            **({"sort": "asc", "sortIndex": 0} if col["id"] == "accession_tag" else {}),
        }
        for col in table_columns
    ]

    grid_options = {
        **DEFAULT_GRID_OPTIONS,
        "rowSelection": "multiple",
        "suppressRowClickSelection": True,  # Prevent row selection on click
        "rowMultiSelectWithClick": False,  # Require checkbox for selection
        "suppressRowDeselection": True,  # Maintain selection when clicking elsewhere
        "paginationPageSize": 25,
        # Add these options for filtered selection
        "suppressHeaderCheckboxSelection": False,  # Enable header checkbox
        "headerCheckboxSelectionFilteredOnly": True,  # Only select filtered rows
        "headerCheckboxSelection": True,  # Enable header checkbox selection
    }

    return html.Div(
        dag.AgGrid(
            id=id,
            columnDefs=columns,
            rowData=row_data,
            defaultColDef=DEFAULT_COL_DEF,
            dashGridOptions=grid_options,
            className="ag-theme-alpine",
            style={"width": "100%", "height": "100%"},
            persistence=True,
            persistence_type="memory",
        )
    )


def make_wiki_table(n_ships, max_size, min_size):
    """Create a summary table for a Starship family."""
    data = [
        {
            "Metric": "Total Number of Starships",
            "Value": f"{n_ships:,.0f}",
        },
        {
            "Metric": "Maximum Size (bp)",
            "Value": f"{max_size:,.0f}",
        },
        {
            "Metric": "Minimum Size (bp)",
            "Value": f"{min_size:,.0f}",
        },
    ]

    columns = [
        {"field": "Metric", "headerName": "Metric", "flex": 2},
        {"field": "Value", "headerName": "Value", "flex": 1},
    ]

    return dag.AgGrid(
        id="wiki-summary-table",
        columnDefs=columns,
        rowData=data,
        defaultColDef=DEFAULT_COL_DEF,
        dashGridOptions={**DEFAULT_GRID_OPTIONS, "paginationPageSize": 10},
        className="ag-theme-alpine",
        style={"width": "100%", "height": "100%"},
        persistence=True,
        persistence_type="memory",
    )

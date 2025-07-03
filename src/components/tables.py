from dash import html
import pandas as pd
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag

from src.config.logging import get_logger
from src.config.cache import cache
from src.database.sql_manager import fetch_paper_data

logger = get_logger(__name__)


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

    # If we have accession_display column, use it for display but keep accession_tag for functionality
    if isinstance(df, pd.DataFrame) and "accession_display" in df.columns:
        # Create a display copy where accession_tag is replaced with accession_display
        display_df = df.copy()
        display_df["accession_tag"] = display_df["accession_display"]
    else:
        display_df = df

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
        grid_columns = [
            {"field": col} for col in display_df.columns if col != "accession_display"
        ]

    # Simplified stable grid options
    grid_options = {
        "pagination": True,
        "paginationPageSize": pg_sz or 10,
        "suppressPropertyNamesCheck": True,
        "rowHeight": 48,
        "headerHeight": 48,
    }

    # Only add selection-related options if selection is enabled
    if select_rows:
        grid_options.update(
            {
                "rowSelection": "multiple",
                "suppressRowClickSelection": True,  # Only select via checkbox
            }
        )

    return dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=display_df.to_dict("records")
        if isinstance(display_df, pd.DataFrame)
        else display_df,
        defaultColDef={
            "resizable": True,
            "sortable": True,
            "filter": True,
            "minWidth": 100,
        },
        dashGridOptions=grid_options,
        className="ag-theme-alpine",
        style={
            "width": "100%",
            "height": "500px",  # Fixed height for stability
            "border": "1px solid #dee2e6",
            "borderRadius": "4px",
        },
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

    # If we have accession_display column, use it for display but keep accession_tag for functionality
    if isinstance(df, pd.DataFrame) and "accession_display" in df.columns:
        # Create a display copy where accession_tag is replaced with accession_display
        display_df = df.copy()
        display_df["accession_tag"] = display_df["accession_display"]
    else:
        display_df = df

    # Convert column format to AG Grid format
    if columns:
        grid_columns = []
        for i, col in enumerate(columns):
            col_def = {
                "field": col.get("id") or col.get("field"),
                "headerName": col.get("name")
                or col.get("headerName")
                or col.get("field", "").replace("_", " ").title(),
                "flex": 1,
                "sortable": True,
                "filter": True,
            }

            # Add checkboxes for the first column if selection is enabled
            if select_rows and i == 0:
                col_def.update(
                    {
                        "checkboxSelection": True,
                        "headerCheckboxSelection": True,
                        "width": 180,
                        "minWidth": 150,
                        "flex": 0,
                    }
                )

            # Add special styling for accession_tag
            if col.get("id") == "accession_tag" or col.get("field") == "accession_tag":
                col_def.update(
                    {
                        "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
                    }
                )

            grid_columns.append(col_def)
    else:
        grid_columns = [
            {
                "field": col,
                "headerName": col.replace("_", " ").title(),
                "sortable": True,
            }
            for col in display_df.columns
            if col != "accession_display"
        ]

    # Simplified grid options
    grid_options = {
        "pagination": True,
        "paginationPageSize": pg_sz or 10,
        "suppressPropertyNamesCheck": True,
        "rowHeight": 48,
        "headerHeight": 48,
    }

    # Only add selection-related options if selection is enabled
    if select_rows:
        grid_options.update(
            {
                "rowSelection": "multiple",
                "suppressRowClickSelection": True,
            }
        )

    return dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=display_df.to_dict("records")
        if isinstance(display_df, pd.DataFrame)
        else display_df,
        defaultColDef={
            "resizable": True,
            "minWidth": 100,
        },
        dashGridOptions=grid_options,
        className="ag-theme-alpine",
        style={
            "width": "100%",
            "height": "450px",  # Increased height for better usability
            "border": "1px solid #dee2e6",
            "borderRadius": "4px",
        },
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
            defaultColDef={
                "resizable": True,
                "minWidth": 100,
            },
            dashGridOptions={
                "pagination": True,
                "paginationPageSize": 10,
                "suppressPropertyNamesCheck": True,
                "rowHeight": 48,
                "headerHeight": 48,
            },
            className="ag-theme-alpine",
            style={
                "width": "100%",
                "height": "500px",
                "border": "1px solid #dee2e6",
                "borderRadius": "4px",
            },
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
        # If we have accession_display column, use it for display but keep accession_tag for functionality
        if "accession_display" in df.columns:
            display_df = df.copy()
            display_df["accession_tag"] = display_df["accession_display"]
            row_data = display_df.to_dict("records")
        else:
            row_data = df.to_dict("records")
    else:
        row_data = df  # Assume it's already in records format

    columns = []
    for i, col in enumerate(table_columns):
        col_def = {
            "field": col["id"],
            "headerName": col["name"],
            "flex": 1,
            "sortable": True,
            "filter": True,
        }

        # Add checkbox and special styling to accession columns
        if col["id"] in ["accession_tag", "accession_display"]:
            col_def.update(
                {
                    "checkboxSelection": True,
                    "headerCheckboxSelection": True,
                    "cellStyle": {"cursor": "pointer", "color": "#1976d2"},
                    "width": 180,
                    "minWidth": 150,
                    "flex": 0,
                }
            )

        columns.append(col_def)

    # Simplified stable grid options
    grid_options = {
        "pagination": True,
        "paginationPageSize": 25,
        "suppressPropertyNamesCheck": True,
        "rowHeight": 48,
        "headerHeight": 48,
        "rowSelection": "multiple",
        "suppressRowClickSelection": True,  # Only select via checkbox
    }

    return html.Div(
        dag.AgGrid(
            id=id,
            columnDefs=columns,
            rowData=row_data,
            defaultColDef={
                "resizable": True,
                "minWidth": 100,
            },
            dashGridOptions=grid_options,
            className="ag-theme-alpine",
            style={
                "width": "100%",
                "height": "600px",  # Fixed height for better stability
                "border": "1px solid #dee2e6",
                "borderRadius": "4px",
            },
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
        defaultColDef={
            "resizable": True,
            "minWidth": 100,
        },
        dashGridOptions={
            "pagination": False,
            "suppressPropertyNamesCheck": True,
            "rowHeight": 48,
            "headerHeight": 48,
        },
        className="ag-theme-alpine",
        style={
            "width": "100%",
            "height": "192px",
            "border": "1px solid #dee2e6",
            "borderRadius": "4px",
        },
        persistence=True,
        persistence_type="memory",
    )

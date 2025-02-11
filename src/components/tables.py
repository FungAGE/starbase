import warnings

warnings.filterwarnings("ignore")

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag
import dash_core_components as dcc

from src.config.cache import cache
from src.database.sql_manager import fetch_paper_data


def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


# Function to convert URL string to HTML link
def url_to_link(url, label):
    return f"[{label}]({url})"


def create_ag_grid(df, id, columns=None, select_rows=False, pg_sz=10):
    """
    Generic AG Grid constructor for tables.
    
    Args:
        df (pd.DataFrame or list): Data to display
        id (str): Unique identifier for the table
        columns (list): Column definitions
        select_rows (bool): Enable row selection
        pg_sz (int): Number of rows per page
    """
    # Convert input data to proper format and handle empty/null cases
    if df is None:
        row_data = []
    elif isinstance(df, pd.DataFrame):
        row_data = df.fillna('').to_dict("records")
    elif isinstance(df, list):
        row_data = df
    else:
        row_data = []
        
    # Ensure we have valid columns even with empty data
    if columns is None and not row_data:
        grid_columns = []
    else:
        defaultColDef = {
            "resizable": True,
            "sortable": True,
            "filter": True,
            "minWidth": 100,
            "flex": 1
        }
        
        # If no columns provided, create them from data
        if columns is None:
            grid_columns = []
            if select_rows:
                grid_columns.append({
                    "headerCheckboxSelection": True,
                    "checkboxSelection": True,
                    "width": 50,
                    "pinned": "left",
                    "lockPosition": True,
                    "suppressMenu": True,
                    "headerName": "",
                    "flex": 0
                })
                
            # Get column names from DataFrame or first dict in list
            if isinstance(df, pd.DataFrame):
                col_names = df.columns
            elif row_data and isinstance(row_data[0], dict):
                col_names = row_data[0].keys()
            else:
                col_names = []
                
            grid_columns.extend([
                {
                    "field": col,
                    "headerName": col.replace("_", " ").title(),
                    "flex": 1
                }
                for col in col_names
            ])
        else:
            grid_columns = columns
            if select_rows:
                grid_columns.insert(0, {
                    "headerCheckboxSelection": True,
                    "checkboxSelection": True,
                    "width": 50,
                    "pinned": "left",
                    "lockPosition": True,
                    "suppressMenu": True,
                    "headerName": "",
                    "flex": 0
                })
    
    grid = dag.AgGrid(
        id=id,
        columnDefs=grid_columns,
        rowData=row_data,
        defaultColDef=defaultColDef,
        dashGridOptions={
            "pagination": True,
            "paginationPageSize": pg_sz,
            "rowSelection": "multiple" if select_rows else None,
            "domLayout": 'autoHeight',
        },
        className="ag-theme-alpine",
        style={"width": "100%"}
    )
    return html.Div([
        html.Div(id=f"{id}-click-data", style={"display": "none"}),
        grid,
        dcc.Store(id=f"{id}-cell-clicked")
    ])

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
    if columns is not None:
        # Add clickable styling for accession_tag columns
        grid_columns = [
            {
                "field": col.get("field"),
                "headerName": col.get("headerName"),
                "flex": col.get("flex", 1),
                "cellStyle": {"cursor": "pointer", "color": "#1976d2"} 
                    if col.get("field") == "accession_tag" else None
            }
            for col in columns
        ]
    else:
        grid_columns = None
        
    return create_ag_grid(
        df=df,
        id=id,
        columns=grid_columns,
        select_rows=select_rows,
        pg_sz=pg_sz or 10
    )

def make_paper_table():
    """Table for displaying paper data."""
    df = cache.get("paper_data")
    if df is None:
        df = fetch_paper_data()
    if df is not None:
        df_summary = (
            df.groupby("Title")
            .agg({
                "familyName": lambda x: ", ".join(sorted(filter(None, x.unique()))),
                "Author": "first",
                "PublicationYear": "first",
                "DOI": "first",
                "Url": "first",
            })
            .reset_index()
        )

        columns = [
            {"field": "Title", "headerName": "Title", "flex": 2},
            {"field": "PublicationYear", "headerName": "Publication Year", "flex": 1},
            {"field": "Author", "headerName": "Authors", "flex": 1.5},
            {"field": "DOI", "headerName": "DOI", "flex": 1, "cellRenderer": "markdown"},
        ]
        
        return create_ag_grid(df_summary.sort_values(by="PublicationYear", ascending=False), 
                            "papers-table", 
                            columns=columns)
    return html.Div("No data available")

def make_ship_blast_table(ship_blast_results, id, df_columns):
    """Table for displaying BLAST results."""
    columns = [
        {
            "field": col["id"],
            "headerName": col["name"],
            "flex": 1,
            "cellStyle": {"cursor": "pointer", "color": "#1976d2"} 
                if col["id"] == "accession_tag" else None
        }
        for col in df_columns
    ]
    
    return create_ag_grid(ship_blast_results, 
                         id, 
                         columns=columns, 
                         select_rows=True)

def make_dl_table(df, id, table_columns):
    """Table for displaying download data."""
    # Ensure we have valid data to work with
    if df is None or (isinstance(df, list) and len(df) == 0):
        df = pd.DataFrame(columns=[col["id"] for col in table_columns])
    
    columns = [
        {
            "field": col["id"],
            "headerName": col["name"],
            "flex": 1
        }
        for col in table_columns
    ]
        
    return create_ag_grid(df, 
                         id, 
                         columns=columns, 
                         select_rows=True)

def make_wiki_table(n_ships, max_size, min_size):
    """Create a summary table for a Starship family."""
    
    data = [
        {
            "metric": "Total Number of Starships",
            "value": f"{n_ships:,.0f}",
        },
        {
            "metric": "Maximum Size (bp)",
            "value": f"{max_size:,.0f}",
        },
        {
            "metric": "Minimum Size (bp)",
            "value": f"{min_size:,.0f}",
        },
    ]

    return html.Div(
        dash_table.DataTable(
            data=data,
            columns=[
                {"name": col, "id": col} for col in ["metric", "value"]
            ],
            style_table={"overflowX": "auto", "maxWidth": "500px"},
            style_cell={"textAlign": "left"},
            style_header={"fontWeight": "bold"},
        )
    )
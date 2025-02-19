import warnings
warnings.filterwarnings("ignore")

import logging

from dash import dash_table, html
import dash_bootstrap_components as dbc
import pandas as pd
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_ag_grid as dag
import dash_core_components as dcc

from src.config.cache import cache
from src.database.sql_manager import fetch_paper_data

logger = logging.getLogger(__name__)
def truncate_string(s, length=40):
    return s if len(s) <= length else s[:length] + "..."


def create_ag_grid(df, id, columns=None, select_rows=False, pg_sz=10):
    """
    Creates an AG Grid component with error handling and consistent styling.
    
    Args:
        df (pd.DataFrame or list): Data to display
        id (str): Unique identifier for the grid
        columns (list): Column definitions
        select_rows (bool): Enable row selection
        pg_sz (int): Number of rows per page
    """
    try:
        # Handle different input types
        if df is None:
            row_data = []
        elif isinstance(df, pd.DataFrame):
            if df.empty:
                row_data = []
            else:
                if 'familyName' in df.columns:
                    df['familyName'] = df['familyName'].fillna('Unclassified')
                    df.loc[df['familyName'].str.strip() == '', 'familyName'] = 'Unclassified'
                row_data = df.fillna('').to_dict("records")
        elif isinstance(df, list):
            row_data = df
        else:
            raise ValueError(f"Unexpected data type for df: {type(df)}")

        # Generate columns if not provided
        if columns is None:
            if not row_data:
                grid_columns = []
            else:
                if isinstance(df, pd.DataFrame):
                    col_names = df.columns
                elif row_data and isinstance(row_data[0], dict):
                    col_names = row_data[0].keys()
                else:
                    col_names = []
                
                grid_columns = []
                for col in col_names:
                    column_def = {
                        "field": col,
                        "headerName": col.replace("_", " ").title(),
                        "flex": 1
                    }
                    
                    # Add styling for accession_tag
                    if col == "accession_tag":
                        column_def["cellStyle"] = {"cursor": "pointer", "color": "#1976d2"}
                        if select_rows:  # Add checkbox to accession_tag column
                            column_def["checkboxSelection"] = True
                            column_def["headerCheckboxSelection"] = True
                            column_def["pinned"] = "left"
                    
                    grid_columns.append(column_def)
        else:
            grid_columns = columns

        # Create grid component with essential configuration
        grid = dag.AgGrid(
            id=id,
            columnDefs=grid_columns,
            rowData=row_data,
            rowSelection="multiple" if select_rows else None,
            pagination=True,
            paginationPageSize=pg_sz,
            getRowId="params.data.accession_tag",
            persistence=True,
            persistence_type="session",
            className="ag-theme-alpine",
            style={"width": "100%"},
            dashGridOptions={
                "domLayout": 'autoHeight',
                "enableCellTextSelection": True,
                "ensureDomOrder": True,
                "onGridReady": """
                    function(params) {
                        params.api.sizeColumnsToFit();
                        window.gridApi = params.api;
                    }
                """,
                "onFirstDataRendered": """
                    function(params) {
                        params.api.sizeColumnsToFit();
                    }
                """,
                "rowHeight": 48,
                "headerHeight": 48,
                "rowMultiSelectWithClick": True if select_rows else False,
                "suppressRowClickSelection": False,
            },
            defaultColDef={
                "resizable": True,
                "sortable": True,
                "filter": True,
                "minWidth": 100,
            }
        )
        
        return grid
        
    except Exception as e:
        return html.Div(
            dmc.Alert(
                title="Error",
                children=f"Failed to create table: {str(e)}",
                color="red",
                variant="filled"
            ),
            style={"padding": "20px"}
        )
    

def make_ship_table(df, id, columns=None, select_rows=False, pg_sz=None):
    """
    Specific table constructor for ship data using DataTable.
    """
    if df is None or (isinstance(df, pd.DataFrame) and df.empty):
        if columns:
            # Handle both "field" and "id" column formats
            df = pd.DataFrame(columns=[col.get("field") or col.get("id") for col in columns])
        else:
            df = pd.DataFrame()
    
    # Convert AG Grid column format to DataTable format
    if columns:
        data_columns = [
            {
                "name": col.get("name") or col.get("headerName") or col.get("field", "").replace("_", " ").title(),
                "id": col.get("id") or col.get("field"),
                "selectable": True,
            }
            for col in columns
        ]
    else:
        data_columns = [
            {"name": col.replace("_", " ").title(), "id": col}
            for col in df.columns
        ]

    return dash_table.DataTable(
        id=id,
        columns=data_columns,
        data=df.to_dict('records'),
        page_size=pg_sz or 10,
        page_current=0,
        page_action='native',
        sort_action='native',
        sort_mode='multi',
        sort_by=[{'column_id': 'familyName', 'direction': 'asc'}],
        row_selectable='multi' if select_rows else None,
        selected_rows=[],
        style_table={'overflowX': 'auto'},
        style_cell={
            'padding': '10px',
            'textAlign': 'left'
        },
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
        style_data_conditional=[
            {
                'if': {'column_id': 'accession_tag'},
                'color': '#1976d2',
                'cursor': 'pointer'
            }
        ]
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
            .agg({
                "familyName": lambda x: ", ".join(sorted(filter(None, x.unique()))),
                "Author": "first",
                "PublicationYear": "first",
                "DOI": "first",
                "Url": "first",
            })
            .reset_index()
        )
        df = df_summary.sort_values(by="PublicationYear", ascending=False)
        
        df['DOI'] = df['DOI'].apply(
            lambda x: f'[{x}](https://doi.org/{x})'
            if pd.notnull(x) else ''
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
        {
            "field": "DOI", 
            "headerName": "DOI", 
            "flex": 1,
            "cellRenderer": "markdown"
        },
    ]
    
    if isinstance(df, pd.DataFrame):
        row_data = df.fillna('').to_dict("records")
    else:
        row_data = []
        
    return html.Div(
        create_ag_grid(
            df=row_data,
            id="papers-table", 
            columns=columns
        ),
        style={"width": "100%"}
    )

def make_dl_table(df, id, table_columns):
    """Table for displaying download data."""
    if df is None or (isinstance(df, list) and len(df) == 0):
        df = pd.DataFrame(columns=[col["id"] for col in table_columns])
    
    columns = [
        {
            "field": col["id"],
            "headerName": col["name"],
            "flex": 1,
            **({"cellStyle": {"cursor": "pointer", "color": "#1976d2"}}
               if col["id"] == "accession_tag" else {}),
        }
        for col in table_columns
    ]
        
    return create_ag_grid(
        df=df, 
        id=id, 
        columns=columns, 
        select_rows=True,
        pg_sz=25
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

    return create_ag_grid(
        df=data,
        id="wiki-summary-table",
        columns=columns,
        select_rows=False,
        pg_sz=10
    )
from dash import Output, Input, State, callback, callback_context
from dash.exceptions import PreventUpdate
import dash_mantine_components as dmc
import functools
import traceback

from src.config.logging import get_logger

logger = get_logger(__name__)


def handle_callback_error(callback_func):
    """
    Decorator to handle callback errors gracefully

    Args:
        callback_func: The callback function to wrap

    Returns:
        Wrapped function with error handling
    """

    @functools.wraps(callback_func)
    def wrapper(*args, **kwargs):
        try:
            return callback_func(*args, **kwargs)
        except PreventUpdate:
            raise
        except Exception as e:
            # Log the error with full traceback
            logger.error(f"Callback error in {callback_func.__name__}: {str(e)}")
            logger.error(f"Inputs: args={args}, kwargs={kwargs}")
            logger.error(traceback.format_exc())

            # Return detailed error alert
            return dmc.Alert(
                title="Error Loading Data",
                children=[
                    "We encountered a problem processing your request.",
                    dmc.Space(h=10),
                    dmc.Text(f"Error: {str(e)}", size="sm"),
                    dmc.Space(h=10),
                    dmc.Code(str(traceback.format_exc()), block=True),
                ],
                color="var(--mantine-color-red-6)",
                variant="filled",
            )

    return wrapper


def _extract_accession(raw_value, prefix):
    """Pull the accession token (e.g. SSB0001) out of a cell value that may be
    bracket-wrapped or path-like (e.g. '[SSB0001]', 'some/path/SSB0001')."""
    parts = [p.strip() for p in str(raw_value).strip("[]").split("/") if p.strip()]
    match = next((p for p in parts if p.upper().startswith(prefix)), None)
    if match:
        return match
    return parts[-1] if parts else None


def create_modal_callback(
    table_id,
    modal_id,
    content_id,
    title_id,
    ship_cols=("ship_accession_tag", "ship_accession_display"),
    group_cols=("accession_tag", "accession_display"),
):
    """Open a native dmc.Modal populated from the ship/group accession data
    when a user clicks an accession cell in an AG Grid or Dash DataTable."""
    from src.components.data import (
        create_ship_accession_modal_data,
        create_accession_modal_data,
        render_ship_accession_modal,
        render_group_accession_modal,
    )

    @callback(
        Output(modal_id, "opened"),
        Output(content_id, "children"),
        Output(title_id, "children"),
        [
            Input(table_id, "cellClicked"),  # AG Grid
            Input(table_id, "active_cell"),
        ],  # Dash DataTable
        [
            State(table_id, "derived_virtual_data"),
            State(table_id, "page_current"),
            State(table_id, "page_size"),
        ],
        prevent_initial_call=True,
    )
    def toggle_modal(cell_clicked, active_cell, table_data, page_current, page_size):
        ctx = callback_context
        if not ctx.triggered or not (cell_clicked or active_cell):
            raise PreventUpdate

        triggered_id = ctx.triggered[0]["prop_id"]
        col_id = None
        value = None

        if f"{table_id}.cellClicked" in triggered_id and cell_clicked:
            col_id = cell_clicked.get("colId")
            value = cell_clicked.get("value")
        elif f"{table_id}.active_cell" in triggered_id and active_cell:
            col_id = active_cell.get("column_id")
            row_idx = (page_current or 0) * (page_size or 0) + active_cell.get("row", 0)
            if table_data and row_idx < len(table_data):
                value = table_data[row_idx].get(col_id)

        if col_id is None or not value:
            raise PreventUpdate

        if col_id in ship_cols:
            prefix = "SSB"
            fetch_data, render_content = create_ship_accession_modal_data, render_ship_accession_modal
        elif col_id in group_cols:
            prefix = "SSA"
            fetch_data, render_content = create_accession_modal_data, render_group_accession_modal
        else:
            raise PreventUpdate

        accession = _extract_accession(value, prefix)
        if not accession:
            raise PreventUpdate

        try:
            data = fetch_data(accession)
            return True, render_content(data), data.get("title", accession)
        except Exception as e:
            logger.error(f"Error building accession modal for {accession}: {str(e)}")
            logger.error(traceback.format_exc())
            return (
                True,
                dmc.Alert(
                    f"Error loading details for {accession}.",
                    color="red",
                    variant="light",
                ),
                "Error",
            )

    return toggle_modal

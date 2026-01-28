from dash import html, Output, Input, State, callback, no_update, callback_context
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
                color="red",
                variant="filled",
            )

    return wrapper


def create_modal_callback(table_id, modal_id, content_id, title_id, column_check=None):
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
        try:
            ctx = callback_context
            triggered_id = ctx.triggered[0]["prop_id"]

            if not (cell_clicked or active_cell):
                return False, no_update, no_update

            accession = None
            if f"{table_id}.cellClicked" in triggered_id and cell_clicked:
                # AG Grid format - handle both accession_tag and accession_display
                if cell_clicked["colId"] in [
                    "ship_accession_tag",
                    "ship_accession_display",
                ]:
                    accession = str(cell_clicked["value"])
            elif f"{table_id}.active_cell" in triggered_id and active_cell:
                # Dash DataTable format - handle both accession_tag and accession_display
                if active_cell["column_id"] in [
                    "ship_accession_tag",
                    "ship_accession_display",
                ]:
                    # Calculate the actual row index based on pagination
                    actual_row_idx = (page_current or 0) * page_size + active_cell[
                        "row"
                    ]
                    if table_data and actual_row_idx < len(table_data):
                        accession = str(
                            table_data[actual_row_idx][active_cell["column_id"]]
                        )

            if accession:
                # Clean and standardize the accession tag
                accession = accession.strip("[]").split("/")[-1].strip()
                logger.debug(f"Looking for accession in cache: {accession}")

                # Instead of opening a Dash modal, trigger the universal modal via JavaScript
                return (
                    False,  # Don't open the Dash modal
                    no_update,
                    no_update,
                )

            return False, no_update, no_update
        except Exception as e:
            logger.error(f"Error in toggle_modal: {str(e)}")
            logger.error(traceback.format_exc())
            error_content = html.Div(
                [
                    html.P("Error loading modal content"),
                    html.P(f"Details: {str(e)}"),
                ]
            )
            error_title = dmc.Title("Error", order=3)
            return True, error_content, error_title

    return toggle_modal

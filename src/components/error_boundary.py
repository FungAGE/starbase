from dash import html
import dash_mantine_components as dmc
from functools import wraps
import traceback
import logging

logger = logging.getLogger(__name__)

def error_boundary(children, id=None):
    """Wraps content in a div with error handling"""
    # Create the Paper component without id if none provided
    paper = dmc.Paper(
        children=children,
        p="md",
        withBorder=True,
        style={"position": "relative"}  # Add relative positioning
    )
    
    # Only wrap in a div with id if id is provided
    if id:
        return html.Div(paper, id=id, style={"position": "relative"})  # Add relative positioning
    return paper

def handle_callback_error(func):
    """Decorator to handle callback errors gracefully"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Callback error in {func.__name__}: {str(e)}")
            logger.error(traceback.format_exc())
            return dmc.Alert(
                title="Error",
                children="An error occurred while loading this component. Please try refreshing the page.",
                color="red",
                variant="filled",
                style={
                    "position": "relative",  # Add relative positioning
                    "zIndex": 1,  # Ensure proper stacking
                    "marginTop": "1rem"
                }
            )
    return wrapper 
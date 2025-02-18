from dash import html
import dash_mantine_components as dmc
from functools import wraps
import traceback
import logging

logger = logging.getLogger(__name__)

def error_boundary(children, id=None):
    """Wraps content in a div with error handling"""
    paper = dmc.Paper(
        children=children,
        p="md",
        withBorder=True,
        style={"position": "relative"}
    )
    
    # Only wrap in a div with id if id is provided
    if id:
        return html.Div(paper, id=id, style={"position": "relative"})
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
                children="An error occurred while loading this component. This is likely due to a temporary issue with the server. Please try refreshing the page in a few minutes.",
                color="red",
                variant="filled",
                style={
                    "position": "relative",
                    "zIndex": 1,
                    "marginTop": "1rem"
                }
            )
    return wrapper 
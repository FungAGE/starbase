from dash import html
import dash_mantine_components as dmc
from functools import wraps
import traceback
import logging
import functools
from dash.exceptions import PreventUpdate

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
            
            # Return user-friendly error component
            return dmc.Alert(
                title="An Error Occurred",
                children=[
                    "We encountered a problem processing your request.",
                    dmc.Space(h=10),
                    dmc.Text("The error has been logged and we'll look into it.", size="sm"),
                    dmc.Space(h=10),
                    dmc.Button(
                        "Refresh Page",
                        variant="outline",
                        color="red",
                        onClick="window.location.reload(true)"
                    )
                ],
                color="red",
                variant="filled",
            )
    return wrapper

def create_error_boundary(page_content):
    """
    Wrap page content with error boundary
    
    Args:
        page_content: The page layout to wrap
        
    Returns:
        Error boundary wrapped content
    """
    error_wrapped = html.Div([
        # Wrap the page content directly without LoadingOverlay
        html.Div(
            page_content,
            id="page-content-wrapper"
        ),
        html.Div(id="error-boundary-content")
    ])
    return error_wrapped 
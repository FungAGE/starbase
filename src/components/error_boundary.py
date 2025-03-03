from dash import html
import dash_mantine_components as dmc
from functools import wraps
import traceback
import logging
import functools
from dash.exceptions import PreventUpdate
from dash import callback, Input, Output, callback_context
import dash_core_components as dcc

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
            # Add validation for input args
            for arg in args:
                if arg is None:
                    logger.warning(f"Received None argument in {callback_func.__name__}")
                    error_component = dmc.Alert(
                        title="Invalid Input",
                        children="Please try refreshing the page.",
                        color="yellow",
                        variant="filled"
                    )
                    # Get number of outputs from callback context
                    num_outputs = len(callback_context.outputs_list) if callback_context.outputs_list else 1
                    return tuple([error_component] * num_outputs)

            result = callback_func(*args, **kwargs)
            
            # Validate callback output
            if result is None:
                logger.error(f"Callback {callback_func.__name__} returned None")
                raise ValueError("Callback returned None")
                
            return result

        except PreventUpdate:
            raise
        except Exception as e:
            logger.error(f"Callback error in {callback_func.__name__}: {str(e)}")
            logger.error(f"Inputs: args={args}, kwargs={kwargs}")
            logger.error(traceback.format_exc())
            
            error_component = dmc.Alert(
                title="An Error Occurred",
                children=[
                    dmc.Text(
                        f"Error type: {type(e).__name__}",
                        size="sm",
                        c="red"
                    ),
                    dmc.Space(h=10),
                    dmc.Text(
                        "We encountered a problem processing your request.",
                        size="sm"
                    ),
                    dmc.Space(h=10),
                    dmc.Text(
                        "The error has been logged and we'll look into it.",
                        size="xs",
                        c="gray"
                    ),
                    dmc.Space(h=10),
                    dmc.Group([
                        dmc.Button(
                            "Refresh Page",
                            variant="outline",
                            color="red",
                            id="refresh-page-button"
                        ),
                        dmc.Button(
                            "Clear Cache",
                            variant="subtle",
                            color="gray",
                            id="clear-cache-button"
                        )
                    ])
                ],
                color="red",
                variant="filled",
                style={"maxWidth": "500px", "margin": "auto"}
            )
            
            # Get number of outputs from callback context
            num_outputs = len(callback_context.outputs_list) if callback_context.outputs_list else 1
            return tuple([error_component] * num_outputs)
            
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
        html.Div(
            page_content,
            style={
                "marginTop": "60px",  # Add margin to avoid navbar overlap
                "minHeight": "calc(100vh - 60px)"  # Adjust height to account for navbar
            }
        ),
        dmc.LoadingOverlay(
            id="page-loading-overlay",
            overlayProps={
                "color": "rgba(255, 255, 255, 0.8)",
                "children": page_content,
                "style": {
                    "marginTop": "60px",  # Add margin to avoid navbar overlap
                    "minHeight": "calc(100vh - 60px)"  # Adjust height to account for navbar
            }},
            loaderProps={"variant": "dots", "color": "blue"},
            visible=False,
        ),
        html.Div(id="error-boundary-content"),
        dmc.Notification(
            id="global-error-notification",
            title="Error",
            message="",
            color="red",
            autoClose=5000,
            action="hide",
            style={"position": "fixed", "top": 80, "right": 20}
        )
    ])
    return error_wrapped

# Add a global error handler callback
@callback(
    Output("global-error-notification", "message"),
    Output("global-error-notification", "opened"),
    Input("page-content-wrapper", "children"),
    prevent_initial_call=True
)
@handle_callback_error
def handle_global_errors(content):
    if isinstance(content, dict) and content.get("props", {}).get("title") == "An Error Occurred":
        return "An error occurred while loading the page. Please try refreshing.", True
    return "", False

# Add callbacks for button actions
@callback(
    Output("refresh-page-button", "n_clicks"),
    Input("refresh-page-button", "n_clicks"),
    prevent_initial_call=True
)
def refresh_page(n_clicks):
    if n_clicks:
        return dcc.Location(href=None, refresh=True)
    raise PreventUpdate

@callback(
    Output("clear-cache-button", "n_clicks"),
    Input("clear-cache-button", "n_clicks"),
    prevent_initial_call=True
)
def clear_cache(n_clicks):
    if n_clicks:
        return dcc.Location(href="/", refresh=True)
    raise PreventUpdate 
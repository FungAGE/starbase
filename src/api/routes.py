from flask import Blueprint, jsonify, request
from src.config.limiter import limiter
from src.components.callbacks import (
    create_ship_accession_modal_data,
    create_accession_modal_data,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

# Create Blueprints for different route groups
blast_routes = Blueprint("blast", __name__, url_prefix="/api/blast")
accession_routes = Blueprint("accession", __name__, url_prefix="/api/accession")
error_handlers = Blueprint("errors", __name__)


@blast_routes.route("/blast-submit", methods=["POST"])
@limiter.limit("10 per hour")
def check_blast_limit():
    """Check BLAST submission limit."""
    remote_addr = request.remote_addr
    logger.info(f"BLAST submission from IP: {remote_addr}")
    return jsonify({"allowed": True})


@accession_routes.route("/ship_accession_details/<ship_accession_id>", methods=["GET"])
def get_ship_accession_details(ship_accession_id):
    """Get details for a specific ship accession."""
    try:
        modal_data = create_ship_accession_modal_data(ship_accession_id)

        return jsonify(modal_data)
    except Exception as e:
        logger.error(f"Error in get_ship_accession_details: {str(e)}")
        return jsonify({"error": str(e)}), 500


@accession_routes.route("/accession_details/<accession_id>", methods=["GET"])
def get_accession_details(accession_id):
    """Get details for a specific ship accession."""
    try:
        modal_data = create_accession_modal_data(accession_id)

        return jsonify(modal_data)
    except Exception as e:
        logger.error(f"Error in get_ship_accession_details: {str(e)}")
        return jsonify({"error": str(e)}), 500


def dash_to_html(component):
    """Convert Dash components directly to HTML string."""
    import dash_mantine_components as dmc
    from dash.development.base_component import Component

    # Base cases
    if component is None:
        return ""
    elif isinstance(component, str):
        return component
    elif isinstance(component, (int, float)):
        return str(component)
    elif isinstance(component, list):
        return "".join([dash_to_html(c) for c in component if c is not None])

    # Not a Dash component
    if not isinstance(component, Component):
        return str(component)

    # Get component type
    component_type = (
        component._type if hasattr(component, "_type") else type(component).__name__
    )

    # Ship Details section
    if isinstance(component, dmc.Paper):
        children = component.children if hasattr(component, "children") else []

        # Extract the section title (first child is usually the title)
        section_title = ""
        section_content = ""

        for child in children if isinstance(children, list) else [children]:
            if isinstance(child, dmc.Title):
                section_title = dash_to_html(child)
            else:
                section_content += dash_to_html(child)

        return f"""
        <div class="modal-section" style="margin-bottom:20px; border:1px solid #dee2e6; border-radius:4px; padding:16px; background-color:white;">
            <div class="section-title" style="font-weight:600; font-size:18px; margin-bottom:16px;">{section_title}</div>
            <div class="section-content">{section_content}</div>
        </div>
        """

    # Section titles
    elif isinstance(component, dmc.Title):
        children = component.children if hasattr(component, "children") else ""
        return f"<h3 style='margin:0; font-weight:600;'>{dash_to_html(children)}</h3>"

    # Grid layouts
    elif isinstance(component, dmc.SimpleGrid):
        children = component.children if hasattr(component, "children") else []
        return f"""
        <div style="display:grid; grid-template-columns:repeat(2, 1fr); gap:16px; margin-bottom:16px;">
            {dash_to_html(children)}
        </div>
        """

    # Label-value groups
    elif isinstance(component, dmc.Group):
        children = component.children if hasattr(component, "children") else []

        # Typically a Group contains a label and a value
        items = [
            dash_to_html(child)
            for child in (children if isinstance(children, list) else [children])
        ]

        if len(items) >= 2:
            # First item is usually the label, rest is the value
            label = items[0]
            value = "".join(items[1:])
            return f"""
            <div style="display:flex; justify-content:space-between; margin-bottom:8px;">
                <div>{label}</div>
                <div>{value}</div>
            </div>
            """
        else:
            return "".join(items)

    # Text fields - bold for labels
    elif isinstance(component, dmc.Text):
        children = component.children if hasattr(component, "children") else ""
        fw = component.fw if hasattr(component, "fw") else None

        if fw and fw >= 700:
            return f"<strong>{dash_to_html(children)}</strong>"
        return dash_to_html(children)

    # Badges
    elif isinstance(component, dmc.Badge):
        children = component.children if hasattr(component, "children") else ""
        color = component.color if hasattr(component, "color") else "gray"

        colors = {
            "green": ("#e9faf0", "#099268"),
            "yellow": ("#fff8e1", "#f59f00"),
            "blue": ("#e7f5ff", "#1c7ed6"),
            "gray": ("#f1f3f5", "#495057"),
        }
        bg_color, text_color = colors.get(color, colors["gray"])

        return f"""
        <span style="background-color:{bg_color}; color:{text_color}; padding:2px 8px; 
                     border-radius:4px; font-size:12px; font-weight:500;">
            {dash_to_html(children)}
        </span>
        """

    # Tables
    elif isinstance(component, dmc.Table) or component_type == "Table":
        children = component.children if hasattr(component, "children") else []
        return f"""
        <table style="width:100%; border-collapse:collapse; margin-bottom:16px;">
            {dash_to_html(children)}
        </table>
        """

    # Table rows
    elif component_type == "Tr":
        children = component.children if hasattr(component, "children") else []
        return f"<tr>{dash_to_html(children)}</tr>"

    # Table headers
    elif component_type == "Th":
        children = component.children if hasattr(component, "children") else []
        return f"<th style='text-align:left; padding:8px; border-bottom:1px solid #dee2e6; background-color:#f8f9fa;'>{dash_to_html(children)}</th>"

    # Table cells
    elif component_type == "Td":
        children = component.children if hasattr(component, "children") else []
        return f"<td style='padding:8px; border-bottom:1px solid #dee2e6;'>{dash_to_html(children)}</td>"

    # Table header
    elif component_type == "Thead":
        children = component.children if hasattr(component, "children") else []
        return f"<thead>{dash_to_html(children)}</thead>"

    # Table body
    elif component_type == "Tbody":
        children = component.children if hasattr(component, "children") else []
        return f"<tbody>{dash_to_html(children)}</tbody>"

    # Links
    elif isinstance(component, dmc.Anchor):
        children = component.children if hasattr(component, "children") else ""
        href = component.href if hasattr(component, "href") else "#"
        target = component.target if hasattr(component, "target") else "_blank"

        return f"""<a href="{href}" target="{target}" style="color:#228be6; text-decoration:none;">
                    {dash_to_html(children)}
                   </a>"""

    # Stacks - vertical containers
    elif isinstance(component, dmc.Stack):
        children = component.children if hasattr(component, "children") else []
        return f"""
        <div style="display:flex; flex-direction:column; gap:16px;">
            {dash_to_html(children)}
        </div>
        """

    # Alerts - error/warning messages
    elif isinstance(component, dmc.Alert):
        children = component.children if hasattr(component, "children") else ""
        title = component.title if hasattr(component, "title") else ""
        color = component.color if hasattr(component, "color") else "blue"

        colors = {
            "red": ("#fef2f2", "#dc2626", "#fecaca"),
            "yellow": ("#fffbeb", "#d97706", "#fed7aa"),
            "blue": ("#f0f9ff", "#0284c7", "#bae6fd"),
        }
        bg_color, text_color, border_color = colors.get(color, colors["blue"])

        return f"""
        <div style="background-color:{bg_color}; color:{text_color}; border:1px solid {border_color}; 
                    border-radius:4px; padding:12px 16px; margin-bottom:16px;">
            {f'<div style="font-weight:600; margin-bottom:8px;">{title}</div>' if title else ""}
            <div>{dash_to_html(children)}</div>
        </div>
        """

    # Flex containers
    elif isinstance(component, dmc.Flex):
        children = component.children if hasattr(component, "children") else []
        wrap = component.wrap if hasattr(component, "wrap") else "nowrap"
        gap = component.gap if hasattr(component, "gap") else "8px"

        return f"""
        <div style="display:flex; gap:{gap}; flex-wrap:{wrap};">
            {dash_to_html(children)}
        </div>
        """

    # Default fallback - get children if available
    children = component.children if hasattr(component, "children") else ""
    if (
        not children
        and hasattr(component, "props")
        and component.props
        and "children" in component.props
    ):
        children = component.props["children"]

    return dash_to_html(children)


@error_handlers.app_errorhandler(500)
def handle_500(e):
    """Handle internal server errors."""
    logger.error(f"Internal Server Error: {str(e)}")
    return jsonify(
        {
            "error": "Internal Server Error",
            "message": "The server encountered an error. Please try again later.",
        }
    ), 500


@error_handlers.app_errorhandler(429)
def handle_429(e):
    """Handle too many requests errors."""
    return jsonify(
        {
            "error": "Too Many Requests",
            "message": "Please wait before making more requests.",
        }
    ), 429

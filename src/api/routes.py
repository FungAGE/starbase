from flask import Blueprint, jsonify, request
from src.config.limiter import limiter
from src.components.data import (
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

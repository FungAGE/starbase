from flask import Blueprint, jsonify, request
from sqlalchemy import text
from src.utils.telemetry import log_request
from src.config.cache import cache, cleanup_old_cache, cache_dir
from src.config.database import TelemetrySession
from src.components.callbacks import create_accession_modal
from src.config.limiter import limiter

from src.config.logging import get_logger
logger = get_logger(__name__)

# Create Blueprints for different route groups
blast_routes = Blueprint('blast', __name__)
accession_routes = Blueprint('accession', __name__)
error_handlers = Blueprint('errors', __name__)

@blast_routes.route('/blast-submit', methods=['POST'])
@limiter.limit("10 per hour")
def check_blast_limit():
    """Check BLAST submission limit."""
    remote_addr = request.remote_addr
    logger.info(f"BLAST submission from IP: {remote_addr}")
    return {"allowed": True}


@accession_routes.route('/accession/<accession_id>', methods=['GET'])
def get_accession_details(accession_id):
    """Get details for a specific accession."""
    try:
        content, title = create_accession_modal(accession_id)
        return jsonify({
            'content': content,
            'title': title
        })
    except Exception as e:
        logger.error(f"Error in get_accession_details: {str(e)}")
        return jsonify({'error': str(e)}), 500


@error_handlers.app_errorhandler(500)
def handle_500(e):
    """Handle internal server errors."""
    logger.error(f"Internal Server Error: {str(e)}")
    return jsonify({
        "error": "Internal Server Error",
        "message": "The server encountered an error. Please try again later."
    }), 500


@error_handlers.app_errorhandler(429)
def handle_429(e):
    """Handle too many requests errors."""
    return jsonify({
        "error": "Too Many Requests",
        "message": "Please wait before making more requests."
    }), 429
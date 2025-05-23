"""
Telemetry routes module.
Contains all Flask routes related to telemetry.
"""

from flask import Blueprint, jsonify
from sqlalchemy import text
from src.config.database import TelemetrySession
from src.config.settings import IPSTACK_API_KEY
from src.config.cache import cache
from src.telemetry.utils import log_request, update_ip_locations
from src.config.logging import get_logger

logger = get_logger(__name__)

# Create a Blueprint for telemetry routes
telemetry_routes = Blueprint("telemetry", __name__, url_prefix="/api/telemetry")


@telemetry_routes.route("/health", methods=["GET"])
def health():
    """Check telemetry system health."""
    try:
        session = TelemetrySession()
        test_ip = "127.0.0.1"
        test_endpoint = "/"
        log_request(test_ip, test_endpoint)
        result = session.execute(text("SELECT COUNT(*) FROM request_logs")).scalar()
        return jsonify({"status": "healthy", "record_count": result}), 200
    except Exception as e:
        return jsonify({"status": "unhealthy", "error": str(e)}), 503
    finally:
        session.close()


@telemetry_routes.route("/refresh", methods=["POST"])
def refresh():
    """Endpoint to refresh telemetry data."""
    try:
        update_ip_locations(IPSTACK_API_KEY)
        cache.delete("telemetry_data")
        return jsonify({"status": "success"}), 200
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500

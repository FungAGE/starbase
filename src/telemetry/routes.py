"""
Telemetry routes module.
Contains all Flask routes related to telemetry.
"""

from flask import Blueprint, jsonify
from sqlalchemy import text
from src.config.settings import IPSTACK_API_KEY
from src.config.cache import cache
from src.telemetry.utils import log_request, update_ip_locations
from src.config.logging import get_logger

logger = get_logger(__name__)

# Create a Blueprint for telemetry routes
telemetry_routes = Blueprint("telemetry", __name__, url_prefix="/api/telemetry")


@telemetry_routes.route("/health", methods=["GET"])
def health():
    """Check telemetry system health - verifies DB and task system."""
    from src.database.sql_engine import get_telemetry_session
    from src.config.celery_config import CELERY_AVAILABLE
    
    health_data = {"status": "healthy"}
    status_code = 200
    
    try:
        # Test database connectivity
        with get_telemetry_session() as session:
            count = session.execute(
                text("SELECT COUNT(*) FROM request_logs")
            ).scalar()
            health_data["database"] = "connected"
            health_data["record_count"] = count
    except Exception as e:
        logger.error(f"Database health check failed: {str(e)}")
        health_data["status"] = "unhealthy"
        health_data["database"] = f"error: {str(e)}"
        status_code = 503
    
    # Check task system status
    health_data["task_system"] = "celery" if CELERY_AVAILABLE else "direct"
    
    return jsonify(health_data), status_code

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

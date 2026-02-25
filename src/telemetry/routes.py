"""
Telemetry routes module.
Contains all Flask routes related to telemetry.
"""

from flask import Blueprint, jsonify, request
from functools import wraps
from sqlalchemy import text
from src.config.cache import cache
from src.telemetry.tasks import update_ip_locations_task
from src.database.sql_engine import get_telemetry_session
from src.config.logging import get_logger
from src.telemetry.utils import (
    get_client_ip,
    is_development_ip,
    get_blast_limit_info,
)

logger = get_logger(__name__)

telemetry_routes = Blueprint("telemetry", __name__, url_prefix="/api/telemetry")


@telemetry_routes.route("/health", methods=["GET"])
def health():
    """Check telemetry system health - verifies DB and task system."""
    from src.config.celery_config import CELERY_AVAILABLE

    health_data = {"status": "healthy"}
    status_code = 200

    try:
        # Test database connectivity
        with get_telemetry_session() as session:
            count = session.execute(text("SELECT COUNT(*) FROM request_logs")).scalar()
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
        update_ip_locations_task()
        cache.delete("telemetry_data")
        return jsonify({"status": "success"}), 200
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500


def blast_limit_decorator(f):
    """Decorator to limit BLAST operations per user"""
    from dash.exceptions import PreventUpdate

    @wraps(f)
    def wrapped(*args, **kwargs):
        try:
            client_ip = get_client_ip()

            # Skip rate limiting for development IPs
            if is_development_ip(client_ip):
                return f(*args, **kwargs)

            # Only check limits for BLAST-related endpoints
            if request.path == "/api/blast-submit":
                limit_info = get_blast_limit_info(client_ip)

                if limit_info["remaining"] <= 0:
                    logger.warning(f"BLAST limit exceeded for IP: {client_ip}")
                    raise PreventUpdate("Hourly BLAST limit exceeded")

                # Log the BLAST request
                logger.debug(f"BLAST submission from IP: {client_ip}")
                from src.telemetry.tasks import log_request_task
                from src.config.celery_config import run_task

                run_task(log_request_task, client_ip, "/api/blast-submit")

            return f(*args, **kwargs)

        except PreventUpdate:
            raise
        except Exception as e:
            logger.error(f"Error in blast_limit_decorator: {str(e)}")
            raise

    return wrapped

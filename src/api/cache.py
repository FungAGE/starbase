from flask import Blueprint, jsonify
from src.config.cache import cache, cleanup_old_cache, cache_dir
from src.config.limiter import limiter
import logging

from src.config.logging import get_logger
logger = get_logger(__name__)

# Create a Blueprint for cache routes
cache_routes = Blueprint('cache', __name__)

@cache_routes.route('/status', methods=['GET'])
def cache_status():
    """Check cache status."""
    try:
        return jsonify({
            "status": "success",
            "cache_dir": cache_dir,
            "cache_type": cache.config['CACHE_TYPE']
        }), 200
    except Exception as e:
        logger.error(f"Cache status check failed: {str(e)}")
        return jsonify({"status": "error", "error": str(e)}), 500


@cache_routes.route('/refresh', methods=['POST'])
@limiter.limit("1 per minute")
def refresh_cache():
    """Force refresh of all cached data."""
    try:
        cache.clear()
        cleanup_old_cache()
        return jsonify({
            "status": "success",
            "message": "Cache cleared and old files cleaned up"
        }), 200
    except Exception as e:
        logger.error(f"Cache refresh failed: {str(e)}")
        return jsonify({"status": "error", "error": str(e)}), 500


@cache_routes.route('/cleanup', methods=['POST'])
@limiter.limit("1 per minute")
def force_cache_cleanup():
    """Force cleanup of old cache files."""
    try:
        cleanup_old_cache()
        return jsonify({
            "status": "success",
            "message": "Cache cleanup completed"
        }), 200
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return jsonify({"status": "error", "error": str(e)}), 500
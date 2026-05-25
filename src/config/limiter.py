import os

from flask import request
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from src.config.settings import IS_DEV
from src.config.logging import get_logger

logger = get_logger(__name__)

# Use Redis for rate limiting in production when available (shared across workers)
_redis_url = os.getenv("REDIS_URL") or os.getenv("CELERY_BROKER_URL")
_storage_uri = _redis_url if (_redis_url and "redis" in _redis_url) else "memory://"
if _storage_uri != "memory://":
    logger.info("Using Redis for rate limiting")

limiter = Limiter(
    get_remote_address,
    storage_uri=_storage_uri,
    default_limits=["500 per minute", "10000 per hour"],
    strategy="fixed-window-elastic-expiry",
)


@limiter.request_filter
def exclude_dev_endpoints():
    """Exclude Dash development endpoints from rate limiting in development mode only."""
    # Only bypass rate limits for dev endpoints if IS_DEV
    if IS_DEV and request.path.startswith(("/_dash-", "/_reload-hash")):
        return True
    return False


# TODO: Commenting this out so that rate limits can be tested in development
# Consider using a flag to turn this on and off
# @limiter.request_filter
# def log_rate_limit():
#    """Custom rate limiter filter to bypass limits for development IPs."""
#    remote_addr = get_client_ip()
#    if is_development_ip(remote_addr):
#        return True
#    logger.info(f"Rate limit hit by IP: {remote_addr}")
#    return False

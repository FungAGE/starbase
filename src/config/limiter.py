from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from src.config.logging import get_logger

logger = get_logger(__name__)

# TODO: Consider using a different storage_uri for production
limiter = Limiter(
    get_remote_address,
    storage_uri="memory://",
    default_limits=["500 per minute", "10000 per hour"],
    strategy="fixed-window-elastic-expiry",
)

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

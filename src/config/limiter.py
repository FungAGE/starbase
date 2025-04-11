from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from src.utils.telemetry import log_request, get_client_ip, is_development_ip

from src.config.logging import get_logger
logger = get_logger(__name__)

# TODO: Consider using a different storage_uri for production
limiter = Limiter(
    get_remote_address,
    storage_uri="memory://",
    default_limits=["200 per minute"],
    strategy="fixed-window-elastic-expiry"
)

@limiter.request_filter
def log_rate_limit():
    """Custom rate limiter filter to bypass limits for development IPs."""
    remote_addr = get_client_ip()
    if is_development_ip(remote_addr):
        return True
    logger.info(f"Rate limit hit by IP: {remote_addr}")
    return False

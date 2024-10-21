import logging
from flask_limiter.util import get_remote_address
from app import limiter

# Setup logging
logging.basicConfig(level=logging.INFO)


# Log rate-limited events
@limiter.request_filter
def log_rate_limit():
    remote_addr = get_remote_address()
    logging.info(f"Rate limit hit by IP: {remote_addr}")
    # Alternatively, store the data in a database or monitoring service
    return False  # Don't exclude the request from rate limiting, just log it

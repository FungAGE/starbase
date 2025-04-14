from flask_caching import Cache
import os
import time
import logging

from src.config.settings import DATA_DIR

logger = logging.getLogger(__name__)

cache_dir = os.path.join(DATA_DIR, "cache")
os.makedirs(cache_dir, exist_ok=True)

CACHE_TIMESTAMP = str(int(time.time()))

cache = Cache(
    config={
        "CACHE_TYPE": "filesystem",
        "CACHE_DIR": cache_dir,
        "CACHE_DEFAULT_TIMEOUT": 0,  # cache indefinitely
        "CACHE_THRESHOLD": 1000,
        "CACHE_KEY_PREFIX": CACHE_TIMESTAMP,
    }
)


def cleanup_old_cache():
    """Remove all cache files"""
    try:
        cleanup_count = 0
        for filename in os.listdir(cache_dir):
            filepath = os.path.join(cache_dir, filename)
            os.remove(filepath)
            cleanup_count += 1
        logger.info(f"Cleaned up {cleanup_count} cache files")
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")

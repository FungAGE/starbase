from flask_caching import Cache
import os
import time
from src.config.settings import PROJECT_ROOT
from src.config.logging import get_logger

logger = get_logger(__name__)

cache_dir = os.path.join(PROJECT_ROOT, "src", "database", "cache")
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
    """Remove all cache files and directories"""
    try:
        cleanup_count = 0
        for filename in os.listdir(cache_dir):
            filepath = os.path.join(cache_dir, filename)
            if os.path.isfile(filepath):
                os.remove(filepath)
            elif os.path.isdir(filepath):
                import shutil

                shutil.rmtree(filepath)  # remove old directories if they exist
            cleanup_count += 1
        logger.debug(f"Cleaned up {cleanup_count} cache items")
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")

from flask_caching import Cache
import os
import time
from src.config.settings import PROJECT_ROOT, IS_DEV
from src.config.logging import get_logger
from functools import wraps
import pandas as pd

logger = get_logger(__name__)

cache_dir = os.path.join(PROJECT_ROOT, "src", "database", "cache")
os.makedirs(cache_dir, exist_ok=True)

CACHE_TIMESTAMP = str(int(time.time()))

DEFAULT_CACHE_TIMEOUT = 3600
LONG_CACHE_TIMEOUT = 86400
MAX_CACHE_SIZE = 512 * 1024 * 1024

cache = Cache(
    config={
        "CACHE_TYPE": "filesystem",
        "CACHE_DIR": cache_dir,
        "CACHE_DEFAULT_TIMEOUT": DEFAULT_CACHE_TIMEOUT,
        "CACHE_THRESHOLD": 1000,
        "CACHE_KEY_PREFIX": "starbase",
    }
)


def cleanup_old_cache(max_age_days=None):
    """Optional cache cleanup for persistent database data.

    Since our cached data is persistent with database versions and remains
    stable once generated at startup, cleanup is primarily for disk space management.

    Args:
        max_age_days: If provided, remove files older than this many days.
                     If None, only enforce size limits.
    """
    try:
        from pathlib import Path
        import time

        cleanup_count = 0
        total_size = 0
        current_time = time.time()
        cache_files = []

        # Collect all cache files with their metadata
        for filepath in Path(cache_dir).rglob("*"):
            if filepath.is_file() and filepath.suffix not in [
                ".lock"
            ]:  # Skip lock files
                try:
                    stats = filepath.stat()
                    age_days = (current_time - stats.st_mtime) / (24 * 3600)
                    cache_files.append(
                        (filepath, stats.st_size, stats.st_mtime, age_days)
                    )
                    total_size += stats.st_size
                except OSError:
                    continue

        # Remove files older than max_age_days if specified
        if max_age_days is not None:
            for filepath, size, _, age_days in cache_files[
                :
            ]:  # Copy list to avoid modification during iteration
                if age_days > max_age_days:
                    try:
                        filepath.unlink()
                        total_size -= size
                        cleanup_count += 1
                        logger.debug(
                            f"Removed old cache file (age: {age_days:.1f} days): {filepath}"
                        )
                    except OSError:
                        continue

        # If we're over size limit, remove oldest files until under limit
        if total_size > MAX_CACHE_SIZE:
            # Sort by modification time (oldest first)
            cache_files.sort(key=lambda x: x[2])

            for filepath, size, _, _ in cache_files:
                if total_size <= MAX_CACHE_SIZE:
                    break
                try:
                    filepath.unlink()
                    total_size -= size
                    cleanup_count += 1
                    logger.debug(f"Removed old cache file (size limit): {filepath}")
                except OSError:
                    continue

        if cleanup_count > 0:
            logger.info(
                f"Cleaned up {cleanup_count} cache items. Current cache size: {total_size / 1024 / 1024:.2f}MB"
            )

    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}", exc_info=True)


def cache_key_builder(*args, **kwargs):
    """Build a cache key from function arguments"""
    key_parts = [str(arg) for arg in args]
    key_parts.extend(f"{k}:{v}" for k, v in sorted(kwargs.items()))
    return "|".join(key_parts)


def smart_cache(timeout=3600, unless=None):
    """Smart caching decorator that handles pandas DataFrames"""

    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            if unless and unless(*args, **kwargs):
                return f(*args, **kwargs)
            actual_timeout = timeout
            if IS_DEV and timeout is None:
                actual_timeout = 300

            cache_key = f"{f.__name__}:{cache_key_builder(*args, **kwargs)}"
            result = cache.get(cache_key)

            if result is not None:
                if isinstance(result, dict) and "pandas_df" in result:
                    # Reconstruct DataFrame from cached dict
                    return pd.DataFrame.from_dict(result["pandas_df"])
                return result

            result = f(*args, **kwargs)

            if isinstance(result, pd.DataFrame):
                # Cache DataFrame as dictionary
                cache.set(
                    cache_key, {"pandas_df": result.to_dict()}, timeout=actual_timeout
                )
            else:
                cache.set(cache_key, result, timeout=actual_timeout)

            return result

        return wrapper

    return decorator

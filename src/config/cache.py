from flask_caching import Cache
import os
import time
from src.config.settings import PROJECT_ROOT
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


def cleanup_old_cache():
    """Smart cache cleanup that only removes expired items and maintains size limits"""
    try:
        import time
        import json
        from pathlib import Path
        
        cleanup_count = 0
        total_size = 0
        current_time = time.time()
        cache_files = []

        # Collect cache file information
        for filepath in Path(cache_dir).rglob('*'):
            if filepath.is_file():
                try:
                    stats = filepath.stat()
                    # Check if file is a cache file by attempting to read it
                    with open(filepath, 'rb') as f:
                        cache_data = cache._read_cache_file(f)
                        if cache_data:
                            expiry = cache_data.get('exp', 0)
                            if expiry and current_time >= expiry:
                                filepath.unlink()
                                cleanup_count += 1
                            else:
                                cache_files.append((filepath, stats.st_size, stats.st_mtime))
                                total_size += stats.st_size
                except (json.JSONDecodeError, IOError):
                    # Not a valid cache file or can't read it
                    continue

        # If we're over size limit, remove oldest files until under limit
        if total_size > MAX_CACHE_SIZE:
            # Sort by modification time (oldest first)
            cache_files.sort(key=lambda x: x[2])
            
            for filepath, size, _ in cache_files:
                if total_size <= MAX_CACHE_SIZE:
                    break
                try:
                    filepath.unlink()
                    total_size -= size
                    cleanup_count += 1
                except IOError:
                    continue

        logger.info(f"Cleaned up {cleanup_count} cache items. Current cache size: {total_size / 1024 / 1024:.2f}MB")
        
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
            
            cache_key = f"{f.__name__}:{cache_key_builder(*args, **kwargs)}"
            result = cache.get(cache_key)
            
            if result is not None:
                logger.debug(f"Cache hit for {cache_key}")
                if isinstance(result, dict) and "pandas_df" in result:
                    # Reconstruct DataFrame from cached dict
                    return pd.DataFrame.from_dict(result["pandas_df"])
                return result
            
            logger.debug(f"Cache miss for {cache_key}")
            result = f(*args, **kwargs)
            
            if isinstance(result, pd.DataFrame):
                # Cache DataFrame as dictionary
                cache.set(cache_key, {"pandas_df": result.to_dict()}, timeout=timeout)
            else:
                cache.set(cache_key, result, timeout=timeout)
            
            return result
        return wrapper
    return decorator
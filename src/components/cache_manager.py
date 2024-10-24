import os
import pickle
from typing import Any
import logging

logger = logging.getLogger(__name__)


CACHE_DIR = "cache"  # Directory where all cached objects will be saved

# Ensure the cache directory exists
if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR)
    logger.info(f"Cache directory created at {CACHE_DIR}")


def get_cache_filepath(cache_key: str) -> str:
    """Generate the full path for a cache file based on a cache key."""
    filepath = os.path.join(CACHE_DIR, f"{cache_key}.pkl")
    logger.debug(f"Cache file path for key '{cache_key}': {filepath}")
    return filepath


def save_to_cache(obj: Any, cache_key: str):
    """Save an object to a cache file."""
    filepath = get_cache_filepath(cache_key)
    try:
        with open(filepath, "wb") as f:
            pickle.dump(obj, f)
        logger.info(f"Cached object under key '{cache_key}' at {filepath}")
    except Exception as e:
        logger.error(f"Failed to save cache for key '{cache_key}': {str(e)}")


def load_from_cache(cache_key: str) -> Any:
    """Load an object from the cache."""
    filepath = get_cache_filepath(cache_key)
    if os.path.exists(filepath):
        try:
            with open(filepath, "rb") as f:
                logger.info(f"Loaded cached object for key '{cache_key}'")
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Failed to load cache for key '{cache_key}': {str(e)}")
    else:
        logger.warning(f"No cache file found for key '{cache_key}'")
    return None


def cache_exists(cache_key: str) -> bool:
    """Check if a cached file exists for a given cache key."""
    exists = os.path.exists(get_cache_filepath(cache_key))
    if exists:
        logger.debug(f"Cache exists for key '{cache_key}'")
    else:
        logger.debug(f"Cache does not exist for key '{cache_key}'")
    return exists


def invalidate_cache(cache_key: str):
    """Delete the cached file for a specific cache key."""
    filepath = get_cache_filepath(cache_key)
    if os.path.exists(filepath):
        os.remove(filepath)
        logger.info(
            f"Cache invalidated for key '{cache_key}' (file deleted: {filepath})"
        )
    else:
        logger.warning(f"No cache file found to invalidate for key '{cache_key}'")

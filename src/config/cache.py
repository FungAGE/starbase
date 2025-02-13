from flask_caching import Cache
import os
from src.config.settings import DATA_DIR
import logging

logger = logging.getLogger(__name__)

cache_dir = os.path.join(DATA_DIR, "cache")
os.makedirs(cache_dir, exist_ok=True)

cache = Cache(config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': cache_dir,
    'CACHE_DEFAULT_TIMEOUT': 86400,  # 24 hour cache
    'CACHE_THRESHOLD': 1000
})

def initialize_cache():
    """Initialize and warm up the cache with frequently accessed data"""
    # Import here to avoid circular dependency
    from src.config.precompute import precompute_tasks
    
    results = {}
    for key, fetch_func in precompute_tasks.items():
        try:
            data = cache.get(key)
            if data is None:
                data = fetch_func()
                if data is not None:
                    cache.set(key, data, timeout=86400)  # 24 hour cache
                    results[key] = True
                else:
                    results[key] = False
        except Exception as e:
            logger.error(f"Failed to cache {key}: {str(e)}")
            results[key] = False
    
    return results
from flask_caching import Cache
import os
from src.database.sql_manager import fetch_meta_data, fetch_paper_data, fetch_ship_table, fetch_all_ships, get_database_stats
import logging
from src.config.settings import DATA_DIR

logger = logging.getLogger(__name__)

cache_dir = os.path.join(DATA_DIR, "cache")
os.makedirs(cache_dir, exist_ok=True)

cache = Cache(config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': cache_dir,
    'CACHE_DEFAULT_TIMEOUT': 86400,  # 24 hours
    'CACHE_THRESHOLD': 1000,
    'CACHE_OPTIONS': {
        'mode': 0o777
    }
})

def initialize_cache():
    """Initialize and warm up the cache with frequently accessed data"""
    cache_keys = {
        "meta_data": fetch_meta_data,
        "paper_data": fetch_paper_data,
        "ship_table": fetch_ship_table,
        "all_ships": fetch_all_ships,
        "database_stats": get_database_stats
    }
    
    results = {}
    for key, fetch_func in cache_keys.items():
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
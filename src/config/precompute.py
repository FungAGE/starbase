import logging
from src.config.cache import cache
from src.database.sql_manager import (
    fetch_meta_data,
    fetch_paper_data,
    fetch_ship_table,
    fetch_all_ships,
    get_database_stats,
    fetch_download_data
)

logger = logging.getLogger(__name__)

precompute_tasks = {
    "meta_data": lambda: fetch_meta_data(curated=True),
    "paper_data": fetch_paper_data,
    "ship_table": fetch_ship_table,
    "all_ships": fetch_all_ships,
    "database_stats": get_database_stats,
    "download_data_curated_true_derep_false": lambda: fetch_download_data(curated=True, dereplicate=False),
    "download_data_curated_true_derep_true": lambda: fetch_download_data(curated=True, dereplicate=True),
    "download_data_curated_false_derep_false": lambda: fetch_download_data(curated=False, dereplicate=False),
    "download_data_curated_false_derep_true": lambda: fetch_download_data(curated=False, dereplicate=True),
}

def precompute_all():
    """Precompute and cache all necessary data and figures."""
    cache_status = {}
    
    try:
        for key, func in precompute_tasks.items():
            try:
                logger.info(f"Precomputing {key}...")
                result = func()
                if result is not None:
                    cache.set(key, result)
                    cache_status[key] = True
                    logger.info(f"Successfully cached {key}")
                else:
                    cache_status[key] = False
                    logger.error(f"Failed to precompute {key}: returned None")
            except Exception as e:
                logger.error(f"Failed to precompute {key}: {str(e)}")
                cache_status[key] = False
                    
    except Exception as e:
        logger.error(f"Critical error during precomputation: {str(e)}")
    
    # Log final status
    success_rate = sum(1 for v in cache_status.values() if v) / len(cache_status)
    logger.info(f"Precomputation complete. Success rate: {success_rate:.1%}")
    logger.info(f"Cache status: {cache_status}")
    
    return cache_status
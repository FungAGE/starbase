import logging
from src.config.cache import cache

logger = logging.getLogger(__name__)

from src.config.constants import precompute_tasks

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
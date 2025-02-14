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
from flask_caching import Cache
import os

cache_dir = os.path.join(os.environ.get('HOME', '/tmp'), '.cache')
if not os.path.exists(cache_dir):
    os.makedirs(cache_dir)

cache = Cache(config={
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': cache_dir,
    'CACHE_DEFAULT_TIMEOUT': 300,
    'CACHE_THRESHOLD': 1000
})
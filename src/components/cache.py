from flask_caching import Cache

# Initialize the cache (the app will bind to it later)
cache = Cache(config={"CACHE_TYPE": "SimpleCache", "CACHE_DEFAULT_TIMEOUT": 300})

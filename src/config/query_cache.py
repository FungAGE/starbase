"""In-process query cache for the compute backend (no Flask app context required)."""

from __future__ import annotations

import os
from functools import wraps
from typing import Any, Callable

import pandas as pd
from cachetools import LRUCache, TTLCache

from src.config.logging import get_logger

logger = get_logger(__name__)

_MAXSIZE = int(os.getenv("BACKEND_QUERY_CACHE_MAXSIZE", "32"))
_DEFAULT_TTL = int(os.getenv("BACKEND_QUERY_CACHE_TTL", "3600"))

_persistent_cache: LRUCache[str, Any] = LRUCache(maxsize=_MAXSIZE)
_ttl_caches: dict[int, TTLCache] = {}


def _cache_key(name: str, args: tuple, kwargs: dict) -> str:
    parts = [name] + [str(a) for a in args]
    parts.extend(f"{k}:{v}" for k, v in sorted(kwargs.items()))
    return "|".join(parts)


def _get_ttl_cache(ttl: int) -> TTLCache:
    if ttl not in _ttl_caches:
        _ttl_caches[ttl] = TTLCache(maxsize=_MAXSIZE, ttl=ttl)
    return _ttl_caches[ttl]


def query_cache(timeout: int | None = 3600, key_prefix: str | None = None):
    """Cache function results in memory; serializes pandas DataFrames to dicts."""

    def decorator(f: Callable):
        @wraps(f)
        def wrapper(*args, **kwargs):
            prefix = key_prefix or f.__name__
            cache_key = _cache_key(prefix, args, kwargs)
            store = _persistent_cache if timeout is None else _get_ttl_cache(timeout)

            cached = store.get(cache_key)
            if cached is not None:
                if isinstance(cached, dict) and "pandas_df" in cached:
                    return pd.DataFrame.from_dict(cached["pandas_df"])
                return cached

            result = f(*args, **kwargs)
            if isinstance(result, pd.DataFrame):
                store[cache_key] = {"pandas_df": result.to_dict()}
            else:
                store[cache_key] = result
            return result

        return wrapper

    return decorator

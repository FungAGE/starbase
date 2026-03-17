#!/usr/bin/env python3
"""
Cache clearing utility for starbase development.
Run this when you need to clear cached database queries.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from pathlib import Path
from src.config.settings import CACHE_DIR


def clear_cache():
    """Clear all cache files"""
    cache_path = Path(CACHE_DIR)
    if cache_path.exists():
        cleared_count = 0
        for cache_file in cache_path.glob("*"):
            if cache_file.is_file():
                try:
                    cache_file.unlink()
                    cleared_count += 1
                    print(f"Cleared: {cache_file.name}")
                except Exception as e:
                    print(f"Failed to clear {cache_file.name}: {e}")

        print(f"\n✅ Cleared {cleared_count} cache files")
    else:
        print(f"Cache directory does not exist: {CACHE_DIR}")


if __name__ == "__main__":
    clear_cache()

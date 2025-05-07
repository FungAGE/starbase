from flask_caching import Cache
import os
import time


from src.config.logging import get_logger

logger = get_logger(__name__)

# Use /dev/shm for temporary storage
cache_dir = "/dev/shm/starbase_cache"
os.makedirs(cache_dir, exist_ok=True)

# Create tmp directory
tmp_dir = os.path.join(cache_dir, "tmp")
os.makedirs(tmp_dir, exist_ok=True)

CACHE_TIMESTAMP = str(int(time.time()))

cache = Cache(
    config={
        "CACHE_TYPE": "filesystem",
        "CACHE_DIR": cache_dir,
        "CACHE_DEFAULT_TIMEOUT": 0,  # cache indefinitely
        "CACHE_THRESHOLD": 1000,
        "CACHE_KEY_PREFIX": CACHE_TIMESTAMP,
    }
)


def cleanup_old_cache(max_age_hours=24):
    """Remove cache files older than max_age_hours"""
    try:
        now = time.time()
        cleanup_count = 0
        skipped_count = 0

        # Walk through all subdirectories in cache_dir
        for root, dirs, files in os.walk(cache_dir):
            for filename in files:
                filepath = os.path.join(root, filename)

                # Skip if it's a directory or not a regular file
                if not os.path.isfile(filepath):
                    continue

                # Check file age
                try:
                    file_age_hours = (now - os.path.getmtime(filepath)) / 3600

                    if file_age_hours > max_age_hours:
                        # This file is older than our threshold
                        os.remove(filepath)
                        cleanup_count += 1
                    else:
                        skipped_count += 1
                except Exception as e:
                    # Handle errors for individual files without stopping the process
                    logger.warning(f"Error processing file {filepath}: {str(e)}")

        logger.info(
            f"Cache cleanup: removed {cleanup_count} files older than {max_age_hours} hours, kept {skipped_count} recent files"
        )
        return {"removed": cleanup_count, "kept": skipped_count}
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return {"error": str(e)}


def cleanup_by_directory_policy():
    """Clean cache files with directory-specific age policies"""
    policies = {
        os.path.join(cache_dir, "tmp"): 24,  # Temp files: 1 day
    }

    total_removed = 0
    for directory, hours in policies.items():
        if os.path.exists(directory):
            result = cleanup_directory(directory, max_age_hours=hours)
            total_removed += result.get("removed", 0)

    return {"removed": total_removed}


def cleanup_directory(directory, max_age_hours=24):
    """
    Clean files in a specific directory that are older than max_age_hours.
    Does not recurse into subdirectories.

    Args:
        directory (str): Directory path to clean
        max_age_hours (int): Maximum age of files to keep in hours

    Returns:
        dict: Statistics about the cleanup operation
    """
    try:
        now = time.time()
        cleanup_count = 0
        skipped_count = 0

        # Check if directory exists
        if not os.path.exists(directory):
            logger.warning(f"Directory does not exist: {directory}")
            return {"removed": 0, "kept": 0, "error": "Directory does not exist"}

        # Process only files in this directory (no subdirectories)
        for filename in os.listdir(directory):
            filepath = os.path.join(directory, filename)

            # Skip if it's not a regular file
            if not os.path.isfile(filepath):
                continue

            # Check file age
            try:
                file_age_hours = (now - os.path.getmtime(filepath)) / 3600

                if file_age_hours > max_age_hours:
                    # This file is older than our threshold
                    os.remove(filepath)
                    cleanup_count += 1
                else:
                    skipped_count += 1
            except Exception as e:
                # Handle errors for individual files without stopping the process
                logger.warning(f"Error processing file {filepath}: {str(e)}")

        logger.info(
            f"Directory cleanup ({directory}): removed {cleanup_count} files older than {max_age_hours} hours, kept {skipped_count} recent files"
        )
        return {"removed": cleanup_count, "kept": skipped_count, "directory": directory}
    except Exception as e:
        logger.error(f"Directory cleanup failed for {directory}: {str(e)}")
        return {"removed": 0, "kept": 0, "error": str(e), "directory": directory}

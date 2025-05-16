"""
Simple scheduler to replace Celery Beat.
Uses threading and schedule to run periodic tasks.
"""

import threading
import time
import schedule
from src.config.settings import IPSTACK_API_KEY
from src.config.logging import get_logger
from src.telemetry.tasks import (
    update_ip_locations,
    refresh_telemetry,
    check_cache_status,
)
from src.tasks import cleanup_cache_task

logger = get_logger(__name__)


def run_scheduler():
    """Run the scheduler in a separate thread."""
    scheduler_thread = threading.Thread(target=scheduler_worker, daemon=True)
    scheduler_thread.start()
    logger.info("Scheduler started")
    return scheduler_thread


def scheduler_worker():
    """Worker function for the scheduler thread."""
    # Define scheduled tasks

    # Every hour (at minute 0)
    schedule.every().hour.at(":00").do(run_update_ip_locations)

    # Every 15 minutes
    schedule.every(15).minutes.do(run_refresh_telemetry)

    # Every 5 minutes
    schedule.every(5).minutes.do(run_check_cache_status)

    # Every hour
    schedule.every().hour.do(run_cleanup_cache)

    logger.info("Scheduled tasks initialized")

    # Run the scheduler loop
    while True:
        schedule.run_pending()
        time.sleep(30)  # Sleep for 30 seconds between checks


def run_update_ip_locations():
    """Run the update_ip_locations task."""
    try:
        logger.info("Running update_ip_locations task")
        result = update_ip_locations(IPSTACK_API_KEY)
        logger.info(f"update_ip_locations completed: {result}")
        return result
    except Exception as e:
        logger.error(f"Error running update_ip_locations: {str(e)}")
        return {"status": "error", "message": str(e)}


def run_refresh_telemetry():
    """Run the refresh_telemetry task."""
    try:
        logger.info("Running refresh_telemetry task")
        result = refresh_telemetry()
        logger.info(f"refresh_telemetry completed: {result}")
        return result
    except Exception as e:
        logger.error(f"Error running refresh_telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}


def run_check_cache_status():
    """Run the check_cache_status task."""
    try:
        logger.info("Running check_cache_status task")
        result = check_cache_status()
        logger.info(f"check_cache_status completed: {result}")
        return result
    except Exception as e:
        logger.error(f"Error running check_cache_status: {str(e)}")
        return {"status": "error", "message": str(e)}


def run_cleanup_cache():
    """Run the cleanup_cache task."""
    try:
        logger.info("Running cleanup_cache task")
        result = cleanup_cache_task()
        logger.info(f"cleanup_cache completed: {result}")
        return result
    except Exception as e:
        logger.error(f"Error running cleanup_cache: {str(e)}")
        return {"status": "error", "message": str(e)}

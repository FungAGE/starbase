#!/usr/bin/env python3
"""
Script to start Celery worker and beat scheduler for the Starbase application.
This script can be used to run Celery processes in a Kubernetes pod or locally.
"""

import os
import sys
import subprocess
import signal
import time

# Add the project root to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.config.logging import get_logger

logger = get_logger(__name__)

# Global variables to track processes
worker_process = None
beat_process = None


def start_celery_worker():
    """Start the Celery worker process."""
    global worker_process

    try:
        logger.info("Starting Celery worker...")
        worker_process = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "celery",
                "-A",
                "src.config.celery_config",
                "worker",
                "--loglevel=info",
                "--concurrency=2",
            ]
        )
        logger.info(f"Celery worker started with PID: {worker_process.pid}")
        return worker_process
    except Exception as e:
        logger.error(f"Failed to start Celery worker: {e}")
        return None


def start_celery_beat():
    """Start the Celery beat scheduler process."""
    global beat_process

    try:
        logger.info("Starting Celery beat scheduler...")
        beat_process = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "celery",
                "-A",
                "src.config.celery_config",
                "beat",
                "--loglevel=info",
            ]
        )
        logger.info(f"Celery beat scheduler started with PID: {beat_process.pid}")
        return beat_process
    except Exception as e:
        logger.error(f"Failed to start Celery beat scheduler: {e}")
        return None


def stop_processes():
    """Stop all Celery processes gracefully."""
    global worker_process, beat_process

    logger.info("Stopping Celery processes...")

    if beat_process:
        logger.info(f"Stopping beat scheduler (PID: {beat_process.pid})")
        beat_process.terminate()
        try:
            beat_process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            logger.warning("Beat scheduler didn't stop gracefully, forcing...")
            beat_process.kill()

    if worker_process:
        logger.info(f"Stopping worker (PID: {worker_process.pid})")
        worker_process.terminate()
        try:
            worker_process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            logger.warning("Worker didn't stop gracefully, forcing...")
            worker_process.kill()

    logger.info("All Celery processes stopped")


def signal_handler(signum, frame):
    """Handle shutdown signals."""
    logger.info(f"Received signal {signum}, shutting down...")
    stop_processes()
    sys.exit(0)


def main():
    """Main function to start Celery processes."""
    # Set up signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    logger.info("Starting Starbase Celery services...")

    # Start the worker
    worker = start_celery_worker()
    if not worker:
        logger.error("Failed to start worker, exiting")
        sys.exit(1)

    # Start the beat scheduler
    beat = start_celery_beat()
    if not beat:
        logger.error("Failed to start beat scheduler, exiting")
        stop_processes()
        sys.exit(1)

    logger.info("All Celery processes started successfully")

    try:
        # Keep the main process alive and monitor child processes
        while True:
            # Check if processes are still running
            if worker.poll() is not None:
                logger.error("Worker process died unexpectedly")
                break

            if beat.poll() is not None:
                logger.error("Beat scheduler process died unexpectedly")
                break

            time.sleep(5)  # Check every 5 seconds

    except KeyboardInterrupt:
        logger.info("Received keyboard interrupt")
    finally:
        stop_processes()


if __name__ == "__main__":
    main()

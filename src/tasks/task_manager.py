"""
Task manager for handling async Celery tasks with UI integration.
Provides a bridge between synchronous UI callbacks and asynchronous Celery tasks.
"""

import uuid
import time
from typing import Dict, Any, Optional
from celery.result import AsyncResult
from src.config.celery_config import celery
from src.config.logging import get_logger

logger = get_logger(__name__)

# Global task status store (in production, this would be Redis or a database)
task_status_store: Dict[str, Dict[str, Any]] = {}


class TaskManager:
    """Manager for handling async task execution and status tracking."""

    @staticmethod
    def submit_task(task_func, *args, **kwargs) -> str:
        """
        Submit a task for async execution and return a task ID.

        Args:
            task_func: The Celery task function to execute
            *args: Arguments to pass to the task
            **kwargs: Keyword arguments to pass to the task

        Returns:
            str: Task ID that can be used to check status
        """
        task_id = str(uuid.uuid4())

        try:
            # Submit the task to Celery
            async_result = task_func.delay(*args, **kwargs)

            # Store task info
            task_status_store[task_id] = {
                "celery_task_id": async_result.id,
                "status": "PENDING",
                "result": None,
                "error": None,
                "started_at": time.time(),
                "completed_at": None,
                "task_name": task_func.name,
            }

            logger.info(
                f"Submitted task {task_func.name} with ID {task_id} (Celery ID: {async_result.id})"
            )
            return task_id

        except Exception as e:
            # If task submission fails, store the error
            task_status_store[task_id] = {
                "celery_task_id": None,
                "status": "FAILURE",
                "result": None,
                "error": str(e),
                "started_at": time.time(),
                "completed_at": time.time(),
                "task_name": task_func.name,
            }
            logger.error(f"Failed to submit task {task_func.name}: {e}")
            return task_id

    @staticmethod
    def get_task_status(task_id: str) -> Dict[str, Any]:
        """
        Get the current status of a task.

        Args:
            task_id: The task ID returned by submit_task

        Returns:
            Dict containing task status information
        """
        if task_id not in task_status_store:
            return {"status": "NOT_FOUND", "error": f"Task {task_id} not found"}

        task_info = task_status_store[task_id]
        celery_task_id = task_info.get("celery_task_id")

        if celery_task_id:
            try:
                # Get the actual status from Celery
                async_result = AsyncResult(celery_task_id, app=celery)

                # Update our status based on Celery's status
                if async_result.ready():
                    if async_result.successful():
                        task_info["status"] = "SUCCESS"
                        task_info["result"] = async_result.result
                        task_info["completed_at"] = time.time()
                    else:
                        task_info["status"] = "FAILURE"
                        task_info["error"] = str(async_result.info)
                        task_info["completed_at"] = time.time()
                else:
                    # Task is still running
                    if async_result.state == "RETRY":
                        task_info["status"] = "RETRY"
                    else:
                        task_info["status"] = "PENDING"

            except Exception as e:
                logger.error(f"Error checking task status for {task_id}: {e}")
                task_info["status"] = "FAILURE"
                task_info["error"] = f"Error checking task status: {e}"
                task_info["completed_at"] = time.time()

        return task_info.copy()

    @staticmethod
    def wait_for_task(
        task_id: str, timeout: int = 300, poll_interval: float = 1.0
    ) -> Dict[str, Any]:
        """
        Wait for a task to complete (for synchronous usage).

        Args:
            task_id: The task ID to wait for
            timeout: Maximum time to wait in seconds
            poll_interval: How often to check status in seconds

        Returns:
            Dict containing final task status
        """
        start_time = time.time()

        while time.time() - start_time < timeout:
            status = TaskManager.get_task_status(task_id)

            if status["status"] in ["SUCCESS", "FAILURE", "NOT_FOUND"]:
                return status

            time.sleep(poll_interval)

        # Timeout reached
        return {
            "status": "TIMEOUT",
            "error": f"Task {task_id} timed out after {timeout} seconds",
        }

    @staticmethod
    def cleanup_old_tasks(max_age: int = 3600):
        """
        Clean up old task records.

        Args:
            max_age: Maximum age of completed tasks in seconds
        """
        current_time = time.time()
        to_remove = []

        for task_id, task_info in task_status_store.items():
            completed_at = task_info.get("completed_at")
            if completed_at and (current_time - completed_at) > max_age:
                to_remove.append(task_id)

        for task_id in to_remove:
            del task_status_store[task_id]
            logger.debug(f"Cleaned up old task record: {task_id}")


# Convenience functions for specific tasks
def submit_blast_search(
    query_header: str, query_seq: str, query_type: str, eval_threshold: float = 0.01
) -> str:
    """Submit a BLAST search task."""
    from src.tasks import run_blast_search_task

    return TaskManager.submit_task(
        run_blast_search_task, query_header, query_seq, query_type, eval_threshold
    )


def submit_hmmer_search(
    query_header: str, query_seq: str, query_type: str, eval_threshold: float = 0.01
) -> str:
    """Submit an HMMER search task."""
    from src.tasks import run_hmmer_search_task

    return TaskManager.submit_task(
        run_hmmer_search_task, query_header, query_seq, query_type, eval_threshold
    )


def submit_classification_workflow(
    class_dict: Dict[str, Any], meta_dict: Optional[Dict[str, Any]] = None
) -> str:
    """Submit a classification workflow task."""
    from src.tasks import run_classification_workflow_task

    return TaskManager.submit_task(
        run_classification_workflow_task, class_dict, meta_dict
    )


def submit_multi_pgv(
    gff_files: list, seqs: list, tmp_file: str, len_thr: int, id_thr: int
) -> str:
    """Submit a multi PGV task."""
    from src.tasks import run_multi_pgv_task

    return TaskManager.submit_task(
        run_multi_pgv_task, gff_files, seqs, tmp_file, len_thr, id_thr
    )


def submit_single_pgv(gff_file: str, tmp_file: str) -> str:
    """Submit a single PGV task."""
    from src.tasks import run_single_pgv_task

    return TaskManager.submit_task(run_single_pgv_task, gff_file, tmp_file)


def submit_telemetry_refresh(ipstack_api_key: str) -> str:
    """Submit a telemetry refresh task."""
    from src.tasks import refresh_telemetry_task

    return TaskManager.submit_task(refresh_telemetry_task, ipstack_api_key)


def submit_cache_cleanup() -> str:
    """Submit a cache cleanup task."""
    from src.tasks import cleanup_cache_task

    return TaskManager.submit_task(cleanup_cache_task)

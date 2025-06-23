"""
Synchronous adapter for Celery tasks.
Provides backward compatibility by allowing tasks to run synchronously when Celery is not available.
"""

import os
from typing import Any, Dict, Optional
from src.config.logging import get_logger
from src.tasks.task_manager import TaskManager

logger = get_logger(__name__)

# Check if Celery should be used (based on environment variables)
USE_CELERY = bool(os.getenv("CELERY_BROKER_URL") and os.getenv("CELERY_RESULT_BACKEND"))


def run_task_sync(task_func, *args, **kwargs) -> Dict[str, Any]:
    """
    Run a task either via Celery (async) or directly (sync) based on configuration.

    Args:
        task_func: The task function to execute
        *args: Arguments to pass to the task
        **kwargs: Keyword arguments to pass to the task

    Returns:
        Dict containing the task result
    """
    if USE_CELERY:
        try:
            # Try to use Celery
            task_id = TaskManager.submit_task(task_func, *args, **kwargs)

            # Wait for completion with a reasonable timeout
            timeout = kwargs.pop("timeout", 300)  # 5 minutes default
            result = TaskManager.wait_for_task(task_id, timeout=timeout)

            if result["status"] == "SUCCESS":
                return result["result"]
            elif result["status"] == "TIMEOUT":
                logger.warning(
                    f"Task {task_func.name} timed out after {timeout} seconds"
                )
                return None
            else:
                logger.error(
                    f"Task {task_func.name} failed: {result.get('error', 'Unknown error')}"
                )
                return None

        except Exception as e:
            logger.error(f"Failed to execute task {task_func.name} via Celery: {e}")
            logger.info("Falling back to synchronous execution")
            # Fall through to synchronous execution

    # Execute synchronously (either by choice or as fallback)
    try:
        logger.info(f"Executing task {task_func.name} synchronously")
        # For synchronous execution, we need to call the actual function without 'self'
        # Remove the 'self' parameter that Celery tasks expect
        if hasattr(task_func, "delay"):
            # This is a Celery task, we need to call the underlying function
            # Get the original function from the task
            func = (
                task_func.__wrapped__
                if hasattr(task_func, "__wrapped__")
                else task_func
            )
            if hasattr(func, "__func__"):
                func = func.__func__
                # Remove 'self' from args since it's not needed for direct execution
                if args and hasattr(args[0], "__self__"):
                    args = args[1:]
            return func(*args, **kwargs)
        else:
            return task_func(*args, **kwargs)

    except Exception as e:
        logger.error(f"Task {task_func.name} failed in synchronous execution: {e}")
        return None


# Convenience wrappers for the main tasks
def run_blast_search_sync(
    query_header: str, query_seq: str, query_type: str, eval_threshold: float = 0.01
):
    """Run BLAST search task synchronously or asynchronously based on configuration."""
    from src.tasks import run_blast_search_task

    return run_task_sync(
        run_blast_search_task, query_header, query_seq, query_type, eval_threshold
    )


def run_hmmer_search_sync(
    query_header: str, query_seq: str, query_type: str, eval_threshold: float = 0.01
):
    """Run HMMER search task synchronously or asynchronously based on configuration."""
    from src.tasks import run_hmmer_search_task

    return run_task_sync(
        run_hmmer_search_task, query_header, query_seq, query_type, eval_threshold
    )


def run_classification_workflow_sync(
    class_dict: Dict[str, Any], meta_dict: Optional[Dict[str, Any]] = None
):
    """Run classification workflow task synchronously or asynchronously based on configuration."""
    from src.tasks import run_classification_workflow_task

    # For synchronous execution, we need to ensure the class_dict is properly formatted
    # If it's a BlastData object, convert to dict for consistent handling
    if hasattr(class_dict, "__dict__"):
        # Convert BlastData object to dict for serialization compatibility
        class_dict_serialized = {
            "seq_type": getattr(class_dict, "seq_type", "nucl"),
            "fasta_file": getattr(class_dict, "fasta_file", None),
            "blast_df": getattr(class_dict, "blast_df", None),
            "fetch_ship_params": getattr(class_dict.fetch_ship_params, "__dict__", {})
            if hasattr(class_dict, "fetch_ship_params")
            else {},
            "fetch_captain_params": getattr(
                class_dict.fetch_captain_params, "__dict__", {}
            )
            if hasattr(class_dict, "fetch_captain_params")
            else {},
        }
        return run_task_sync(
            run_classification_workflow_task, class_dict_serialized, meta_dict
        )
    else:
        return run_task_sync(run_classification_workflow_task, class_dict, meta_dict)


def run_multi_pgv_sync(
    gff_files: list, seqs: list, tmp_file: str, len_thr: int, id_thr: int
):
    """Run multi PGV task synchronously or asynchronously based on configuration."""
    from src.tasks import run_multi_pgv_task

    return run_task_sync(run_multi_pgv_task, gff_files, seqs, tmp_file, len_thr, id_thr)


def run_single_pgv_sync(gff_file: str, tmp_file: str):
    """Run single PGV task synchronously or asynchronously based on configuration."""
    from src.tasks import run_single_pgv_task

    return run_task_sync(run_single_pgv_task, gff_file, tmp_file)

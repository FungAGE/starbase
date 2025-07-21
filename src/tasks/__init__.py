import tempfile
import os
import json

from src.config.cache import cache, cleanup_old_cache
from src.telemetry.utils import update_ip_locations
from src.utils.seq_utils import write_temp_fasta
from src.utils.blast_utils import run_blast, run_hmmer
from src.config.logging import get_logger

logger = get_logger(__name__)

# Make sure this is at the top level of the module
__all__ = [
    "refresh_telemetry_task",
    "cleanup_cache_task",
    "run_blast_search_task",
    "run_hmmer_search_task",
    "run_multi_pgv_task",
    "run_single_pgv_task",
    "run_classification_workflow_task",
]


def refresh_telemetry_task(ipstack_api_key):
    """Task to refresh telemetry data (formerly Celery task)"""
    try:
        update_ip_locations(ipstack_api_key)
        cache.delete("telemetry_data")
        return {"status": "success"}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}


def cleanup_cache_task():
    """Task to clean up cache files (formerly Celery task)"""
    try:
        cleanup_old_cache()
        return {"status": "success", "message": "Cache cleanup completed"}
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return {"status": "error", "message": str(e)}


# @celery.task(name="run_blast_search", bind=True, max_retries=3, retry_backoff=True)
def run_blast_search_task(query_header, query_seq, query_type, eval_threshold=0.01):
    """Celery task to run BLAST search"""

    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name

        # Run BLAST
        blast_results_file = run_blast(
            query_type=query_type,
            query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=eval_threshold,
            threads=2,
        )

        if blast_results_file is None:
            return None

        # Read the file contents instead of returning the path
        with open(blast_results_file, "r") as f:
            blast_results_content = f.read()

        # Create a temp file name to help consumers know what type of file this was
        temp_name = os.path.basename(blast_results_file)

        # Clean up the temporary file after reading
        try:
            os.unlink(blast_results_file)
        except Exception as e:
            logger.error(f"Error cleaning up BLAST results file: {str(e)}")

        # Return the file contents and metadata instead of the path
        return {
            "content": blast_results_content,
            "original_filename": temp_name,
            "file_type": "blast",
        }

    except Exception as e:
        logger.error(f"BLAST search failed: {str(e)}")
        return None


# @celery.task(name="run_hmmer_search", bind=True, max_retries=3, retry_backoff=True)
def run_hmmer_search_task(query_header, query_seq, query_type, eval_threshold=0.01):
    """Celery task to run HMMER search"""

    try:
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)

        results_dict, protein_filename = run_hmmer(
            query_type=query_type,
            input_gene="tyr",
            input_eval=eval_threshold,
            query_fasta=tmp_query_fasta,
            threads=2,
        )

        # If protein_filename exists, read its content
        protein_content = None
        if protein_filename and os.path.exists(protein_filename):
            with open(protein_filename, "r") as f:
                protein_content = f.read()

            # Clean up the temporary file
            try:
                os.unlink(protein_filename)
            except Exception as e:
                logger.error(f"Error cleaning up HMMER results file: {str(e)}")

        return {
            "results": results_dict,
            "protein_content": protein_content,
            "original_filename": os.path.basename(protein_filename)
            if protein_filename
            else None,
        }

    except Exception as e:
        logger.error(f"HMMER search failed: {str(e)}")
        return None


def run_multi_pgv_task(gff_files, seqs, tmp_file, len_thr, id_thr):
    from src.pages.pgv import multi_pgv

    try:
        return multi_pgv(gff_files, seqs, tmp_file, len_thr, id_thr)
    except Exception as e:
        logger.error(f"Multi PGV failed: {str(e)}")
        return None


def run_single_pgv_task(gff_file, tmp_file):
    """Task to run `single_pgv` (formerly Celery task)"""
    from src.pages.pgv import single_pgv

    try:
        return single_pgv(gff_file, tmp_file)
    except Exception as e:
        logger.error(f"Single PGV failed: {str(e)}")
        return None


# @celery.task(name="run_classification_workflow_task", bind=True)
def run_classification_workflow_task(class_dict, meta_dict=None):
    """Function to run the classification workflow.

    The workflow runs tasks sequentially.
    """
    from src.utils.classification_utils import run_classification_workflow
    from src.utils.blast_data import WorkflowState

    try:
        # Run the sequential workflow
        result = run_classification_workflow(class_dict, meta_dict)

        # If result is None or empty, return a default error state
        if result is None:
            logger.warning("Workflow returned None - returning failed state")
            return {
                "complete": True,
                "error": "Classification workflow returned no results",
                "status": "failed",
                "found_match": False,
                "match_stage": None,
                "match_result": None,
                "workflow_started": True,
                "current_stage": None,
                "current_stage_idx": 0,
                "start_time": 0.0,
                "stages": {},
                "class_dict": {},
                "task_id": "",
            }

        # Convert WorkflowState objects to dictionaries
        if isinstance(result, WorkflowState):
            logger.debug("Converting WorkflowState object to dictionary")
            result = result.to_dict()

        # Ensure result is a dictionary
        if not isinstance(result, dict):
            logger.error(f"Workflow result is not a dictionary: {type(result)}")
            return {
                "complete": True,
                "error": f"Invalid workflow result type: {type(result)}",
                "status": "failed",
                "found_match": False,
                "match_stage": None,
                "match_result": None,
                "workflow_started": True,
                "current_stage": None,
                "current_stage_idx": 0,
                "start_time": 0.0,
                "stages": {},
                "class_dict": {},
                "task_id": "",
            }

        # Test JSON serialization before returning
        try:
            # Test if result is JSON serializable
            json.dumps(result)
            logger.debug("Workflow result is JSON serializable")
        except TypeError as json_error:
            logger.error(f"Result is not JSON serializable: {json_error}")
            # Return a safe version with error information
            return {
                "complete": True,
                "error": f"Classification result is not JSON serializable: {str(json_error)}",
                "status": "failed",
                "found_match": False,
                "match_stage": None,
                "match_result": None,
                "workflow_started": True,
                "current_stage": None,
                "current_stage_idx": 0,
                "start_time": 0.0,
                "stages": {},
                "class_dict": {},
                "task_id": "",
            }

        return result
    except Exception as e:
        logger.error(f"Classification workflow task failed: {str(e)}")
        logger.exception("Full traceback:")
        return {
            "complete": True,
            "error": f"Classification workflow failed: {str(e)}",
            "status": "failed",
            "found_match": False,
            "match_stage": None,
            "match_result": None,
            "workflow_started": True,
            "current_stage": None,
            "current_stage_idx": 0,
            "start_time": 0.0,
            "stages": {},
            "class_dict": {},
            "task_id": "",
        }

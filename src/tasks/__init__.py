import tempfile
import os
import json

from src.config.cache import cache, cleanup_old_cache
from src.telemetry.utils import update_ip_locations
from src.utils.seq_utils import write_temp_fasta
from src.utils.blast_utils import run_blast, run_hmmer
from src.config.logging import get_logger
from src.config.celery_config import celery, CELERY_AVAILABLE

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


# Implementation functions
def _refresh_telemetry_impl(ipstack_api_key):
    """Implementation of refresh telemetry task"""
    try:
        update_ip_locations(ipstack_api_key)
        cache.delete("telemetry_data")
        return {"status": "success"}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}


def _cleanup_cache_impl():
    """Implementation of cleanup cache task"""
    try:
        cleanup_old_cache()
        return {"status": "success", "message": "Cache cleanup completed"}
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return {"status": "error", "message": str(e)}


# Create task wrappers or plain functions depending on Celery availability
if CELERY_AVAILABLE and celery:
    @celery.task(name="refresh_telemetry_task")
    def refresh_telemetry_task(ipstack_api_key):
        """Task to refresh telemetry data (Celery task)"""
        return _refresh_telemetry_impl(ipstack_api_key)

    @celery.task(name="src.tasks.cleanup_cache_task")
    def cleanup_cache_task():
        """Task to clean up cache files (Celery task)"""
        return _cleanup_cache_impl()
else:
    def refresh_telemetry_task(ipstack_api_key):
        """Task to refresh telemetry data (direct call)"""
        return _refresh_telemetry_impl(ipstack_api_key)

    def cleanup_cache_task():
        """Task to clean up cache files (direct call)"""
        return _cleanup_cache_impl()


def _run_blast_search_impl(query_header, query_seq, query_type, eval_threshold=0.01, curated=None):
    """Implementation of BLAST search task"""
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
        curated=curated,
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


if CELERY_AVAILABLE and celery:
    @celery.task(name="run_blast_search", bind=True, max_retries=3, retry_backoff=True)
    def run_blast_search_task(
        self, query_header, query_seq, query_type, eval_threshold=0.01, curated=None
    ):
        """Celery task to run BLAST search"""
        try:
            return _run_blast_search_impl(query_header, query_seq, query_type, eval_threshold, curated)
        except Exception as e:
            logger.error(f"BLAST search failed: {str(e)}")
            # Retry the task if it's a transient error
            if self.request.retries < self.max_retries:
                raise self.retry(countdown=60 * (2 ** self.request.retries))
            return None
else:
    def run_blast_search_task(query_header, query_seq, query_type, eval_threshold=0.01, curated=None):
        """Direct BLAST search (no Celery)"""
        try:
            return _run_blast_search_impl(query_header, query_seq, query_type, eval_threshold, curated)
        except Exception as e:
            logger.error(f"BLAST search failed: {str(e)}")
            return None


def _run_hmmer_search_impl(query_header, query_seq, query_type, eval_threshold=0.01):
    """Implementation of HMMER search task"""
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


if CELERY_AVAILABLE and celery:
    @celery.task(name="run_hmmer_search", bind=True, max_retries=3, retry_backoff=True)
    def run_hmmer_search_task(self, query_header, query_seq, query_type, eval_threshold=0.01):
        """Celery task to run HMMER search"""
        try:
            return _run_hmmer_search_impl(query_header, query_seq, query_type, eval_threshold)
        except Exception as e:
            logger.error(f"HMMER search failed: {str(e)}")
            # Retry the task if it's a transient error
            if self.request.retries < self.max_retries:
                raise self.retry(countdown=60 * (2 ** self.request.retries))
            return None
else:
    def run_hmmer_search_task(query_header, query_seq, query_type, eval_threshold=0.01):
        """Direct HMMER search (no Celery)"""
        try:
            return _run_hmmer_search_impl(query_header, query_seq, query_type, eval_threshold)
        except Exception as e:
            logger.error(f"HMMER search failed: {str(e)}")
            return None


if CELERY_AVAILABLE and celery:

    @celery.task(name="run_classification_workflow_task", bind=True)
    def run_classification_workflow_task(self, workflow_state, blast_data=None, classification_data=None, meta_dict=None):
        """
        Run the main classification workflow.
        This workflow runs the following tasks sequentially:
            1. Exact match check
            2. Contained match check
            3. Similarity match check
            4. Family classification
            5. Navis classification
            6. Haplotype classification

        Args:
            workflow_state_dict: WorkflowState object or dictionary
            blast_data_dict: BlastData object or dictionary
            classification_data_dict: ClassificationData object or dictionary
            meta_dict: MetaData object or dictionary
        """
        return _run_classification_workflow_internal(workflow_state, blast_data, classification_data, meta_dict)
else:
    def run_classification_workflow_task(workflow_state, blast_data=None, classification_data=None, meta_dict=None):
        """
        Run the main classification workflow directly (no Celery).
        This workflow runs the following tasks sequentially:
            1. Exact match check
            2. Contained match check
            3. Similarity match check
            4. Family classification
            5. Navis classification
            6. Haplotype classification

        Args:
            workflow_state_dict: WorkflowState object or dictionary
            blast_data_dict: BlastData object or dictionary
            classification_data_dict: ClassificationData object or dictionary
            meta_dict: MetaData object or dictionary
        """
        return _run_classification_workflow_internal(workflow_state, blast_data, classification_data, meta_dict)


def run_classification_workflow_sync(workflow_state, blast_data=None, classification_data=None, meta_dict=None):
    """
    Synchronous version of the classification workflow for direct calls.
    This is the same as the Celery task but without the task decorator.
    """
    return _run_classification_workflow_internal(workflow_state, blast_data, classification_data, meta_dict)


def _run_classification_workflow_internal(workflow_state, blast_data=None, classification_data=None, meta_dict=None):
    """
    Internal implementation of the classification workflow.
    """
    from src.utils.classification_utils import run_classification_workflow

    try:
        # Convert dictionaries to objects for the workflow function
        from src.utils.blast_data import WorkflowState, BlastData, ClassificationData
        
        workflow_state_obj = WorkflowState.from_dict(workflow_state) if isinstance(workflow_state, dict) else workflow_state
        blast_data_obj = BlastData.from_dict(blast_data) if isinstance(blast_data, dict) else blast_data
        classification_data_obj = ClassificationData.from_dict(classification_data) if isinstance(classification_data, dict) and classification_data else ClassificationData()
        
        # Run the sequential workflow
        result = run_classification_workflow(workflow_state_obj, blast_data_obj, classification_data_obj, meta_dict)

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

        # Ensure result is a dictionary
        if not isinstance(result, dict):
            # Convert WorkflowState object to dictionary if needed
            if hasattr(result, 'to_dict'):
                result = result.to_dict()
                logger.debug("Converted WorkflowState object to dictionary")
            else:
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

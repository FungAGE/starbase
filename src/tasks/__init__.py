from src.config.celery_config import celery
from src.config.cache import cache, cleanup_old_cache
from src.telemetry.utils import update_ip_locations
from src.utils.seq_utils import write_temp_fasta
from src.utils.blast_utils import run_blast, run_hmmer
import tempfile
import logging
import pandas as pd
from typing import Optional
import os
import json

logger = logging.getLogger(__name__)

# Make sure this is at the top level of the module
__all__ = [
    "refresh_telemetry_task",
    "cleanup_cache_task",
    "run_blast_search_task",
    "run_hmmer_search_task",
    "check_exact_matches_task",
    "check_contained_matches_task",
    "check_similar_matches_task",
    "run_family_classification_task",
    "run_navis_classification_task",
    "run_haplotype_classification_task",
    "run_metaeuk_easy_predict_task",
    "run_multi_pgv_task",
    "run_single_pgv_task",
    "run_classification_workflow_task",
]


@celery.task(name="refresh_telemetry")
def refresh_telemetry_task(ipstack_api_key):
    """Celery task to refresh telemetry data"""
    try:
        update_ip_locations(ipstack_api_key)
        cache.delete("telemetry_data")
        return {"status": "success"}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}


@celery.task(name="cleanup_cache")
def cleanup_cache_task():
    """Celery task to clean up cache files"""
    try:
        cleanup_old_cache()
        return {"status": "success", "message": "Cache cleanup completed"}
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return {"status": "error", "message": str(e)}


@celery.task(name="run_blast_search")
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


@celery.task(name="run_hmmer_search")
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


@celery.task(name="check_exact_matches_task")
def check_exact_matches_task(fasta: str, ships_dict: list) -> Optional[str]:
    from src.utils.classification_utils import check_exact_match

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        # Convert the list of dictionaries back to DataFrame
        existing_ships = pd.DataFrame.from_dict(ships_dict)
        return check_exact_match(fasta_path, existing_ships)
    except Exception as e:
        logger.error(f"Exact match check failed: {str(e)}")
        return None


@celery.task(name="check_contained_matches_task")
def check_contained_matches_task(fasta: str, ships_dict: list) -> Optional[str]:
    from src.utils.classification_utils import check_contained_match

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        # Convert the list of dictionaries back to DataFrame
        existing_ships = pd.DataFrame.from_dict(ships_dict)

        # Run the check
        result = check_contained_match(fasta_path, existing_ships=existing_ships)
        return result

    except Exception as e:
        logger.error(f"Contained match check failed: {str(e)}")
        return None


@celery.task(name="check_similar_matches_task")
def check_similar_matches_task(
    fasta: str, ships_dict: list, threshold: float = 0.9
) -> Optional[str]:
    from src.utils.classification_utils import check_similar_match

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        logger.debug(f"Starting similar match task with {len(ships_dict)} ships")
        # Convert list of dicts back to DataFrame
        existing_ships = pd.DataFrame.from_dict(ships_dict)
        logger.debug(f"Converted to DataFrame with shape {existing_ships.shape}")

        result = check_similar_match(fasta_path, existing_ships, threshold)
        logger.debug(f"Similar match result: {result}")
        return result
    except Exception as e:
        logger.error(f"Similar match check failed: {str(e)}")
        return None


@celery.task(name="run_family_classification_task")
def run_family_classification_task(fasta, seq_type):
    from src.utils.classification_utils import classify_family

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        family_dict, protein_file = classify_family(
            fasta=fasta_path,
            seq_type=seq_type,
            threads=1,
        )

        if family_dict:
            # Convert numpy types to Python native types
            family_dict = {
                "family": family_dict["family"],
                "aln_length": int(family_dict["aln_length"]),  # Convert np.int64 to int
                "evalue": float(family_dict["evalue"]),  # Convert np.float64 to float
            }

            # If protein file exists, read its content
            if protein_file and os.path.exists(protein_file):
                with open(protein_file, "r") as f:
                    protein_content = f.read()
                family_dict["protein_content"] = protein_content

                # Clean up temporary file
                try:
                    os.unlink(protein_file)
                except Exception as e:
                    logger.error(f"Error cleaning up HMMER results file: {str(e)}")

        return family_dict
    except Exception as e:
        logger.error(f"Family classification failed: {str(e)}")
        return None


@celery.task(name="run_navis_classification_task")
def run_navis_classification_task(fasta, existing_captains):
    from src.utils.classification_utils import classify_navis

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        # Convert existing_captains to DataFrame if it's a list of dicts
        if isinstance(existing_captains, list):
            existing_captains = pd.DataFrame.from_dict(existing_captains)

        return classify_navis(fasta_path, existing_captains)
    except Exception as e:
        logger.error(f"Navis classification failed: {str(e)}")
        return None


@celery.task(name="run_haplotype_classification_task")
def run_haplotype_classification_task(fasta, existing_ships, navis):
    from src.utils.classification_utils import classify_haplotype

    try:
        # Check if fasta is content rather than a file path
        if isinstance(fasta, dict) and "content" in fasta:
            # Create a temporary file with the content
            tmp_file = tempfile.NamedTemporaryFile(suffix=".fa", delete=False).name
            with open(tmp_file, "w") as f:
                f.write(fasta["content"])
            fasta_path = tmp_file
        else:
            fasta_path = fasta

        # Convert existing_ships to DataFrame if it's a list of dicts
        if isinstance(existing_ships, list):
            existing_ships = pd.DataFrame.from_dict(existing_ships)

        return classify_haplotype(fasta_path, existing_ships, navis)
    except Exception as e:
        logger.error(f"Haplotype classification failed: {str(e)}")
        return None


@celery.task(name="run_metaeuk_easy_predict_task")
def run_metaeuk_easy_predict_task(fasta, ref_db, output_prefix, threads):
    from src.utils.classification_utils import metaeuk_easy_predict

    try:
        return metaeuk_easy_predict(
            query_fasta=fasta,
            ref_db=ref_db,
            output_prefix=output_prefix,
            threads=threads,
        )
    except Exception as e:
        logger.error(f"Metaeuk failed: {str(e)}")
        return None


@celery.task(name="run_multi_pgv_task")
def run_multi_pgv_task(gff_files, seqs, tmp_file, len_thr, id_thr):
    from src.pages.pgv import multi_pgv

    try:
        return multi_pgv(gff_files, seqs, tmp_file, len_thr, id_thr)
    except Exception as e:
        logger.error(f"Multi PGV failed: {str(e)}")
        return None


@celery.task(name="run_single_pgv_task")
def run_single_pgv_task(gff_file, tmp_file):
    """Celery task to run `single_pgv`"""
    from src.pages.pgv import single_pgv

    try:
        return single_pgv(gff_file, tmp_file)
    except Exception as e:
        logger.error(f"Single PGV failed: {str(e)}")
        return None


@celery.task(name="run_classification_workflow_task", bind=True)
def run_classification_workflow_task(self, upload_data):
    """Celery task to run the classification workflow in the background.

    The workflow itself runs tasks sequentially, but the entire workflow
    runs as a background task to not block the UI.
    """
    from src.utils.classification_utils import run_classification_workflow

    try:
        # Run the sequential workflow
        result = run_classification_workflow(upload_data)

        # Test JSON serialization before returning
        try:
            # Test if result is JSON serializable
            json.dumps(result)
        except TypeError as json_error:
            logger.error(f"Result is not JSON serializable: {json_error}")
            # Return a safe version with error information
            return {
                "complete": True,
                "error": f"Classification result is not JSON serializable: {str(json_error)}",
                "status": "failed",
            }

        return result
    except Exception as e:
        logger.error(f"Classification workflow task failed: {str(e)}")
        logger.exception("Full traceback:")
        raise

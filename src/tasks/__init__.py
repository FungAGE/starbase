from src.config.celery_config import celery
from src.config.cache import cache, cleanup_old_cache
from src.utils.telemetry import maintain_ip_locations
from src.utils.seq_utils import write_temp_fasta
from src.utils.blast_utils import run_blast, run_hmmer
from src.config.settings import BLAST_DB_PATHS
import tempfile
import logging
import pandas as pd
from typing import Optional

logger = logging.getLogger(__name__)

# Make sure this is at the top level of the module
__all__ = [
    'refresh_telemetry_task',
    'cleanup_cache_task',
    'run_blast_search_task',
    'run_hmmer_search_task',
    'check_exact_matches_task',
    'check_contained_matches_task',
    'check_similar_matches_task',
    'check_exact_matches_task',
    'run_family_classification_task',
    'run_navis_classification_task',
    'run_haplotype_classification_task',
    'run_metaeuk_easy_predict_task',
    'run_multi_pgv_task',
    'run_single_pgv_task',
]

@celery.task(name='refresh_telemetry')
def refresh_telemetry_task(ipstack_api_key):
    """Celery task to refresh telemetry data"""
    try:
        maintain_ip_locations(ipstack_api_key)
        cache.delete('telemetry_data')
        return {"status": "success"}
    except Exception as e:
        logger.error(f"Error refreshing telemetry: {str(e)}")
        return {"status": "error", "message": str(e)}

@celery.task(name='cleanup_cache')
def cleanup_cache_task():
    """Celery task to clean up cache files"""
    try:
        cleanup_old_cache()
        return {"status": "success", "message": "Cache cleanup completed"}
    except Exception as e:
        logger.error(f"Cache cleanup failed: {str(e)}")
        return {"status": "error", "message": str(e)}

@celery.task(name='run_blast_search')
def run_blast_search_task(query_header, query_seq, query_type, eval_threshold=0.01):
    """Celery task to run BLAST search"""
    try:
        # Write sequence to temporary FASTA file
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        tmp_blast = tempfile.NamedTemporaryFile(suffix=".blast", delete=True).name

        # Run BLAST
        blast_results = run_blast(
            db_list=BLAST_DB_PATHS,
            query_type=query_type,
            query_fasta=tmp_query_fasta,
            tmp_blast=tmp_blast,
            input_eval=eval_threshold,
            threads=2,
        )
        
        if blast_results is None:
            return None
            
        return blast_results.to_dict("records")
        
    except Exception as e:
        logger.error(f"BLAST search failed: {str(e)}")
        return None

@celery.task(name='run_hmmer_search')
def run_hmmer_search_task(query_header, query_seq, query_type, eval_threshold=0.01):
    """Celery task to run HMMER search"""
    try:
        tmp_query_fasta = write_temp_fasta(query_header, query_seq)
        
        results_dict, protein_filename = run_hmmer(
            db_list=BLAST_DB_PATHS,
            query_type=query_type,
            input_gene="tyr",
            input_eval=eval_threshold,
            query_fasta=tmp_query_fasta,
            threads=2,
        )
        
        return {
            'results': results_dict,
            'protein_file': protein_filename
        }
        
    except Exception as e:
        logger.error(f"HMMER search failed: {str(e)}")
        return None 

@celery.task(name='check_exact_matches_task')
def check_exact_matches_task(fasta, ships_dict):
    """Celery task to run `check_exact_match`"""
    from src.utils.classification_utils import check_exact_match
    try:
        # Convert the list of dictionaries back to DataFrame
        existing_ships = pd.DataFrame(ships_dict)
        return check_exact_match(fasta, existing_ships)
    except Exception as e:
        logger.error(f"Exact match check failed: {str(e)}")
        return None

@celery.task(name='check_contained_matches_task')
def check_contained_matches_task(fasta: str, ships_dict: list) -> Optional[str]:
    """Celery task to run `check_contained_match`"""
    from src.utils.classification_utils import check_contained_match
    try:
        # Convert the list of dictionaries back to DataFrame
        existing_ships = pd.DataFrame.from_dict(ships_dict)
        
        # Run the check
        result = check_contained_match(fasta, existing_ships)
        return result
        
    except Exception as e:
        logger.error(f"Contained match check failed: {str(e)}")
        return None

@celery.task(name='check_similar_matches_task')
def check_similar_matches_task(fasta: str, ships_dict: list, threshold: float = 0.9) -> Optional[str]:
    """Celery task to run `check_similar_match`"""
    from src.utils.classification_utils import check_similar_match
    try:
        logger.info(f"Starting similar match task with {len(ships_dict)} ships")
        # Convert list of dicts back to DataFrame
        existing_ships = pd.DataFrame.from_dict(ships_dict)
        logger.info(f"Converted to DataFrame with shape {existing_ships.shape}")
        
        result = check_similar_match(fasta, existing_ships, threshold)
        logger.info(f"Similar match result: {result}")
        return result
    except Exception as e:
        logger.error(f"Similar match check failed: {str(e)}")
        return None

@celery.task(name='run_family_classification_task')
def run_family_classification_task(fasta, existing_ships):
    """Celery task to run `classify_family`"""
    from src.utils.classification_utils import classify_family
    try:
        return classify_family(fasta, existing_ships)
    except Exception as e:
        logger.error(f"Family classification failed: {str(e)}")
        return None

@celery.task(name='run_navis_classification_task')
def run_navis_classification_task(fasta, existing_ships):
    """Celery task to run `classify_navis`"""
    from src.utils.classification_utils import classify_navis
    try:
        return classify_navis(fasta, existing_ships)
    except Exception as e:
        logger.error(f"Navis classification failed: {str(e)}")
        return None

@celery.task(name='run_haplotype_classification_task')
def run_haplotype_classification_task(fasta, existing_ships):
    """Celery task to run `classify_haplotype`"""
    from src.utils.classification_utils import classify_haplotype
    try:
        return classify_haplotype(fasta, existing_ships)
    except Exception as e:
        logger.error(f"Haplotype classification failed: {str(e)}")
        return None

@celery.task(name="run_metaeuk_easy_predict_task")
def run_metaeuk_easy_predict_task(fasta, ref_db, output_prefix, threads):
    """Celery task to run `metaeuk`"""
    from src.utils.classification_utils import metaeuk_easy_predict
    try:
        return metaeuk_easy_predict(query_fasta=fasta, ref_db=ref_db, output_prefix=output_prefix, threads=threads)
    except Exception as e:
        logger.error(f"Metaeuk failed: {str(e)}")
        return None

@celery.task(name="run_multi_pgv_task")
def run_multi_pgv_task(gff_files, seqs, tmp_file, len_thr, id_thr):
    """Celery task to run `multi_pgv`"""
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

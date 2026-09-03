"""Scheduled BLAST database rebuild - compute backend only.

Rebuilds the ships (all/curated) nucleotide BLAST databases, sourmash
signatures and the minimap2 reference FASTA from the current ship DB, and
creates the captain protein DB when missing. This mirrors the monolith's
startup rebuild (app.py initialize_app), which the split frontend no longer
runs - without it the backend's BLAST dirs go stale as ships are added.
"""

from __future__ import annotations

import shutil
import threading
import time

from src.config.celery_config import CELERY_AVAILABLE, celery
from src.config.logging import get_logger
from src.config.settings import IS_DEV

logger = get_logger(__name__)

# A full ships rebuild (makeblastdb + sourmash signatures over the whole
# ship DB) can exceed the 300 s global task_time_limit; give it its own.
REBUILD_TIME_LIMIT = 3600
REBUILD_SOFT_TIME_LIMIT = 3300


def _databases_to_verify() -> dict:
    """Check the core BLAST index files for every managed database."""
    from src.config.settings import BLAST_DB_PATHS
    from src.database.blastdb import blast_db_exists

    return {
        "ships_all": blast_db_exists(BLAST_DB_PATHS["ship"]["all"]["nucl"]),
        "ships_curated": blast_db_exists(BLAST_DB_PATHS["ship"]["curated"]["nucl"]),
        "captains": blast_db_exists(BLAST_DB_PATHS["gene"]["tyr"]["prot"]),
    }


def rebuild_blast_dbs() -> dict:
    """Run the full BLAST DB rebuild. Returns a status dict, never raises."""
    if IS_DEV:
        logger.info("BLAST DB rebuild skipped (dev mode)")
        return {"rebuilt": False, "skipped": "dev mode"}

    if shutil.which("makeblastdb") is None:
        logger.warning("BLAST DB rebuild skipped: makeblastdb not on PATH")
        return {"rebuilt": False, "skipped": "makeblastdb not available"}

    from src.database.blastdb import create_dbs

    started = time.time()
    try:
        create_dbs()
    except Exception as exc:
        logger.error(f"BLAST DB rebuild failed: {exc}")
        logger.exception("Full traceback:")
        return {
            "rebuilt": False,
            "error": str(exc),
            "duration_s": round(time.time() - started, 1),
        }

    databases = _databases_to_verify()
    ok = all(databases.values())
    result = {
        "rebuilt": ok,
        "databases": databases,
        "duration_s": round(time.time() - started, 1),
    }
    if not ok:
        result["error"] = f"Rebuild finished but verification failed: {databases}"
        logger.error(f"BLAST DB rebuild verification failed: {databases}")
    else:
        logger.info(f"BLAST DB rebuild complete in {result['duration_s']}s")
    return result


if CELERY_AVAILABLE and celery:

    @celery.task(
        name="rebuild_blast_dbs",
        time_limit=REBUILD_TIME_LIMIT,
        soft_time_limit=REBUILD_SOFT_TIME_LIMIT,
    )
    def rebuild_blast_dbs_task() -> dict:
        return rebuild_blast_dbs()

    def dispatch_blastdb_rebuild() -> dict:
        async_result = rebuild_blast_dbs_task.delay()
        return {"dispatched": True, "task_id": async_result.id}

else:

    def rebuild_blast_dbs_task() -> dict:
        return rebuild_blast_dbs()

    def dispatch_blastdb_rebuild() -> dict:
        # No Celery worker to hand this off to; run in a background thread
        # like the starfish pipeline does (a rebuild can take minutes).
        threading.Thread(target=rebuild_blast_dbs_task, daemon=True).start()
        return {"dispatched": True, "task_id": None}

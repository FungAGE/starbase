"""Local implementation of the admin submission process/promote workflow.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/admin.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.

promote_submission still opens two separate sessions (get_submissions_session
for the staging row, SubmissionProcessor's own StarbaseSession for the main-DB
insert) even though both now point at the same physical starbase.sqlite file
(see settings.py DB_PATHS unification) -- not a single transaction yet.
SubmissionProcessor is called from several other places (batch CLI tools),
so making it accept an externally-managed session is a separate, wider change,
not bundled into this move.
"""

from src.config.logging import get_logger
from src.database.sql_engine import get_submissions_session
from sqlalchemy import text

logger = get_logger(__name__)


def process_submissions(sub_ids: list) -> list:
    """Run the staging process step (accession + classification) for each
    submission id, sequentially so later ones see earlier ones' results."""
    from src.utils.web_submission_adapter import process_staging_submissions_ordered

    return process_staging_submissions_ordered(sub_ids)


def promote_submission(sub_id: int):
    """
    Load a staging submission and insert it into the main starbase DB via
    perform_database_insertion (SubmissionProcessor).

    Returns (success: bool, accession: str|None, error: str|None).
    """
    from src.utils.web_submission_adapter import (
        perform_database_insertion,
        parse_submission_fasta,
    )

    try:
        with get_submissions_session() as session:
            row = session.execute(
                text("SELECT * FROM submissions WHERE id = :id"), {"id": sub_id}
            ).fetchone()

        if not row:
            return False, None, f"Submission {sub_id} not found"

        row = dict(row._mapping)

        if row.get("processing_status") != "processed":
            return (
                False,
                None,
                f"Submission {sub_id} must be processed before promotion "
                f"(current status: {row.get('processing_status') or 'pending'})",
            )

        if not row.get("seq_contents"):
            return False, None, "No sequence data stored in this submission"

        staging_accession = str(row.get("accession_tag") or "").strip()
        if not staging_accession:
            return (
                False,
                None,
                f"Submission {sub_id} has no accession_tag; run Process first",
            )

        seq_id, sequence = parse_submission_fasta(row["seq_contents"])

        strand = row.get("shipstrand") or "+"
        processed_data = {
            "sequence": sequence,
            "starshipID": seq_id,
            "evidence": row.get("evidence") or "",
            "source": "web_submission_promoted",
            "genus": row.get("genus") or "",
            "species": row.get("species") or "",
            "strain": row.get("strain") or None,
            "contig_id": row.get("hostchr") or "",
            "assembly_accession": row.get("assembly_accession") or None,
            "element_start": row.get("shipstart"),
            "element_end": row.get("shipend"),
            "element_strand": strand,
            "curator": row.get("uploader") or "",
            "curated_status": "curated",
            "notes": row.get("comment") or "",
            "trust_staging": True,
            "staging_accession_tag": staging_accession,
        }

        if row.get("classification_family") or row.get("classification_navis"):
            processed_data["classification"] = {
                "family": row.get("classification_family"),
                "navis": row.get("classification_navis"),
                "haplotype": row.get("classification_haplotype"),
                "source": row.get("classification_source"),
                "closest_match": row.get("closest_match"),
                "confidence": row.get("classification_confidence"),
            }
            if row.get("classification_family"):
                processed_data["ship_family"] = row["classification_family"]
            if row.get("classification_navis"):
                processed_data["ship_navis"] = row["classification_navis"]
            if row.get("classification_haplotype"):
                processed_data["ship_haplotype"] = row["classification_haplotype"]

        result = perform_database_insertion(
            processed_data,
            anno_contents=row.get("anno_contents"),
            anno_filename=row.get("anno_filename"),
            anno_date=None,
            seq_date=0,
        )

        with get_submissions_session() as session:
            session.execute(
                text(
                    """
                    UPDATE submissions
                    SET needs_review = 0, processing_status = 'promoted'
                    WHERE id = :id
                    """
                ),
                {"id": sub_id},
            )
            session.commit()

        logger.info(
            "Promoted submission %s → accession %s", sub_id, result["accession"]
        )
        return True, result["accession"], None

    except Exception as exc:
        logger.error("Promotion failed for submission %s: %s", sub_id, exc)
        return False, None, str(exc)

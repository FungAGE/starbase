"""Local SQLAlchemy implementation of the curation/annotation-review workflow.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/curation.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.
"""

from src.config.logging import get_logger
from src.database.models.schema import (
    Annotation,
    AnnotationHistory,
    BlastpResult,
    GeneFeature,
    HhsearchResult,
    InterproResult,
    RpsblastResult,
)
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

FLAG_LABELS = {
    0: "GREEN",
    1: "YELLOW",
    2: "RED",
    3: "REVIEW_NAME",
    4: "N_A",
    5: "ORANGE",
    7: "UNANNOTATED",
}

_EDITABLE_FIELDS = {
    "annotation",
    "flag",
    "public_notes",
    "private_notes",
    "assigned_to",
}

_RESULT_MODELS = {
    "blastp": BlastpResult,
    "rpsblast": RpsblastResult,
    "hhsearch": HhsearchResult,
    "interpro": InterproResult,
}


def _annotation_to_dict(a: Annotation) -> dict:
    return {
        "id": a.id,
        "sequence": a.sequence,
        "annotation": a.annotation,
        "public_notes": a.public_notes,
        "private_notes": a.private_notes,
        "flag": a.flag,
        "flag_label": FLAG_LABELS.get(a.flag, "UNKNOWN"),
        "assigned_to": a.assigned_to,
        "created_at": a.created_at.isoformat() if a.created_at else None,
        "updated_at": a.updated_at.isoformat() if a.updated_at else None,
    }


def fetch_annotation_queue(
    flag: int = None, assigned_to: str = None, limit: int = 50
) -> list:
    """List annotations for the review queue, oldest-first (stable ordering
    for "next unreviewed" navigation), optionally filtered by flag/assignee."""
    with get_starbase_session() as session:
        query = session.query(Annotation)
        if flag is not None:
            query = query.filter(Annotation.flag == flag)
        if assigned_to:
            query = query.filter(Annotation.assigned_to == assigned_to)
        rows = query.order_by(Annotation.id).limit(limit).all()
        return [_annotation_to_dict(a) for a in rows]


def fetch_annotation(annotation_id: int) -> dict:
    """Single annotation plus its linked gene features, result rows, and
    edit history -- everything the review page needs in one call."""
    with get_starbase_session() as session:
        a = session.query(Annotation).filter_by(id=annotation_id).first()
        if not a:
            raise ValueError(f"Annotation {annotation_id} not found")

        features = (
            session.query(GeneFeature).filter_by(annotation_id=annotation_id).all()
        )
        history = (
            session.query(AnnotationHistory)
            .filter_by(annotation_id=annotation_id)
            .order_by(AnnotationHistory.changed_at)
            .all()
        )
        results = {}
        for tool, model in _RESULT_MODELS.items():
            rows = session.query(model).filter_by(annotation_id=annotation_id).all()
            results[tool] = [
                {
                    "id": r.id,
                    "database": r.database,
                    "result": r.result,
                    "run_date": r.run_date.isoformat() if r.run_date else None,
                    "status": r.status,
                }
                for r in rows
            ]

        return {
            **_annotation_to_dict(a),
            "gene_features": [
                {
                    "id": f.id,
                    "joined_ship_id": f.joined_ship_id,
                    "start": f.start,
                    "stop": f.stop,
                    "type": f.type,
                    "strand": f.strand,
                }
                for f in features
            ],
            "history": [
                {
                    "id": h.id,
                    "changed_by": h.changed_by,
                    "changed_at": h.changed_at.isoformat() if h.changed_at else None,
                    "old_flag": h.old_flag,
                    "new_flag": h.new_flag,
                    "old_annotation": h.old_annotation,
                    "new_annotation": h.new_annotation,
                    "old_public_notes": h.old_public_notes,
                    "new_public_notes": h.new_public_notes,
                    "old_private_notes": h.old_private_notes,
                    "new_private_notes": h.new_private_notes,
                }
                for h in history
            ],
            "results": results,
        }


def update_annotation(annotation_id: int, changes: dict, changed_by: str) -> dict:
    """Apply an edit to an annotation and record it in annotation_history.

    Only annotation/flag/public_notes/private_notes/assigned_to are
    editable -- anything else in `changes` is ignored.
    """
    changes = {k: v for k, v in changes.items() if k in _EDITABLE_FIELDS}
    if not changes:
        raise ValueError("No editable fields in changes")

    with get_starbase_session() as session:
        a = session.query(Annotation).filter_by(id=annotation_id).first()
        if not a:
            raise ValueError(f"Annotation {annotation_id} not found")

        history_kwargs = {"annotation_id": annotation_id, "changed_by": changed_by}
        if "flag" in changes:
            history_kwargs["old_flag"] = a.flag
            history_kwargs["new_flag"] = changes["flag"]
        if "annotation" in changes:
            history_kwargs["old_annotation"] = a.annotation
            history_kwargs["new_annotation"] = changes["annotation"]
        if "public_notes" in changes:
            history_kwargs["old_public_notes"] = a.public_notes
            history_kwargs["new_public_notes"] = changes["public_notes"]
        if "private_notes" in changes:
            history_kwargs["old_private_notes"] = a.private_notes
            history_kwargs["new_private_notes"] = changes["private_notes"]

        for field, value in changes.items():
            setattr(a, field, value)

        from datetime import datetime

        a.updated_at = datetime.utcnow()
        session.add(AnnotationHistory(**history_kwargs))
        session.commit()

        logger.info(
            "Updated annotation %s (%s) by %s", annotation_id, list(changes), changed_by
        )
        return _annotation_to_dict(a)

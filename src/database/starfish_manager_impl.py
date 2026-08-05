"""Local implementation of Starfish pipeline run management.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/starfish.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.

Phase 4 slice 4a: run definition + listing/detail only. Actually launching
the nextflow subprocess (start/cancel/rerun/resume) is slice 4b -- those
endpoints are meaningless without it and are not implemented here yet.

Ported from MAS4starships' starship.tasks/views, fixing a real bug found in
the source: StarfishRunCreateView and run_starfish_pipeline() there compute
samplesheet_path/output_dir independently and disagree (the task
unconditionally overwrites what create wrote). Here there is exactly one
place paths get computed -- at create time -- and nothing else touches them.
"""

import csv
import os

from src.config.logging import get_logger
from src.config.settings import STARFISH_RUNS_DIR
from src.database.models.schema import StarfishElement, StarfishRun, StarfishRunGenome
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

_REQUIRED_GENOME_FIELDS = {"genome_id", "fna_path", "gff3_path"}
_SAMPLESHEET_COLS = ["genomeID", "taxID", "fna", "gff3", "emapper", "cds", "faa"]


def _run_to_dict(r: StarfishRun) -> dict:
    duration = None
    if r.started_at and r.completed_at:
        duration = (r.completed_at - r.started_at).total_seconds()
    return {
        "id": r.id,
        "run_name": r.run_name,
        "description": r.description,
        "status": r.status,
        "created_at": r.created_at.isoformat() if r.created_at else None,
        "started_at": r.started_at.isoformat() if r.started_at else None,
        "completed_at": r.completed_at.isoformat() if r.completed_at else None,
        "created_by": r.created_by,
        "model": r.model,
        "threads": r.threads,
        "missing": r.missing,
        "maxcopy": r.maxcopy,
        "pid_threshold": r.pid_threshold,
        "hsp": r.hsp,
        "flank": r.flank,
        "neighbourhood": r.neighbourhood,
        "samplesheet_path": r.samplesheet_path,
        "output_dir": r.output_dir,
        "log_file": r.log_file,
        "error_message": r.error_message,
        "num_genomes": r.num_genomes,
        "num_elements_found": r.num_elements_found,
        "duration_seconds": duration,
    }


def _genome_to_dict(g: StarfishRunGenome) -> dict:
    return {
        "id": g.id,
        "genome_id": g.genome_id,
        "tax_id": g.tax_id,
        "fna_path": g.fna_path,
        "gff3_path": g.gff3_path,
        "emapper_path": g.emapper_path,
        "cds_path": g.cds_path,
        "faa_path": g.faa_path,
        "num_elements": g.num_elements,
        "status": g.status,
        "error_message": g.error_message,
    }


def _element_to_dict(e: StarfishElement) -> dict:
    return {
        "id": e.id,
        "element_id": e.element_id,
        "genome_id": e.genome_id,
        "contig_id": e.contig_id,
        "start": e.start,
        "end": e.end,
        "strand": e.strand,
        "length": e.end - e.start + 1,
        "family": e.family,
        "navis": e.navis,
        "haplotype": e.haplotype,
        "confidence": e.confidence,
        "notes": e.notes,
        "imported_submission_id": e.imported_submission_id,
    }


def list_runs(status: str = None, limit: int = 50) -> list:
    with get_starbase_session() as session:
        query = session.query(StarfishRun)
        if status:
            query = query.filter(StarfishRun.status == status)
        rows = query.order_by(StarfishRun.created_at.desc()).limit(limit).all()
        return [_run_to_dict(r) for r in rows]


def get_run(run_id: int) -> dict:
    with get_starbase_session() as session:
        r = session.query(StarfishRun).filter_by(id=run_id).first()
        if not r:
            raise ValueError(f"StarfishRun {run_id} not found")
        genomes = session.query(StarfishRunGenome).filter_by(run_id=run_id).all()
        elements = session.query(StarfishElement).filter_by(run_id=run_id).all()
        return {
            **_run_to_dict(r),
            "genomes": [_genome_to_dict(g) for g in genomes],
            "elements": [_element_to_dict(e) for e in elements],
        }


def create_run(
    run_name: str,
    genomes: list,
    description: str = "",
    created_by: str = "curator",
    model: str = "tyr",
    threads: int = 20,
    missing: int = 1,
    maxcopy: int = 5,
    pid_threshold: int = 90,
    hsp: int = 1000,
    flank: int = 6,
    neighbourhood: int = 10000,
) -> dict:
    """Create a StarfishRun + its genome rows and write the samplesheet CSV.

    Does not start the pipeline -- status stays 'pending' (see slice 4b).
    """
    if not run_name or not run_name.strip():
        raise ValueError("run_name is required")
    if not genomes:
        raise ValueError("At least one genome is required")

    clean_genomes = []
    for g in genomes:
        missing_fields = _REQUIRED_GENOME_FIELDS - set(g.keys())
        if missing_fields or not all(g.get(f) for f in _REQUIRED_GENOME_FIELDS):
            raise ValueError(
                f"Genome row missing required field(s): genome_id, fna_path, gff3_path "
                f"(got: {g})"
            )
        clean_genomes.append(g)

    with get_starbase_session() as session:
        existing = session.query(StarfishRun).filter_by(run_name=run_name).first()
        if existing:
            raise ValueError(f"A run named {run_name!r} already exists")

        run = StarfishRun(
            run_name=run_name,
            description=description,
            status="pending",
            created_by=created_by,
            model=model,
            threads=threads,
            missing=missing,
            maxcopy=maxcopy,
            pid_threshold=pid_threshold,
            hsp=hsp,
            flank=flank,
            neighbourhood=neighbourhood,
            num_genomes=len(clean_genomes),
        )
        session.add(run)
        session.flush()  # need run.id for the directory name

        base_dir = os.path.join(STARFISH_RUNS_DIR, f"{run.id}_{run_name}")
        run.samplesheet_path = os.path.join(base_dir, "samplesheet.csv")
        run.output_dir = os.path.join(base_dir, "results")
        run.log_file = os.path.join(base_dir, "starfish.log")

        os.makedirs(base_dir, exist_ok=True)
        with open(run.samplesheet_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=_SAMPLESHEET_COLS)
            writer.writeheader()
            for g in clean_genomes:
                writer.writerow(
                    {
                        "genomeID": g["genome_id"],
                        "taxID": g.get("tax_id", ""),
                        "fna": g["fna_path"],
                        "gff3": g["gff3_path"],
                        "emapper": g.get("emapper_path", ""),
                        "cds": g.get("cds_path", ""),
                        "faa": g.get("faa_path", ""),
                    }
                )

        for g in clean_genomes:
            session.add(
                StarfishRunGenome(
                    run_id=run.id,
                    genome_id=g["genome_id"],
                    tax_id=g.get("tax_id"),
                    fna_path=g["fna_path"],
                    gff3_path=g["gff3_path"],
                    emapper_path=g.get("emapper_path"),
                    cds_path=g.get("cds_path"),
                    faa_path=g.get("faa_path"),
                    status="pending",
                )
            )

        session.commit()
        logger.info(
            "Created starfish run %s (%s) with %d genome(s)",
            run.id,
            run_name,
            len(clean_genomes),
        )
        return _run_to_dict(run)

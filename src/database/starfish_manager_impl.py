"""Local implementation of Starfish pipeline run management.

Runs on backend only (or monolith local-debug). Imported directly by
backend/routers/starfish.py -- no HTTP hop server-side. Mirrors the
sql_manager.py / sql_manager_impl.py split.

Slice 4b adds start/cancel/rerun/resume, which dispatch actual execution
to backend/tasks/starfish.py (Celery task or background thread -- see
there for why not the same inline-call pattern the short BLAST/HMMER
tasks use elsewhere in this codebase).

Ported from MAS4starships' starship.tasks/views, fixing a real bug found in
the source: StarfishRunCreateView and run_starfish_pipeline() there compute
samplesheet_path/output_dir independently and disagree (the task
unconditionally overwrites what create wrote). Here there is exactly one
place paths get computed -- at create time -- and nothing else touches them.
"""

import csv
import os
from datetime import datetime

from src.config.logging import get_logger
from src.config.settings import STARFISH_RUNS_DIR
from src.database.models.schema import (
    StarfishElement,
    StarfishRun,
    StarfishRunGenome,
    Submission,
)
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

        # abspath, not just join: STARFISH_RUNS_DIR inherits DATA_DIR's
        # relativity, and the pipeline subprocess's cwd is set to this same
        # base_dir (see backend/tasks/starfish.py) -- a relative
        # samplesheet_path would then get resolved against that cwd too,
        # silently doubling the path (base_dir/base_dir/samplesheet.csv).
        base_dir = os.path.abspath(os.path.join(STARFISH_RUNS_DIR, f"{run.id}_{run_name}"))
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


def start_run(run_id: int) -> dict:
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            raise ValueError(f"StarfishRun {run_id} not found")
        if run.status != "pending":
            raise ValueError(
                f"StarfishRun {run_id} is {run.status!r}, must be 'pending' to start"
            )
        genome_count = session.query(StarfishRunGenome).filter_by(run_id=run_id).count()
        if genome_count == 0:
            raise ValueError(f"StarfishRun {run_id} has no genomes")

    from backend.tasks.starfish import dispatch_pipeline_run

    dispatch_pipeline_run(run_id, resume=False)
    return get_run(run_id)


def cancel_run(run_id: int) -> dict:
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            raise ValueError(f"StarfishRun {run_id} not found")
        if run.status != "running":
            raise ValueError(
                f"StarfishRun {run_id} is {run.status!r}, can only cancel a 'running' run"
            )
        run.status = "cancelled"
        session.commit()

    from backend.tasks.starfish import cancel_run_process

    cancel_run_process(run_id)
    return get_run(run_id)


def rerun_run(run_id: int) -> dict:
    """Full restart: reset state, delete previously-found elements, and
    dispatch a fresh (non-resumed) pipeline run."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            raise ValueError(f"StarfishRun {run_id} not found")
        if run.status not in ("failed", "completed", "cancelled"):
            raise ValueError(
                f"StarfishRun {run_id} is {run.status!r}, cannot rerun "
                "(must be failed, completed, or cancelled)"
            )
        run.status = "pending"
        run.started_at = None
        run.completed_at = None
        run.celery_task_id = None
        run.process_pid = None
        run.error_message = None
        run.num_elements_found = None
        session.query(StarfishElement).filter_by(run_id=run_id).delete()
        for g in session.query(StarfishRunGenome).filter_by(run_id=run_id).all():
            g.status = "pending"
            g.num_elements = None
        session.commit()

    from backend.tasks.starfish import dispatch_pipeline_run

    dispatch_pipeline_run(run_id, resume=False)
    return get_run(run_id)


def resume_run(run_id: int) -> dict:
    """Resume via nextflow's own -resume/work-dir caching -- no MAS-side
    state diffing, matches the ported behavior exactly."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            raise ValueError(f"StarfishRun {run_id} not found")
        if run.status not in ("failed", "cancelled"):
            raise ValueError(
                f"StarfishRun {run_id} is {run.status!r}, can only resume "
                "'failed' or 'cancelled'"
            )
        run.error_message = None
        session.commit()

    from backend.tasks.starfish import dispatch_pipeline_run

    dispatch_pipeline_run(run_id, resume=True)
    return get_run(run_id)


def _extract_element_sequence(
    fna_path: str, contig_id: str, start: int, end: int, strand: str
) -> str:
    """Extract a starship element's sequence from its source genome FASTA.

    BED coordinates (what start/end come from, see backend/tasks/starfish.py
    parse_bed_file) are 0-based half-open -- matching Python/pyfaidx slicing
    directly, no off-by-one adjustment needed.
    """
    from pyfaidx import Fasta

    from src.utils.seq_utils import revcomp

    if not os.path.exists(fna_path):
        raise ValueError(f"Genome FASTA not found: {fna_path}")

    fasta = Fasta(fna_path)
    if contig_id not in fasta:
        raise ValueError(f"Contig {contig_id!r} not found in {fna_path}")

    seq = str(fasta[contig_id][start:end])
    if strand == "-":
        seq = revcomp(seq)
    return seq


def import_element_to_submission(
    element_id: int, uploader: str = "starfish-pipeline"
) -> dict:
    """Feed one found StarfishElement into the existing submission queue --
    it becomes a regular pending submission, reviewed/processed/promoted
    through the exact same admin workflow as any manually-uploaded ship.

    Fixes MAS4starships' broken import_starfish_elements_to_mas: that wrote
    directly to JoinedShips with field names that don't exist on the real
    model (contigID/elementBegin/elementEnd/ship_family/genome -- would
    TypeError), and never actually populated StarfishElement.sequence in
    the first place (BED files are coordinates only). Here the sequence is
    extracted for real from the source genome FASTA at import time, and
    also backfilled onto the StarfishElement row itself.
    """
    with get_starbase_session() as session:
        element = session.query(StarfishElement).filter_by(id=element_id).first()
        if not element:
            raise ValueError(f"StarfishElement {element_id} not found")
        if element.imported_submission_id:
            raise ValueError(
                f"StarfishElement {element_id} was already imported as "
                f"submission {element.imported_submission_id}"
            )
        genome = (
            session.query(StarfishRunGenome).filter_by(id=element.genome_id).first()
        )
        run = session.query(StarfishRun).filter_by(id=element.run_id).first()
        elem = {
            "element_id": element.element_id,
            "contig_id": element.contig_id,
            "start": element.start,
            "end": element.end,
            "strand": element.strand,
        }
        fna_path = genome.fna_path
        run_name = run.run_name

    sequence = _extract_element_sequence(
        fna_path, elem["contig_id"], elem["start"], elem["end"], elem["strand"]
    )
    if not sequence:
        raise ValueError(
            f"Extracted empty sequence for element {elem['element_id']} "
            f"({elem['contig_id']}:{elem['start']}-{elem['end']})"
        )

    # Not src.pages.submit.insert_submission: that module runs
    # dash.register_page() at import time, which requires a live Dash app --
    # fine inside the monolith, but this code also runs on the backend,
    # a pure FastAPI process with no Dash app at all. Build the row
    # directly against the same Submission model instead.
    seq_contents = f">{elem['element_id']}\n{sequence}\n"
    with get_starbase_session() as session:
        submission = Submission(
            seq_contents=seq_contents,
            seq_filename=f"{elem['element_id']}.fasta",
            seq_date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            uploader=uploader,
            evidence="computational",
            hostchr=elem["contig_id"],
            shipstart=elem["start"],
            shipend=elem["end"],
            shipstrand=elem["strand"],
            comment=f"Found by starfish-nextflow run {run_name!r}, element {elem['element_id']}",
            processing_status="pending",
        )
        session.add(submission)
        session.commit()
        session.refresh(submission)
        submission_id = submission.id

    with get_starbase_session() as session:
        e = session.query(StarfishElement).filter_by(id=element_id).first()
        e.imported_submission_id = submission_id
        e.sequence = sequence
        session.commit()

    logger.info(
        "Imported starfish element %s (run %s) as submission %s",
        elem["element_id"],
        run_name,
        submission_id,
    )
    return {"submission_id": submission_id, "element_id": elem["element_id"]}


def get_run_log(run_id: int) -> str:
    """Raw log file content for the Logs tab. Empty string if not written yet
    (pending run, or file not flushed) rather than an error -- matches
    MAS4starships' "no log file available yet" UI state."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            raise ValueError(f"StarfishRun {run_id} not found")
        log_file = run.log_file

    if not log_file or not os.path.exists(log_file):
        return ""
    with open(log_file) as f:
        return f.read()

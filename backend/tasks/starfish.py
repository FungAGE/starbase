"""Starfish/Nextflow pipeline execution (backend only).

Phase 4 slice 4b: actual pipeline launch, ported from MAS4starships'
starship.tasks.run_starfish_pipeline with real fixes:

- Path handling: MAS4starships computes samplesheet_path/output_dir in TWO
  places (StarfishRunCreateView and this task) that disagree; here paths are
  set once at creation time (starfish_manager_impl.create_run) and never
  recomputed. _results_dir()'s assumption ({output_dir}/{run_name}/) is our
  own clean convention, not verified against the real pipeline -- the actual
  starfish-nextflow source isn't available anywhere on this machine (see
  investigation notes), so this may need adjusting once tested for real.
- verify_starfish_outputs() globbed relative patterns never joined to
  results_dir (a no-op glob against cwd) -- fixed here to actually search
  results_dir.
- Cancel: MAS4starships relies on Celery revoke+terminate, which doesn't
  reliably kill a `bash -c nextflow ...` child process. Here the subprocess
  runs in its own process group (preexec_fn=os.setsid) and its real OS pid
  is stored on the run row, so cancel_run_process() can os.killpg() it
  directly regardless of whether Celery is even in use.
- Dropped MAS4starships' dead "activates starfish environment" comment/code
  that saved-and-noop'd $JAVA_CMD/$JAVA_HOME without ever actually
  activating anything -- conda/module activation is an ops concern for
  whatever invokes this process, not something to fake here.

Cannot be tested against a real pipeline in this environment -- no
`nextflow` binary, no starfish-nextflow checkout, nothing in this repo's
conda envs installs Nextflow/Java. The subprocess/process-group/signal
mechanics ARE tested (see tests), just not against a real `nextflow run`.
"""

import glob
import os
import signal
import subprocess
import threading
from datetime import datetime

from src.config.celery_config import CELERY_AVAILABLE, celery
from src.config.logging import get_logger
from src.config.settings import STARFISH_NEXTFLOW_PATH
from src.database.models.schema import StarfishElement, StarfishRun, StarfishRunGenome
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

_REQUIRED_OUTPUT_PATTERNS = [
    "*.insert.bed",
    "*.insert.stats",
    "*.flank.bed",
    "*.flank.singleDR.stats",
    "*.elements.bed",
    "*.elements.feat",
    "*.elements.named.stats",
]


def build_nextflow_command(run: StarfishRun, resume: bool = False) -> list:
    """Build the nextflow argv list (no shell string -- Popen gets a list,
    no injection risk from run_name/paths)."""
    cmd = ["nextflow", "run", os.path.join(STARFISH_NEXTFLOW_PATH, "main.nf")]
    if resume:
        cmd.append("-resume")
    cmd += [
        "-profile",
        "local",
        "--samplesheet",
        run.samplesheet_path,
        "--run_name",
        run.run_name,
        "--model",
        str(run.model),
        "--threads",
        str(run.threads),
        "--missing",
        str(run.missing),
        "--maxcopy",
        str(run.maxcopy),
        "--pid",
        str(run.pid_threshold),
        "--hsp",
        str(run.hsp),
        "--flank",
        str(run.flank),
        "--neighbourhood",
        str(run.neighbourhood),
        "-w",
        os.path.join(run.output_dir, "work"),
        "--outdir",
        run.output_dir,
    ]
    return cmd


def _results_dir(run: StarfishRun) -> str:
    return os.path.join(run.output_dir, run.run_name)


def verify_starfish_outputs(run: StarfishRun) -> bool:
    results_dir = _results_dir(run)
    for pattern in _REQUIRED_OUTPUT_PATTERNS:
        matches = glob.glob(os.path.join(results_dir, "**", pattern), recursive=True)
        if not matches or not any(os.path.getsize(m) > 0 for m in matches):
            logger.warning(
                "Missing/empty starfish output %s under %s", pattern, results_dir
            )
            return False
    pair_viz = os.path.join(results_dir, "pairViz")
    if not os.path.isdir(pair_viz) or not os.listdir(pair_viz):
        logger.warning("pairViz missing/empty under %s", results_dir)
        return False
    return True


def parse_bed_file(bed_file_path: str, run_id: int, genome_row_id: int) -> int:
    """Parse a BED6 file into StarfishElement rows. element_id (col 4) is
    globally unique (see schema) -- duplicate ids are skipped, not errored,
    since resume/rerun can re-parse the same file. Returns count created."""
    basename = os.path.basename(bed_file_path)
    created = 0
    with get_starbase_session() as session, open(bed_file_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            element_id = parts[3]
            if session.query(StarfishElement).filter_by(element_id=element_id).first():
                continue
            session.add(
                StarfishElement(
                    element_id=element_id,
                    run_id=run_id,
                    genome_id=genome_row_id,
                    contig_id=parts[0],
                    start=int(parts[1]),
                    end=int(parts[2]),
                    strand=parts[5],
                    notes=f"Parsed from {basename}",
                )
            )
            created += 1
        session.commit()
    return created


def parse_starfish_results(run_id: int) -> None:
    """For each genome in the run, find its *.elements.bed under
    regionFinder/ and parse it into StarfishElement rows."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        genomes = session.query(StarfishRunGenome).filter_by(run_id=run_id).all()
        region_dir = os.path.join(_results_dir(run), "regionFinder")
        genome_specs = [(g.id, g.genome_id) for g in genomes]

    total_elements = 0
    for genome_row_id, genome_id in genome_specs:
        matches = [
            p
            for p in glob.glob(os.path.join(region_dir, "*.elements.bed"))
            if genome_id in os.path.basename(p)
        ]
        genome_elements = sum(
            parse_bed_file(bed_file, run_id, genome_row_id) for bed_file in matches
        )
        total_elements += genome_elements
        with get_starbase_session() as session:
            g = session.query(StarfishRunGenome).filter_by(id=genome_row_id).first()
            g.num_elements = genome_elements
            g.status = "completed" if matches else "pending"
            session.commit()

    with get_starbase_session() as session:
        r = session.query(StarfishRun).filter_by(id=run_id).first()
        r.num_elements_found = total_elements
        session.commit()


def _run_starfish_pipeline_impl(
    run_id: int, resume: bool = False, celery_task_id: str = None
) -> None:
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run:
            logger.error("StarfishRun %s not found", run_id)
            return
        if run.status == "running" and not resume:
            logger.warning("StarfishRun %s already running, skipping", run_id)
            return
        if resume and run.status not in ("failed", "cancelled", "running"):
            logger.error(
                "StarfishRun %s not in a resumable state (%s)", run_id, run.status
            )
            return

        run.status = "running"
        run.started_at = datetime.utcnow()
        run.celery_task_id = celery_task_id
        run.error_message = None
        session.commit()
        cmd = build_nextflow_command(run, resume=resume)
        base_dir = os.path.dirname(run.samplesheet_path)
        log_path = run.log_file

    try:
        with open(log_path, "a" if resume else "w") as log_f:
            process = subprocess.Popen(
                cmd,
                cwd=base_dir,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                preexec_fn=os.setsid,
            )
            with get_starbase_session() as session:
                r = session.query(StarfishRun).filter_by(id=run_id).first()
                r.process_pid = process.pid
                session.commit()

            returncode = process.wait()

        with get_starbase_session() as session:
            r = session.query(StarfishRun).filter_by(id=run_id).first()
            if r.status == "cancelled":
                # cancel_run() already finalized this while we were waiting
                logger.info("StarfishRun %s exited after being cancelled", run_id)
                return
            r.process_pid = None
            outputs_ok = returncode == 0 and verify_starfish_outputs(run)
            if outputs_ok:
                r.status = "completed"
                r.completed_at = datetime.utcnow()
                session.commit()
                parse_starfish_results(run_id)
                return
            r.status = "failed"
            r.error_message = (
                "Pipeline completed but expected output files are missing or empty"
                if returncode == 0
                else f"Pipeline failed with return code {returncode}"
            )
            session.commit()

    except Exception as exc:
        logger.error("StarfishRun %s crashed: %s", run_id, exc)
        with get_starbase_session() as session:
            r = session.query(StarfishRun).filter_by(id=run_id).first()
            if r and r.status != "cancelled":
                r.status = "failed"
                r.error_message = str(exc)
                r.process_pid = None
                session.commit()


def cancel_run_process(run_id: int) -> bool:
    """Send SIGTERM to the run's process group. Returns True if a signal was
    sent (doesn't guarantee the process has exited yet -- the running task's
    own wait() loop notices the exit and skips overwriting 'cancelled')."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        if not run or not run.process_pid:
            return False
        pid = run.process_pid

    try:
        pgid = os.getpgid(pid)
        os.killpg(pgid, signal.SIGTERM)
        logger.info("Sent SIGTERM to process group %s (run %s)", pgid, run_id)
        return True
    except ProcessLookupError:
        logger.warning("Process %s for run %s already gone", pid, run_id)
        return False


if CELERY_AVAILABLE and celery:

    @celery.task(name="run_starfish_pipeline", bind=True)
    def run_starfish_pipeline_task(self, run_id, resume=False):
        return _run_starfish_pipeline_impl(
            run_id, resume=resume, celery_task_id=self.request.id
        )

    def dispatch_pipeline_run(run_id: int, resume: bool = False) -> None:
        run_starfish_pipeline_task.delay(run_id, resume=resume)

else:

    def run_starfish_pipeline_task(run_id, resume=False):
        return _run_starfish_pipeline_impl(run_id, resume=resume)

    def dispatch_pipeline_run(run_id: int, resume: bool = False) -> None:
        # No Celery worker to hand this off to -- a pipeline run can take
        # hours, so run it in a background thread rather than blocking the
        # HTTP request that triggered it (unlike the short BLAST/HMMER
        # tasks elsewhere in this codebase, which run inline synchronously).
        threading.Thread(
            target=run_starfish_pipeline_task,
            args=(run_id,),
            kwargs={"resume": resume},
            daemon=True,
        ).start()

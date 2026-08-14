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

from __future__ import annotations

import csv
import glob
import os
import signal
import subprocess
import threading
from datetime import datetime

from src.config.celery_config import CELERY_AVAILABLE, celery
from src.config.logging import get_logger
from src.config.settings import (
    STARFISH_AUX_PATH,
    STARFISH_DB_PATH,
    STARFISH_ENV_PATH,
    STARFISH_NEXTFLOW_PATH,
    STARFISH_NEXTFLOW_PROFILE,
)
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
        STARFISH_NEXTFLOW_PROFILE,
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
        # NOTE: --outdir is accepted (nextflow.config declares params.outdir)
        # but no process actually reads it -- every publishDir in
        # starfish-nextflow is hardcoded to the literal "results/${run_name}"
        # relative to the launch cwd. We rely on that literal instead: cwd is
        # set to base_dir below and run.output_dir is itself base_dir/results,
        # so _results_dir() lines up with the real publish path. Passing
        # --outdir here would be a no-op, so it's omitted.
    ]
    if STARFISH_ENV_PATH:
        cmd += ["--starfish_env", STARFISH_ENV_PATH]
    if STARFISH_AUX_PATH:
        cmd += ["--starfish_aux", STARFISH_AUX_PATH]
    if STARFISH_DB_PATH:
        cmd += ["--starfish_db", STARFISH_DB_PATH]
    return cmd


def _results_dir(run: StarfishRun) -> str:
    """Matches starfish-nextflow's hardcoded `publishDir "results/${run_name}"`
    (relative to the launch cwd, which _run_starfish_pipeline_impl sets to
    dirname(samplesheet_path) == the run's base_dir). run.output_dir is
    base_dir/results (see starfish_manager_impl.create_run), so this join
    reproduces the pipeline's real output path -- it is not itself
    configurable via --outdir (see build_nextflow_command)."""
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
    # PAIR_VIZ has errorStrategy 'ignore' in starfish-nextflow -- it can fail
    # (and its output dir end up missing/empty) on a run that otherwise
    # completed successfully, so this is logged but not treated as failure.
    pair_viz = os.path.join(results_dir, "pairViz")
    if not os.path.isdir(pair_viz) or not os.listdir(pair_viz):
        logger.info(
            "pairViz missing/empty under %s (PAIR_VIZ is allowed to fail)",
            results_dir,
        )
    return True


def _none_if_placeholder(value: str) -> str | None:
    """starfish-nextflow writes "." for absent quality/warnings values."""
    return None if not value or value == "." else value


def _find_genome_for_element_id(element_id: str, genome_specs: list) -> tuple | None:
    """genome_specs: [(genome_row_id, genome_id), ...]. starfish-nextflow
    prefixes every feature/element/region id it generates with
    "{genome_id}_" (default --feature_separator, confirmed against real
    output, e.g. "mp040_s00001") -- match longest genome_id first so one
    genome_id can't shadow another that happens to be its prefix
    (e.g. "mp1" vs "mp10")."""
    for genome_row_id, genome_id in sorted(
        genome_specs, key=lambda spec: len(spec[1]), reverse=True
    ):
        if element_id.startswith(f"{genome_id}_"):
            return genome_row_id, genome_id
    return None


def parse_named_stats_file(stats_file_path: str, run_id: int, genome_specs: list) -> dict:
    """Parse a *.elements.named.stats file into StarfishElement rows -- NOT
    *.elements.bed, which is multiple rows per element: element/boundary/
    cargo-gene features sharing a starshipID (see starfish's main/summarize
    Print_element_bed).

    Real header (tab-separated, confirmed against starfish-nextflow's
    main/summarize Print_named_stats -- note the leading "#" on the first
    column, an actual artifact of the real file, not a comment marker):
    #elementID, elementCaptainID, elementContigID, elementBegin, elementEnd,
    elementLength, elementStrand, emptySiteID, emptyContigID, emptyBegin,
    emptyEnd, emptyLength, emptyStrand, emptySiteSeq, totalFlankAlignment,
    quality, warnings

    Also confirmed against real output: this is NOT one row per element --
    each element gets one row per OTHER genome it was compared against to
    place its empty (pre-insertion) site, with quality 'ref' for the
    canonical/best placement and 'dup' for the rest. elementID/contig/
    start/end/strand are identical across all of one element's rows; only
    the empty-site/quality/warnings columns vary. Rows are grouped by
    elementID and the 'ref' row is preferred as the representative (falls
    back to the first row if no 'ref' row exists for that element).

    element_id is globally unique (see schema) -- duplicate ids are
    skipped, not errored, since resume/rerun can re-parse the same file.
    Returns {genome_row_id: count_created}.
    """
    counts = {genome_row_id: 0 for genome_row_id, _ in genome_specs}
    with open(stats_file_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        reader.fieldnames = [(fn or "").lstrip("#") for fn in reader.fieldnames]
        by_element: dict[str, dict] = {}
        for row in reader:
            element_id = row["elementID"]
            existing = by_element.get(element_id)
            if existing is None or row.get("quality") == "ref":
                by_element[element_id] = row

    with get_starbase_session() as session:
        for element_id, row in by_element.items():
            if session.query(StarfishElement).filter_by(element_id=element_id).first():
                continue
            match = _find_genome_for_element_id(element_id, genome_specs)
            if not match:
                logger.warning(
                    "Could not attribute element %s to any genome in run %s",
                    element_id,
                    run_id,
                )
                continue
            genome_row_id, _ = match
            session.add(
                StarfishElement(
                    element_id=element_id,
                    run_id=run_id,
                    genome_id=genome_row_id,
                    contig_id=row["elementContigID"],
                    start=int(row["elementBegin"]),
                    end=int(row["elementEnd"]),
                    strand=row["elementStrand"],
                    confidence=_none_if_placeholder(row.get("quality")),
                    notes=_none_if_placeholder(row.get("warnings")),
                )
            )
            counts[genome_row_id] += 1
        session.commit()
    return counts


def parse_starfish_results(run_id: int) -> None:
    """Parse the run's single *.elements.named.stats file (one per run, all
    genomes merged) into StarfishElement rows, attributing each to a genome
    via its elementID prefix. Outputs sit flat under _results_dir() --
    starfish-nextflow does not write a per-genome subdir."""
    with get_starbase_session() as session:
        run = session.query(StarfishRun).filter_by(id=run_id).first()
        genomes = session.query(StarfishRunGenome).filter_by(run_id=run_id).all()
        results_dir = _results_dir(run)
        genome_specs = [(g.id, g.genome_id) for g in genomes]

    matches = glob.glob(os.path.join(results_dir, "*.elements.named.stats"))
    counts = {genome_row_id: 0 for genome_row_id, _ in genome_specs}
    if matches:
        for stats_file in matches:
            file_counts = parse_named_stats_file(stats_file, run_id, genome_specs)
            for genome_row_id, n in file_counts.items():
                counts[genome_row_id] += n
    else:
        logger.warning("No *.elements.named.stats found under %s", results_dir)

    total_elements = 0
    for genome_row_id, genome_elements in counts.items():
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
            # NOTE: use r (fetched in this still-open session), not the
            # stale `run` captured back before the subprocess launched --
            # its session is long closed by this point, and touching its
            # attributes here raises a DetachedInstanceError. Confirmed via
            # a real end-to-end pipeline run: this is the first time this
            # code path (a real successful completion) has ever executed.
            outputs_ok = returncode == 0 and verify_starfish_outputs(r)
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

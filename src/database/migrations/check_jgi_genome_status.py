"""
Check restriction/embargo status for JGI MycoCosm genomes in starbase.

Queries genomes where genomeSource == "jgi", looks up each biosample
(MycoCosm portal ID) via the JGI Data Portal mycocosm_file_list endpoint
using f=mycocosm_portal_id (see https://files.jgi.doe.gov/apidoc/), with
MycoCosm releases catalog as fallback when JDP has no match.

API docs: https://files.jgi.doe.gov/apidoc/
Tutorial: https://sites.google.com/lbl.gov/data-portal-help/home/tips_tutorials/api-tutorial

Usage:
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --limit 10
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --biosample Zymps1
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --csv
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --output jgi_status.csv
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --status-filter no
    PYTHONPATH=. python src/database/migrations/check_jgi_genome_status.py --public-only
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
import time
from dataclasses import dataclass, asdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

import requests
from bs4 import BeautifulSoup
from sqlalchemy import text

project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

JGI_MYCOCOSM_SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
MYCOCOSM_RELEASES_URL = "https://mycocosm.jgi.doe.gov/mycocosm/home/releases"
REQUEST_TIMEOUT = 60
LEGACY_POLICY_CUTOFF_YEAR = 2019
RESTRICTION_MONTHS = 24
USER_AGENT = "starbase-jgi-status-check"
JDP_PAGE_SIZE = 50

SEQUENCE_EXCLUDE_HINTS = (
    "gff",
    "gtf",
    "gene catalog",
    "genes_",
    "_genes",
    "protein",
    "transcript",
    "cds_",
    "_cds",
    "est_",
    "chip-seq",
    "chipseq",
    "rnaseq",
    "rna-seq",
    "transcriptome",
    "alignment",
    "metabol",
    "hmm",
    " ipr",
    "sigp",
    "functional annotations",
)
SEQUENCE_INCLUDE_HINTS = (
    "genome assembly",
    "whole genome",
    "assemblyscaffolds",
    "assembly_scaffolds",
    "scaffolds.fasta",
    "contigs.fasta",
    "draft assembly",
)

CSV_FIELDNAMES = [
    "genome_id",
    "ome",
    "biosample",
    "assembly_accession",
    "jgi_status",
    "sequence_status",
    "sequence_file_counts",
    "public_ok",
    "calculated_status",
    "visibility",
    "proposal_date",
    "published",
    "portal_name",
    "api_hits",
    "lookup_state",
    "lookup_source",
    "detail",
]


@dataclass
class GenomeStatusResult:
    genome_id: int
    ome: str
    biosample: str
    assembly_accession: str | None
    jgi_status: str | None
    sequence_status: str | None
    sequence_file_counts: str | None
    public_ok: str
    calculated_status: str
    visibility: str | None
    proposal_date: str | None
    published: bool | None
    portal_name: str | None
    api_hits: int
    lookup_state: str
    lookup_source: str | None = None
    detail: str | None = None


def _parse_date(value: str | None) -> date | None:
    if not value:
        return None
    value = value.strip()
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f%z", "%Y-%m-%dT%H:%M:%S%z", "%Y-%m-%d"):
        try:
            parsed = datetime.strptime(value, fmt)
            return parsed.date()
        except ValueError:
            continue
    if "T" in value:
        try:
            return datetime.fromisoformat(value.replace("Z", "+00:00")).date()
        except ValueError:
            return None
    return None


def _normalize_status(value: str | None) -> str | None:
    if not value:
        return None
    normalized = value.strip()
    if not normalized:
        return None
    lower = normalized.lower()
    if lower == "unrestricted":
        return "Unrestricted"
    if lower == "restricted":
        return "Restricted"
    return normalized


def _is_published(organism: dict[str, Any]) -> bool:
    pubs = organism.get("pubs")
    if isinstance(pubs, list) and len(pubs) > 0:
        return True
    if organism.get("doi"):
        return True
    proposal = organism.get("proposal") or {}
    if isinstance(proposal, dict) and proposal.get("doi"):
        return True
    return False


def _file_metadata_blob(file_entry: dict[str, Any]) -> str:
    metadata = file_entry.get("metadata") or {}
    display_location = (metadata.get("portal") or {}).get("display_location") or []
    parts = [
        file_entry.get("file_name") or "",
        metadata.get("file_format") or "",
        metadata.get("template_name") or "",
        metadata.get("jat_label") or "",
        (metadata.get("final_deliv_project") or {}).get("product_search_category")
        or "",
        " ".join(display_location),
        metadata.get("type") or "",
        metadata.get("subtype") or "",
        metadata.get("content") or "",
    ]
    return " ".join(parts).lower()


def _is_sequence_file(file_entry: dict[str, Any]) -> bool:
    metadata = file_entry.get("metadata") or {}
    display = " ".join(
        (metadata.get("portal") or {}).get("display_location") or []
    ).lower()
    fname = (file_entry.get("file_name") or "").lower()

    if any(hint in display for hint in ("genome assembly", "whole genome assembly")):
        return True
    if any(
        hint in fname
        for hint in ("assemblyscaffolds", "assembly_scaffolds", "genome_assembly")
    ):
        return True

    if any(
        hint in display
        for hint in (
            "annotation",
            "transcriptome",
            "functional annotations",
            "filtered models",
            "all models",
        )
    ):
        return False

    blob = _file_metadata_blob(file_entry)
    if any(hint in blob for hint in SEQUENCE_EXCLUDE_HINTS):
        return False
    if any(hint in blob for hint in SEQUENCE_INCLUDE_HINTS):
        return True

    category = (
        (metadata.get("final_deliv_project") or {}).get("product_search_category") or ""
    ).lower()
    if category == "draft assembly" and metadata.get("file_format") == "fasta":
        return True

    return False


def _analyze_sequence_files(organism: dict[str, Any]) -> tuple[str | None, str | None]:
    unrestricted = 0
    restricted = 0
    mixed = 0

    for file_entry in organism.get("files") or []:
        if not _is_sequence_file(file_entry):
            continue
        metadata = file_entry.get("metadata") or {}
        raw_status = metadata.get("data_utilization_status") or organism.get(
            "data_utilization_status"
        )
        status = _normalize_status(raw_status)
        if status == "Unrestricted":
            unrestricted += 1
        elif status == "Restricted":
            restricted += 1
        elif raw_status and raw_status.strip().lower() == "mixed":
            mixed += 1

    total = unrestricted + restricted + mixed
    if total == 0:
        return None, None

    counts = f"{unrestricted}U/{restricted}R"
    if mixed:
        counts += f"/{mixed}M"

    if restricted == 0 and mixed == 0:
        return "Unrestricted", counts
    if unrestricted == 0 and mixed == 0:
        return "Restricted", counts
    return "Mixed", counts


def _determine_public_ok(
    lookup_state: str,
    jgi_status: str | None,
    sequence_status: str | None,
    lookup_source: str | None,
) -> str:
    if lookup_state != "ok":
        return "review"

    if sequence_status == "Unrestricted":
        return "yes"
    if sequence_status == "Restricted":
        return "no"
    if sequence_status == "Mixed":
        return "review"
    if sequence_status == "Unknown":
        return "review"

    if jgi_status == "Unrestricted":
        return "yes"
    if jgi_status in {"Restricted", "Mixed"} and lookup_source == "mycocosm_catalog":
        return "review"
    if jgi_status == "Restricted":
        return "no"
    return "review"


def _apply_sequence_and_public_fields(
    summary: dict[str, Any],
    organism: dict[str, Any] | None,
    lookup_source: str | None,
) -> None:
    jgi_status = summary.get("jgi_status")
    if lookup_source == "jdp" and organism is not None:
        if jgi_status == "Unrestricted":
            summary["sequence_status"] = "Unrestricted"
            summary["sequence_file_counts"] = None
        else:
            sequence_status, sequence_file_counts = _analyze_sequence_files(organism)
            summary["sequence_status"] = sequence_status or "Unknown"
            summary["sequence_file_counts"] = sequence_file_counts
    else:
        summary["sequence_status"] = jgi_status or "Unknown"
        summary["sequence_file_counts"] = None

    summary["public_ok"] = _determine_public_ok(
        "ok",
        jgi_status,
        summary.get("sequence_status"),
        lookup_source,
    )


def calculate_policy_status(
    proposal_date: date | None,
    published: bool | None,
    api_status: str | None,
) -> str:
    """Apply JGI data policy rules from jgi-data.md."""
    normalized = _normalize_status(api_status)
    if normalized == "Unrestricted":
        return "Unrestricted"
    if normalized == "Restricted":
        return "Restricted / Embargoed"

    if published:
        return "Unrestricted"

    if proposal_date is None:
        return "Unknown"

    if proposal_date.year < LEGACY_POLICY_CUTOFF_YEAR:
        return "Unrestricted (legacy policy)"

    today = datetime.now(timezone.utc).date()
    months_since_release = (today.year - proposal_date.year) * 12 + (
        today.month - proposal_date.month
    )
    if months_since_release >= RESTRICTION_MONTHS:
        return "Unrestricted (2yr expired)"

    return "Restricted / Embargoed"


def _organism_matches_portal(organism: dict[str, Any], biosample: str) -> bool:
    portal_id = organism.get("portal_detail_id")
    if portal_id and portal_id.lower() == biosample.lower():
        return True

    for file_entry in organism.get("files") or []:
        metadata = file_entry.get("metadata") or {}
        mycocosm_portal_id = metadata.get("mycocosm_portal_id")
        if mycocosm_portal_id and mycocosm_portal_id.lower() == biosample.lower():
            return True
    return False


def _summarize_organism(organism: dict[str, Any]) -> dict[str, Any]:
    proposal_date_raw = organism.get("proposal_acceptance_date")
    if not proposal_date_raw:
        for file_entry in organism.get("files") or []:
            metadata = file_entry.get("metadata") or {}
            proposal = metadata.get("proposal") or {}
            proposal_date_raw = proposal.get("date_approved")
            if proposal_date_raw:
                break

    published = _is_published(organism)
    if not published:
        for file_entry in organism.get("files") or []:
            metadata = file_entry.get("metadata") or {}
            if metadata.get("proposal", {}).get("doi"):
                published = True
                break

    return {
        "jgi_status": _normalize_status(organism.get("data_utilization_status")),
        "visibility": organism.get("visibility"),
        "proposal_date_raw": proposal_date_raw,
        "published": published,
        "portal_name": organism.get("name") or organism.get("portal_detail_id"),
    }


def parse_status_from_jdp_payload(
    payload: dict[str, Any], biosample: str
) -> dict[str, Any] | None:
    organisms = payload.get("organisms") or []
    matching = [
        organism
        for organism in organisms
        if isinstance(organism, dict) and _organism_matches_portal(organism, biosample)
    ]

    if not matching and len(organisms) == 1 and isinstance(organisms[0], dict):
        matching = [organisms[0]]

    if not matching:
        return None

    organism = matching[0]
    summary = _summarize_organism(organism)
    proposal_date = _parse_date(summary["proposal_date_raw"])
    summary["proposal_date"] = (
        proposal_date.isoformat() if proposal_date else summary["proposal_date_raw"]
    )
    summary["calculated_status"] = calculate_policy_status(
        proposal_date, summary["published"], summary["jgi_status"]
    )
    _apply_sequence_and_public_fields(summary, organism, "jdp")
    return summary


def search_mycocosm_jdp(
    biosample: str,
    session: requests.Session,
) -> dict[str, Any]:
    params = {
        "q": biosample,
        "f": "mycocosm_portal_id",
        "api_version": "2",
        "x": JDP_PAGE_SIZE,
    }
    response = session.get(
        JGI_MYCOCOSM_SEARCH_URL,
        params=params,
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    payload = response.json()

    page = 1
    while payload.get("next_page") and page < 10:
        page += 1
        response = session.get(
            JGI_MYCOCOSM_SEARCH_URL,
            params={**params, "p": page},
            timeout=REQUEST_TIMEOUT,
        )
        response.raise_for_status()
        next_payload = response.json()
        for organism in next_payload.get("organisms") or []:
            if not isinstance(organism, dict):
                continue
            existing = payload.get("organisms") or []
            if not existing:
                payload["organisms"] = [organism]
                continue
            target = existing[0]
            if organism.get("id") == target.get("id"):
                target.setdefault("files", []).extend(organism.get("files") or [])
            else:
                existing.append(organism)
        if not next_payload.get("next_page"):
            break

    return payload


def _catalog_row_matches_portal(row, biosample: str) -> bool:
    link = row.find("a", href=True)
    if not link:
        return False
    href = link["href"].strip("/")
    return href.lower() == biosample.lower()


def parse_status_from_catalog_html(html: str, biosample: str) -> dict[str, Any] | None:
    soup = BeautifulSoup(html, "html.parser")
    table = soup.find("table", id="catalog")
    if table is None:
        return None

    body = table.find("tbody")
    if body is None:
        return None

    for row in body.find_all("tr"):
        if not _catalog_row_matches_portal(row, biosample):
            continue

        cells = row.find_all("td")
        if len(cells) < 5:
            continue

        portal_link = cells[1].find("a", href=True)
        portal_name = portal_link.get_text(strip=True) if portal_link else None
        utilization = cells[4].get_text(strip=True)
        publication_cell = cells[5]
        published = publication_cell.find("a", href=True) is not None

        jgi_status = _normalize_status(utilization)
        result = {
            "jgi_status": jgi_status,
            "visibility": "public",
            "proposal_date": None,
            "published": published,
            "portal_name": portal_name,
            "calculated_status": calculate_policy_status(None, published, jgi_status),
        }
        _apply_sequence_and_public_fields(result, None, "mycocosm_catalog")
        return result

    return None


def search_mycocosm_catalog(
    biosample: str,
    session: requests.Session,
) -> dict[str, Any] | None:
    response = session.get(
        MYCOCOSM_RELEASES_URL,
        params={"flt": biosample},
        timeout=REQUEST_TIMEOUT,
        headers={"User-Agent": USER_AGENT},
    )
    response.raise_for_status()
    return parse_status_from_catalog_html(response.text, biosample)


def check_genome(
    genome_id: int,
    ome: str,
    biosample: str,
    assembly_accession: str | None,
    session: requests.Session,
) -> GenomeStatusResult:
    api_hits = 0
    parsed: dict[str, Any] | None = None
    lookup_source: str | None = None
    detail: str | None = None

    try:
        payload = search_mycocosm_jdp(biosample, session)
        api_hits = int(payload.get("hits_total") or 0)
        parsed = parse_status_from_jdp_payload(payload, biosample)
        if parsed is not None:
            lookup_source = "jdp"
    except requests.RequestException as exc:
        detail = f"JDP API error: {exc}"

    if parsed is None:
        try:
            parsed = search_mycocosm_catalog(biosample, session)
            if parsed is not None:
                lookup_source = "mycocosm_catalog"
                if detail:
                    detail = f"{detail}; used MycoCosm catalog fallback"
        except requests.RequestException as exc:
            catalog_error = f"MycoCosm catalog error: {exc}"
            detail = f"{detail}; {catalog_error}" if detail else catalog_error

    if parsed is None:
        return GenomeStatusResult(
            genome_id=genome_id,
            ome=ome,
            biosample=biosample,
            assembly_accession=assembly_accession,
            jgi_status=None,
            sequence_status=None,
            sequence_file_counts=None,
            public_ok="review",
            calculated_status="Unknown",
            visibility=None,
            proposal_date=None,
            published=None,
            portal_name=None,
            api_hits=api_hits,
            lookup_state="not_found",
            lookup_source=lookup_source,
            detail=detail or "No match in JDP API or MycoCosm catalog",
        )

    return GenomeStatusResult(
        genome_id=genome_id,
        ome=ome,
        biosample=biosample,
        assembly_accession=assembly_accession,
        jgi_status=parsed.get("jgi_status"),
        sequence_status=parsed.get("sequence_status"),
        sequence_file_counts=parsed.get("sequence_file_counts"),
        public_ok=parsed.get("public_ok", "review"),
        calculated_status=parsed["calculated_status"],
        visibility=parsed.get("visibility"),
        proposal_date=parsed.get("proposal_date"),
        published=parsed.get("published"),
        portal_name=parsed.get("portal_name"),
        api_hits=api_hits,
        lookup_state="ok",
        lookup_source=lookup_source,
        detail=detail,
    )


def fetch_jgi_genomes(limit: int | None = None, biosample: str | None = None):
    sql = """
SELECT id, ome, biosample, assembly_accession
FROM genomes
WHERE lower(genomeSource) = 'jgi'
  AND biosample IS NOT NULL
  AND trim(biosample) != ''
"""
    params: dict[str, Any] = {}
    if biosample:
        sql += "  AND lower(biosample) = lower(:biosample)\n"
        params["biosample"] = biosample
    sql += "ORDER BY ome\n"
    if limit is not None:
        sql += "LIMIT :limit"
        params["limit"] = limit

    with get_starbase_session() as db_session:
        return db_session.execute(text(sql), params).fetchall()


def _quiet_third_party_loggers() -> None:
    for name in ("sqlalchemy.engine", "sqlalchemy.pool", "sqlalchemy.orm"):
        logging.getLogger(name).setLevel(logging.WARNING)


def write_csv_results(
    results: list[GenomeStatusResult],
    output: Any,
) -> None:
    writer = csv.DictWriter(output, fieldnames=CSV_FIELDNAMES)
    writer.writeheader()
    for result in results:
        writer.writerow(asdict(result))


def print_results(results: list[GenomeStatusResult]) -> None:
    for result in results:
        status_bits = [
            f"id={result.genome_id}",
            f"ome={result.ome}",
            f"biosample={result.biosample}",
            f"lookup={result.lookup_state}",
        ]
        if result.lookup_source:
            status_bits.append(f"source={result.lookup_source}")
        if result.jgi_status:
            status_bits.append(f"jgi_status={result.jgi_status}")
        if result.sequence_status:
            status_bits.append(f"sequence_status={result.sequence_status}")
        if result.sequence_file_counts:
            status_bits.append(f"sequence_files={result.sequence_file_counts}")
        status_bits.append(f"public_ok={result.public_ok}")
        status_bits.append(f"calculated={result.calculated_status}")
        if result.visibility:
            status_bits.append(f"visibility={result.visibility}")
        if result.proposal_date:
            status_bits.append(f"proposal_date={result.proposal_date}")
        if result.published is not None:
            status_bits.append(f"published={result.published}")
        if result.portal_name:
            status_bits.append(f"portal={result.portal_name}")
        if result.assembly_accession and result.assembly_accession != result.biosample:
            status_bits.append(f"assembly_accession={result.assembly_accession}")
        if result.detail:
            status_bits.append(f"detail={result.detail}")
        print("  ".join(status_bits))


def summarize_results(results: list[GenomeStatusResult]) -> None:
    portal_counts: dict[str, int] = {}
    public_counts: dict[str, int] = {}
    for result in results:
        key = result.lookup_state
        if result.lookup_state == "ok" and result.jgi_status:
            key = f"ok:{result.jgi_status}"
        portal_counts[key] = portal_counts.get(key, 0) + 1
        public_counts[result.public_ok] = public_counts.get(result.public_ok, 0) + 1

    logger.info("Checked %d JGI genome(s)", len(results))
    for key in sorted(portal_counts):
        logger.info("  portal %s: %d", key, portal_counts[key])
    for key in sorted(public_counts):
        logger.info("  public_ok=%s: %d", key, public_counts[key])


def _result_matches_filter(result: GenomeStatusResult, status_filter: str) -> bool:
    filter_value = status_filter.lower()
    return filter_value in {
        result.lookup_state.lower(),
        (result.jgi_status or "").lower(),
        (result.sequence_status or "").lower(),
        result.public_ok.lower(),
        result.calculated_status.lower(),
        (result.lookup_source or "").lower(),
    }


def run_checks(
    limit: int | None = None,
    biosample: str | None = None,
    csv_output: bool = False,
    output_path: Path | None = None,
    status_filter: str | None = None,
    public_only: bool = False,
    delay: float = 0.25,
) -> list[GenomeStatusResult]:
    if csv_output or output_path is not None:
        _quiet_third_party_loggers()

    rows = fetch_jgi_genomes(limit=limit, biosample=biosample)
    if not rows:
        logger.info("No JGI genomes with biosample found.")
        return []

    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": USER_AGENT})

    results: list[GenomeStatusResult] = []
    for index, row in enumerate(rows):
        genome_id, ome, biosample_value, assembly_accession = row
        result = check_genome(
            genome_id, ome, biosample_value, assembly_accession, session
        )
        results.append(result)

        if not csv_output and output_path is None:
            if (
                not status_filter or _result_matches_filter(result, status_filter)
            ) and (not public_only or result.public_ok != "yes"):
                print_results([result])

        if index + 1 < len(rows) and delay > 0:
            time.sleep(delay)

    if status_filter:
        matched = [r for r in results if _result_matches_filter(r, status_filter)]
    elif public_only:
        matched = [r for r in results if r.public_ok != "yes"]
    else:
        matched = results

    if csv_output or output_path is not None:
        if output_path is not None:
            with output_path.open("w", newline="", encoding="utf-8") as output_file:
                write_csv_results(matched, output_file)
            logger.info("Wrote %d row(s) to %s", len(matched), output_path)
        else:
            write_csv_results(matched, sys.stdout)

    if status_filter or public_only:
        summarize_results(matched)
    else:
        summarize_results(results)

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Check JGI MycoCosm restriction status for starbase genomes "
            "with genomeSource='jgi' using biosample (MycoCosm portal ID)."
        )
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of genomes checked (useful for testing)",
    )
    parser.add_argument(
        "--biosample",
        type=str,
        default=None,
        help="Check a single biosample / MycoCosm portal ID",
    )
    parser.add_argument(
        "--accession",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Write results as CSV to stdout (summary stats still go to stderr)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Write results as CSV to this file",
    )
    parser.add_argument(
        "--status-filter",
        type=str,
        default=None,
        help=(
            "Only output genomes matching lookup_state, jgi_status, sequence_status, "
            "public_ok (yes/no/review), or calculated_status"
        ),
    )
    parser.add_argument(
        "--public-only",
        action="store_true",
        help="Only output genomes where public_ok is not yes (no or review)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.25,
        help="Seconds to wait between API requests (default: 0.25)",
    )
    args = parser.parse_args()

    biosample = args.biosample or args.accession

    run_checks(
        limit=args.limit,
        biosample=biosample,
        csv_output=args.csv,
        output_path=args.output,
        status_filter=args.status_filter,
        public_only=args.public_only,
        delay=args.delay,
    )


if __name__ == "__main__":
    main()

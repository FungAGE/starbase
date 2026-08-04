"""
Link JGI MycoCosm genomes in starbase to NCBI assembly accessions.

Uses the JGI Data Portal (mycocosm_file_list) for portal metadata and NCBI
Entrez assembly records for accession resolution. Also reports existing
starbase NCBI genome rows that appear to represent the same assembly.

Usage:
    PYTHONPATH=. python src/utils/link_jgi_ncbi_genomes.py
    PYTHONPATH=. python src/utils/link_jgi_ncbi_genomes.py --limit 20
    PYTHONPATH=. python src/utils/link_jgi_ncbi_genomes.py --biosample Abobie1
    PYTHONPATH=. python src/utils/link_jgi_ncbi_genomes.py --output jgi_ncbi_links.csv
    PYTHONPATH=. python src/utils/link_jgi_ncbi_genomes.py --duplicates-only
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import requests
from sqlalchemy import text

project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from src.config.logging import get_logger
from src.database.sql_engine import get_starbase_session

logger = get_logger(__name__)

JGI_MYCOCOSM_SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
NCBI_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
REQUEST_TIMEOUT = 60
USER_AGENT = "starbase-jgi-ncbi-link"
JDP_PAGE_SIZE = 50
NCBI_ACCESSION_RE = re.compile(r"GCF_\d+(?:\.\d+)?|GCA_\d+(?:\.\d+)?")

CSV_FIELDNAMES = [
    "genome_id",
    "ome",
    "biosample",
    "db_assembly_accession",
    "portal_name",
    "ncbi_taxon_id",
    "ncbi_assembly_accession",
    "ncbi_assembly_name",
    "ncbi_biosample",
    "match_method",
    "match_confidence",
    "starbase_ncbi_genome_id",
    "starbase_ncbi_ome",
    "lookup_state",
    "detail",
]

JGI_GENOMES_SQL = """
SELECT id, ome, biosample, assembly_accession, taxonomy_id
FROM genomes
WHERE lower(genomeSource) = 'jgi'
  AND biosample IS NOT NULL
  AND trim(biosample) != ''
{filters}
ORDER BY ome
{limit_clause}
"""


@dataclass
class GenomeLinkResult:
    genome_id: int
    ome: str
    biosample: str
    db_assembly_accession: str | None
    portal_name: str | None
    ncbi_taxon_id: str | None
    ncbi_assembly_accession: str | None
    ncbi_assembly_name: str | None
    ncbi_biosample: str | None
    match_method: str | None
    match_confidence: str | None
    starbase_ncbi_genome_id: int | None
    starbase_ncbi_ome: str | None
    lookup_state: str
    detail: str | None = None


def _accession_core(accession: str | None) -> str | None:
    if not accession:
        return None
    match = re.match(r"GC[FAR]_(\d+)", accession.strip(), re.IGNORECASE)
    return match.group(1) if match else None


def _is_ncbi_assembly_accession(value: str | None) -> bool:
    return bool(value and NCBI_ACCESSION_RE.fullmatch(value.strip()))


def _organism_matches_portal(organism: dict[str, Any], biosample: str) -> bool:
    portal_id = organism.get("portal_detail_id")
    if portal_id and portal_id.lower() == biosample.lower():
        return True
    for file_entry in organism.get("files") or []:
        metadata = file_entry.get("metadata") or {}
        portal = metadata.get("mycocosm_portal_id")
        if portal and portal.lower() == biosample.lower():
            return True
    return False


def _select_matching_organism(
    payload: dict[str, Any], biosample: str
) -> dict[str, Any] | None:
    organisms = payload.get("organisms") or []
    matching = [
        organism
        for organism in organisms
        if isinstance(organism, dict) and _organism_matches_portal(organism, biosample)
    ]
    if matching:
        return matching[0]
    if len(organisms) == 1 and isinstance(organisms[0], dict):
        return organisms[0]
    return None


def search_mycocosm_jdp(
    biosample: str,
    session: requests.Session,
) -> dict[str, Any] | None:
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

    return _select_matching_organism(payload, biosample)


def _extract_jdp_accessions(organism: dict[str, Any]) -> set[str]:
    accessions: set[str] = set()
    for file_entry in organism.get("files") or []:
        metadata = file_entry.get("metadata") or {}
        blob = " ".join(
            [
                file_entry.get("file_name") or "",
                str(metadata.get("final_deliv_project") or ""),
                str(metadata.get("gold_data") or ""),
            ]
        )
        accessions.update(NCBI_ACCESSION_RE.findall(blob))
    return accessions


def _extract_jdp_taxon_ids(organism: dict[str, Any]) -> set[str]:
    taxon_ids: set[str] = set()
    for file_entry in organism.get("files") or []:
        metadata = file_entry.get("metadata") or {}
        taxon_id = metadata.get("ncbi_taxon_id")
        if taxon_id is not None:
            taxon_ids.add(str(taxon_id))
    return taxon_ids


def _fetch_ncbi_assemblies_for_taxon(
    taxon_id: str,
    session: requests.Session,
    cache: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    if taxon_id in cache:
        return cache[taxon_id]

    response = session.get(
        NCBI_ESEARCH_URL,
        params={
            "db": "assembly",
            "term": f"txid{taxon_id}[Organism:exp]",
            "retmode": "json",
            "retmax": 200,
        },
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    id_list = response.json().get("esearchresult", {}).get("idlist") or []
    if not id_list:
        cache[taxon_id] = []
        return []

    summary_response = session.get(
        NCBI_ESUMMARY_URL,
        params={
            "db": "assembly",
            "retmode": "json",
            "id": ",".join(id_list),
        },
        timeout=REQUEST_TIMEOUT,
    )
    summary_response.raise_for_status()
    result = summary_response.json().get("result") or {}
    assemblies: list[dict[str, str]] = []
    for uid in result.get("uids") or []:
        record = result.get(uid) or {}
        accession = record.get("assemblyaccession")
        if not accession:
            continue
        properties = record.get("propertylist") or []
        assemblies.append(
            {
                "accession": accession,
                "assembly_name": record.get("assemblyname") or "",
                "organism": record.get("organism") or "",
                "submitter": record.get("submitterorganization") or "",
                "biosample": record.get("biosampleaccn") or "",
                "refseq_category": record.get("refseq_category") or "",
                "properties": "|".join(properties),
            }
        )

    cache[taxon_id] = assemblies
    return assemblies


def _score_ncbi_assembly(
    biosample: str,
    assembly: dict[str, str],
) -> tuple[int, list[str]]:
    score = 0
    reasons: list[str] = []
    assembly_name = assembly.get("assembly_name") or ""
    biosample_lower = biosample.lower()
    name_lower = assembly_name.lower()

    if assembly_name == biosample:
        score += 100
        reasons.append("assembly_name_exact")
    elif name_lower.startswith(biosample_lower) or biosample_lower.startswith(
        name_lower
    ):
        score += 60
        reasons.append("assembly_name_prefix")

    if "joint genome institute" in (assembly.get("submitter") or "").lower():
        score += 25
        reasons.append("jgi_submitter")

    properties = assembly.get("properties") or ""
    if "reference" in properties:
        score += 15
        reasons.append("reference_genome")
    if "latest" in properties:
        score += 5
        reasons.append("latest")

    return score, reasons


def _choose_ncbi_assembly(
    biosample: str,
    assemblies: list[dict[str, str]],
) -> tuple[dict[str, str] | None, str | None, str | None]:
    if not assemblies:
        return None, None, None

    ranked = []
    for assembly in assemblies:
        score, reasons = _score_ncbi_assembly(biosample, assembly)
        ranked.append((score, assembly, reasons))
    ranked.sort(key=lambda item: item[0], reverse=True)

    best_score, best, reasons = ranked[0]
    if best_score <= 0:
        return None, None, None

    if "assembly_name_exact" in reasons:
        confidence = "high"
        method = "ncbi_assembly_name"
    elif "jgi_submitter" in reasons:
        confidence = "medium"
        method = "ncbi_jgi_submitter"
    elif len(assemblies) == 1:
        confidence = "low"
        method = "ncbi_taxon_unique"
    else:
        confidence = "medium"
        method = "ncbi_taxon_ranked"

    return best, method, confidence


def _load_ncbi_genome_index() -> dict[str, list[tuple[int, str, str | None]]]:
    sql = """
    SELECT id, ome, assembly_accession
    FROM genomes
    WHERE lower(genomeSource) = 'ncbi'
      AND assembly_accession IS NOT NULL
      AND trim(assembly_accession) != ''
    """
    index: dict[str, list[tuple[int, str, str | None]]] = {}
    with get_starbase_session() as db_session:
        rows = db_session.execute(text(sql)).fetchall()
    for genome_id, ome, assembly_accession in rows:
        core = _accession_core(assembly_accession)
        if not core:
            continue
        index.setdefault(core, []).append((genome_id, ome, assembly_accession))
    return index


def _find_starbase_ncbi_duplicate(
    accession: str | None,
    ncbi_index: dict[str, list[tuple[int, str, str | None]]],
) -> tuple[int | None, str | None]:
    core = _accession_core(accession)
    if not core:
        return None, None
    matches = ncbi_index.get(core) or []
    if not matches:
        return None, None
    genome_id, ome, _ = matches[0]
    return genome_id, ome


def _resolve_accession(
    biosample: str,
    db_assembly_accession: str | None,
    organism: dict[str, Any] | None,
    ncbi_cache: dict[str, list[dict[str, str]]],
    http_session: requests.Session,
) -> tuple[
    str | None,
    str | None,
    str | None,
    str | None,
    str | None,
    str | None,
    str | None,
    str | None,
]:
    if _is_ncbi_assembly_accession(db_assembly_accession):
        return (
            db_assembly_accession,
            None,
            None,
            "db_assembly_accession",
            "high",
            None,
            None,
            None,
        )

    jdp_accessions: set[str] = set()
    taxon_ids: set[str] = set()
    portal_name = None
    if organism is not None:
        jdp_accessions = _extract_jdp_accessions(organism)
        taxon_ids = _extract_jdp_taxon_ids(organism)
        portal_name = organism.get("name") or organism.get("portal_detail_id")

    if len(jdp_accessions) == 1:
        accession = next(iter(jdp_accessions))
        return (
            accession,
            portal_name,
            next(iter(taxon_ids)) if taxon_ids else None,
            "jdp_metadata",
            "high",
            None,
            None,
            None,
        )
    if len(jdp_accessions) > 1:
        preferred = sorted(
            jdp_accessions,
            key=lambda acc: (0 if acc.startswith("GCF_") else 1, acc),
        )[0]
        return (
            preferred,
            portal_name,
            next(iter(taxon_ids)) if taxon_ids else None,
            "jdp_metadata",
            "medium",
            None,
            None,
            f"Multiple JDP accessions; chose {preferred}",
        )

    if not taxon_ids:
        return (
            None,
            portal_name,
            None,
            None,
            None,
            None,
            None,
            "No NCBI taxon ID in JDP metadata",
        )

    taxon_id = sorted(taxon_ids)[0]
    assemblies = _fetch_ncbi_assemblies_for_taxon(taxon_id, http_session, ncbi_cache)
    best, method, confidence = _choose_ncbi_assembly(biosample, assemblies)
    if best is None:
        return (
            None,
            portal_name,
            taxon_id,
            None,
            None,
            None,
            None,
            f"No NCBI assembly match for taxon {taxon_id}",
        )

    detail = None
    if len(assemblies) > 1 and method == "ncbi_taxon_ranked":
        detail = f"Ranked {len(assemblies)} NCBI assemblies for taxon {taxon_id}"

    return (
        best["accession"],
        portal_name,
        taxon_id,
        method,
        confidence,
        best.get("biosample") or None,
        best.get("assembly_name") or None,
        detail,
    )


def link_genome(
    genome_id: int,
    ome: str,
    biosample: str,
    db_assembly_accession: str | None,
    ncbi_index: dict[str, list[tuple[int, str, str | None]]],
    http_session: requests.Session,
    ncbi_cache: dict[str, list[dict[str, str]]],
) -> GenomeLinkResult:
    organism = None
    detail = None
    try:
        organism = search_mycocosm_jdp(biosample, http_session)
    except requests.RequestException as exc:
        detail = f"JDP API error: {exc}"

    if organism is None and not _is_ncbi_assembly_accession(db_assembly_accession):
        return GenomeLinkResult(
            genome_id=genome_id,
            ome=ome,
            biosample=biosample,
            db_assembly_accession=db_assembly_accession,
            portal_name=None,
            ncbi_taxon_id=None,
            ncbi_assembly_accession=None,
            ncbi_assembly_name=None,
            ncbi_biosample=None,
            match_method=None,
            match_confidence=None,
            starbase_ncbi_genome_id=None,
            starbase_ncbi_ome=None,
            lookup_state="not_found",
            detail=detail or "No JDP portal match",
        )

    (
        ncbi_assembly_accession,
        portal_name,
        ncbi_taxon_id,
        match_method,
        match_confidence,
        ncbi_biosample,
        ncbi_assembly_name,
        match_detail,
    ) = _resolve_accession(
        biosample,
        db_assembly_accession,
        organism,
        ncbi_cache,
        http_session,
    )

    if match_detail and detail:
        detail = f"{detail}; {match_detail}"
    elif match_detail:
        detail = match_detail

    duplicate_id, duplicate_ome = _find_starbase_ncbi_duplicate(
        ncbi_assembly_accession,
        ncbi_index,
    )

    lookup_state = "linked" if ncbi_assembly_accession else "unlinked"
    return GenomeLinkResult(
        genome_id=genome_id,
        ome=ome,
        biosample=biosample,
        db_assembly_accession=db_assembly_accession,
        portal_name=portal_name,
        ncbi_taxon_id=ncbi_taxon_id,
        ncbi_assembly_accession=ncbi_assembly_accession,
        ncbi_assembly_name=ncbi_assembly_name,
        ncbi_biosample=ncbi_biosample,
        match_method=match_method,
        match_confidence=match_confidence,
        starbase_ncbi_genome_id=duplicate_id,
        starbase_ncbi_ome=duplicate_ome,
        lookup_state=lookup_state,
        detail=detail,
    )


def fetch_jgi_genomes(
    limit: int | None = None,
    biosample: str | None = None,
) -> list[tuple[int, str, str, str | None, int | None]]:
    filters = ""
    limit_clause = ""
    params: dict[str, Any] = {}
    if biosample:
        filters = "  AND lower(biosample) = lower(:biosample)\n"
        params["biosample"] = biosample
    if limit is not None:
        limit_clause = "LIMIT :limit"
        params["limit"] = limit

    sql = JGI_GENOMES_SQL.format(filters=filters, limit_clause=limit_clause)
    with get_starbase_session() as db_session:
        return db_session.execute(text(sql), params).fetchall()


def _quiet_third_party_loggers() -> None:
    for name in ("sqlalchemy.engine", "sqlalchemy.pool", "sqlalchemy.orm"):
        logging.getLogger(name).setLevel(logging.WARNING)


def write_csv_results(results: list[GenomeLinkResult], output: Any) -> None:
    writer = csv.DictWriter(output, fieldnames=CSV_FIELDNAMES)
    writer.writeheader()
    for result in results:
        writer.writerow(asdict(result))


def print_results(results: list[GenomeLinkResult]) -> None:
    for result in results:
        bits = [
            f"id={result.genome_id}",
            f"ome={result.ome}",
            f"biosample={result.biosample}",
            f"lookup={result.lookup_state}",
        ]
        if result.ncbi_assembly_accession:
            bits.append(f"ncbi={result.ncbi_assembly_accession}")
        if result.match_method:
            bits.append(f"method={result.match_method}")
        if result.match_confidence:
            bits.append(f"confidence={result.match_confidence}")
        if result.starbase_ncbi_genome_id:
            bits.append(f"duplicate_ncbi_id={result.starbase_ncbi_genome_id}")
        if result.detail:
            bits.append(f"detail={result.detail}")
        print("  ".join(bits))


def summarize_results(results: list[GenomeLinkResult]) -> None:
    lookup_counts: dict[str, int] = {}
    method_counts: dict[str, int] = {}
    duplicate_count = 0
    for result in results:
        lookup_counts[result.lookup_state] = (
            lookup_counts.get(result.lookup_state, 0) + 1
        )
        if result.match_method:
            method_counts[result.match_method] = (
                method_counts.get(result.match_method, 0) + 1
            )
        if result.starbase_ncbi_genome_id:
            duplicate_count += 1

    logger.info("Processed %d JGI genome(s)", len(results))
    for key in sorted(lookup_counts):
        logger.info("  lookup %s: %d", key, lookup_counts[key])
    for key in sorted(method_counts):
        logger.info("  method %s: %d", key, method_counts[key])
    logger.info("  starbase NCBI duplicates found: %d", duplicate_count)


def run_links(
    limit: int | None = None,
    biosample: str | None = None,
    csv_output: bool = False,
    output_path: Path | None = None,
    duplicates_only: bool = False,
    delay: float = 0.25,
) -> list[GenomeLinkResult]:
    if csv_output or output_path is not None:
        _quiet_third_party_loggers()

    rows = fetch_jgi_genomes(limit=limit, biosample=biosample)
    if not rows:
        logger.info("No JGI genomes with biosample found.")
        return []

    ncbi_index = _load_ncbi_genome_index()
    http_session = requests.Session()
    http_session.headers.update(
        {"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    ncbi_cache: dict[str, list[dict[str, str]]] = {}

    results: list[GenomeLinkResult] = []
    for index, row in enumerate(rows):
        genome_id, ome, biosample_value, assembly_accession, _taxonomy_id = row
        result = link_genome(
            genome_id,
            ome,
            biosample_value,
            assembly_accession,
            ncbi_index,
            http_session,
            ncbi_cache,
        )
        results.append(result)

        if not csv_output and output_path is None:
            if not duplicates_only or result.starbase_ncbi_genome_id:
                print_results([result])

        if index + 1 < len(rows) and delay > 0:
            time.sleep(delay)

    matched = results
    if duplicates_only:
        matched = [result for result in results if result.starbase_ncbi_genome_id]

    if csv_output or output_path is not None:
        if output_path is not None:
            with output_path.open("w", newline="", encoding="utf-8") as output_file:
                write_csv_results(matched, output_file)
            logger.info("Wrote %d row(s) to %s", len(matched), output_path)
        else:
            write_csv_results(matched, sys.stdout)

    summarize_results(matched if duplicates_only else results)
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Link JGI MycoCosm genomes to NCBI assembly accessions and report "
            "matching starbase NCBI genome rows."
        )
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--biosample", type=str, default=None)
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Write CSV to stdout",
    )
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument(
        "--duplicates-only",
        action="store_true",
        help="Only output rows with a matching starbase NCBI genome",
    )
    parser.add_argument("--delay", type=float, default=0.25)
    args = parser.parse_args()

    run_links(
        limit=args.limit,
        biosample=args.biosample,
        csv_output=args.csv,
        output_path=args.output,
        duplicates_only=args.duplicates_only,
        delay=args.delay,
    )


if __name__ == "__main__":
    main()

"""
Synteny session API: create and load custom synteny sessions from BLAST results.
"""

import json
import uuid
from flask import Blueprint, jsonify, request

from src.config.cache import cache
from src.config.logging import get_logger
from src.utils.synteny_queries import resolve_accessions_to_ship_ids

logger = get_logger(__name__)

synteny_routes = Blueprint("synteny", __name__, url_prefix="/api/synteny")

SYNTENY_SESSION_PREFIX = "synteny_session:"
SYNTENY_SESSION_TTL = 3600  # 1 hour


def _safe_int(val, default=0):
    """Coerce to int for query_start/query_end/evalue from CSV."""
    if val is None:
        return default
    try:
        return int(float(val))
    except (TypeError, ValueError):
        return default


def _safe_float(val, default=0.0):
    """Coerce to float for evalue."""
    if val is None:
        return default
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def _extract_top_blast_hits_by_coverage(
    blast_df,
    query_length,
    min_query_coverage=0.5,
    max_hits=3,
):
    """
    Return accessions for top BLAST hits with query coverage > min_query_coverage.
    Hits must have annotations (resolved later via resolve_accessions_to_ship_ids).
    Sorted by evalue (best first); returns at most max_hits unique accessions.
    """
    if not blast_df or query_length <= 0:
        return []
    filtered = []
    for row in blast_df:
        qstart = _safe_int(row.get("query_start"))
        qend = _safe_int(row.get("query_end"))
        hit_ids = row.get("hit_IDs")
        if hit_ids is None:
            continue
        span = max(0, qend - qstart)
        coverage = span / query_length if query_length else 0
        if coverage < min_query_coverage:
            continue
        evalue = _safe_float(row.get("evalue"), default=1e10)
        if isinstance(hit_ids, str):
            hit_ids = [hit_ids]
        for acc in hit_ids:
            acc = (acc or "").strip().strip("'\"")
            if acc:
                filtered.append((evalue, acc))
                break  # one accession per row
    # Sort by evalue ascending (best first), then take first max_hits unique accessions
    filtered.sort(key=lambda x: x[0])
    seen = set()
    accessions = []
    for _evalue, acc in filtered:
        if acc not in seen and len(accessions) < max_hits:
            seen.add(acc)
            accessions.append(acc)
    return accessions


def _extract_accessions_from_blast_df(blast_df, max_hits=10):
    """Extract unique accession strings from BLAST hit list (blast_df records)."""
    if not blast_df:
        return []
    accessions = []
    seen = set()
    for row in blast_df[: max_hits * 2]:  # allow some duplicates
        hit_ids = row.get("hit_IDs")
        if hit_ids is None:
            continue
        if isinstance(hit_ids, str):
            hit_ids = [hit_ids]
        for acc in hit_ids:
            acc = (acc or "").strip().strip("'\"")
            if acc and acc not in seen:
                seen.add(acc)
                accessions.append(acc)
        if len(accessions) >= max_hits:
            break
    return accessions[:max_hits]


def _parse_gff3_to_entries(gff_text, contig_label="user_sequence"):
    """
    Parse a GFF3 text string into a list of dicts compatible with the synteny viewer.
    Only CDS/gene/exon feature rows with valid start/end are kept.
    """
    if not gff_text or not isinstance(gff_text, str):
        return []
    entries = []
    feature_types_to_keep = {"CDS", "gene", "exon", "mRNA", "transcript"}
    for line in gff_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        feature_type = parts[2]
        if feature_type not in feature_types_to_keep:
            continue
        try:
            start = int(parts[3])
            end = int(parts[4])
        except (ValueError, IndexError):
            continue
        entries.append(
            {
                "id": None,
                "source": parts[1],
                "type": feature_type,
                "start": start,
                "end": end,
                "score": parts[5] if parts[5] != "." else None,
                "strand": parts[6] if parts[6] in ("+", "-") else "+",
                "phase": parts[7] if parts[7] not in (".", "") else None,
                "attributes": parts[8],
                "ship_id": None,
                "ship_accession_display": contig_label,
            }
        )
    return entries


def _run_metaeuk_and_get_gff(fasta_text, label="user_sequence"):
    """
    Write FASTA to temp file, run metaeuk_easy_predict, return parsed GFF entries.
    Returns [] on any error (metaeuk not installed, failure, etc.).
    """
    import tempfile
    import os

    gff_entries = []
    tmp_fasta = None
    try:
        from src.utils.classification_utils import metaeuk_easy_predict

        # Write FASTA to a temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as fh:
            fh.write(fasta_text)
            tmp_fasta = fh.name

        _codon, _prot, gff_path = metaeuk_easy_predict(tmp_fasta)
        if gff_path and os.path.exists(gff_path):
            with open(gff_path, "r") as fh:
                gff_text = fh.read()
            gff_entries = _parse_gff3_to_entries(gff_text, contig_label=label)
            logger.info(
                f"metaeuk produced {len(gff_entries)} GFF entries for custom track"
            )
        else:
            logger.warning("metaeuk_easy_predict returned no GFF file")
    except Exception as e:
        logger.warning(f"metaeuk annotation failed (non-fatal): {e}")
    finally:
        if tmp_fasta and os.path.exists(tmp_fasta):
            try:
                os.unlink(tmp_fasta)
            except OSError:
                pass
    return gff_entries


def _fasta_to_sequence_and_label(fasta_text):
    """Parse FASTA string; return (sequence, label, length)."""
    if not fasta_text or not isinstance(fasta_text, str):
        return None, "User sequence", 0
    lines = [ln.strip() for ln in fasta_text.strip().splitlines() if ln.strip()]
    if not lines:
        return None, "User sequence", 0
    seq_parts = []
    label = "User sequence"
    for line in lines:
        if line.startswith(">"):
            if seq_parts:
                break
            label = line[1:].strip().split()[0] if line[1:].strip() else "User sequence"
        else:
            seq_parts.append(line)
    sequence = "".join(seq_parts).replace(" ", "")
    return sequence, label, len(sequence)


@synteny_routes.route("/session", methods=["POST"])
def create_synteny_session():
    """
    Create a synteny session from BLAST/classification context.
    Body: { "fasta": str, "blast_df": list of dicts (optional), "classification": dict (optional) }
    Returns: { "session_token": str, "redirect_url": str }
    """
    try:
        data = request.get_json(force=True, silent=True) or {}
        fasta = data.get("fasta") or data.get("sequence")
        blast_df = data.get("blast_df")
        # classification kept for future use (e.g. label)
        # classification = data.get("classification") or {}

        sequence, label, seq_len = _fasta_to_sequence_and_label(fasta or "")
        if not sequence:
            return jsonify({"error": "No valid FASTA sequence provided"}), 400

        preselected_ship_ids = []
        preselected_ships = []  # list of { id, label, ... } for client

        if blast_df and seq_len > 0:
            # Top 3 BLAST hits with query coverage > 50%; resolve_accessions_to_ship_ids keeps only ships with GFF
            accessions = _extract_top_blast_hits_by_coverage(
                blast_df,
                query_length=seq_len,
                min_query_coverage=0.5,
                max_hits=3,
            )
            if accessions:
                resolved = resolve_accessions_to_ship_ids(accessions, max_ships=3)
                for ship_id, info in resolved:
                    preselected_ship_ids.append(ship_id)
                    preselected_ships.append(
                        {
                            "id": ship_id,
                            "label": info.get("label", ""),
                            "ship_accession_display": info.get(
                                "ship_accession_display", ""
                            ),
                        }
                    )

        session_payload = {
            "custom_track": {
                "fasta": fasta,
                "sequence": sequence,
                "label": label,
                "length": seq_len,
            },
            "preselected_ship_ids": preselected_ship_ids,
            "preselected_ships": preselected_ships,
        }

        token = str(uuid.uuid4())
        cache_key = f"{SYNTENY_SESSION_PREFIX}{token}"
        cache.set(
            cache_key,
            json.dumps(session_payload),
            timeout=SYNTENY_SESSION_TTL,
        )

        redirect_url = f"/synteny?session={token}"
        logger.info(
            f"Created synteny session {token} with {len(preselected_ship_ids)} ships"
        )
        return jsonify(
            {
                "session_token": token,
                "redirect_url": redirect_url,
            }
        )
    except Exception as e:
        logger.error(f"Create synteny session failed: {e}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@synteny_routes.route("/session/<token>", methods=["GET"])
def get_synteny_session(token):
    """
    Load a synteny session by token.
    Returns: { "custom_track": {...}, "preselected_ship_ids": [...], "preselected_ships": [...] }
    """
    if not token:
        return jsonify({"error": "Missing session token"}), 400
    cache_key = f"{SYNTENY_SESSION_PREFIX}{token}"
    raw = cache.get(cache_key)
    if not raw:
        return jsonify({"error": "Session not found or expired"}), 404
    try:
        payload = json.loads(raw)
        return jsonify(payload)
    except json.JSONDecodeError:
        return jsonify({"error": "Invalid session data"}), 500

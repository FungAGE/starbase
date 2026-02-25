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

        if blast_df:
            accessions = _extract_accessions_from_blast_df(blast_df, max_hits=10)
            if accessions:
                resolved = resolve_accessions_to_ship_ids(accessions, max_ships=10)
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

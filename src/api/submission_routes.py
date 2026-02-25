"""
Submission prefill API: create and load prefill data for the submission form from BLAST.
"""

import base64
import json
from flask import Blueprint, jsonify, request

from src.config.cache import cache
from src.config.logging import get_logger

logger = get_logger(__name__)

submission_routes = Blueprint("submission", __name__, url_prefix="/api/submission")

SUBMISSION_PREFILL_PREFIX = "submission_prefill:"
SUBMISSION_PREFILL_TTL = 3600  # 1 hour


def _fasta_to_base64_contents(fasta_text, filename="query_sequence.fasta"):
    """Turn FASTA string into the format Dash Upload expects: data URL with base64."""
    if not fasta_text or not isinstance(fasta_text, str):
        return None, filename
    encoded = base64.b64encode(fasta_text.strip().encode("utf-8")).decode("ascii")
    return f"data:application/octet-stream;base64,{encoded}", filename


@submission_routes.route("/prefill", methods=["GET"])
def get_submission_prefill():
    """
    Get submission prefill data by token.
    Query: ?token=<uuid>
    Returns: { "seq_contents": str (base64 data URL), "seq_filename": str, ...metadata }
    """
    token = request.args.get("token")
    if not token:
        return jsonify({"error": "Missing token"}), 400
    cache_key = f"{SUBMISSION_PREFILL_PREFIX}{token}"
    raw = cache.get(cache_key)
    if not raw:
        return jsonify({"error": "Prefill not found or expired"}), 404
    try:
        payload = json.loads(raw)
        return jsonify(payload)
    except json.JSONDecodeError:
        return jsonify({"error": "Invalid prefill data"}), 500

from flask import Blueprint, jsonify, request
from src.config.limiter import limiter
from src.components.data import (
    create_ship_accession_modal_data,
    create_accession_modal_data,
)

from src.config.logging import get_logger

logger = get_logger(__name__)

# Create Blueprints for different route groups
blast_routes = Blueprint("blast", __name__, url_prefix="/api/blast")
accession_routes = Blueprint("accession", __name__, url_prefix="/api/accession")
error_handlers = Blueprint("errors", __name__)


@blast_routes.route("/blast-submit", methods=["POST"])
@limiter.limit("10 per hour")
def check_blast_limit():
    """Check BLAST submission limit."""
    remote_addr = request.remote_addr
    logger.info(f"BLAST submission from IP: {remote_addr}")
    return jsonify({"allowed": True})


@blast_routes.route("/export-for-submission", methods=["POST"])
def export_for_submission():
    """
    Create a submission prefill from BLAST result and return redirect URL.
    Body: { "fasta": str, "seq_filename": str (optional), "classification": dict (optional) }
    Returns: { "token": str, "redirect_url": str }
    """
    import json
    import uuid
    from src.config.cache import cache

    SUBMISSION_PREFILL_PREFIX = "submission_prefill:"
    SUBMISSION_PREFILL_TTL = 3600

    def _fasta_to_base64_contents(fasta_text, filename="query_sequence.fasta"):
        if not fasta_text or not isinstance(fasta_text, str):
            return None, filename
        import base64

        encoded = base64.b64encode(fasta_text.strip().encode("utf-8")).decode("ascii")
        return f"data:application/octet-stream;base64,{encoded}", filename

    try:
        data = request.get_json(force=True, silent=True) or {}
        fasta = data.get("fasta") or data.get("sequence")
        seq_filename = data.get("seq_filename") or "query_sequence.fasta"
        classification = data.get("classification") or {}

        seq_contents, seq_filename = _fasta_to_base64_contents(
            fasta or "", seq_filename
        )
        if not seq_contents:
            return jsonify({"error": "No valid FASTA sequence provided"}), 400

        prefill_payload = {
            "seq_contents": seq_contents,
            "seq_filename": seq_filename,
            "classification": classification,
        }

        token = str(uuid.uuid4())
        cache_key = f"{SUBMISSION_PREFILL_PREFIX}{token}"
        cache.set(
            cache_key, json.dumps(prefill_payload), timeout=SUBMISSION_PREFILL_TTL
        )

        redirect_url = f"/submit?token={token}"
        logger.info(f"Created submission prefill {token}")
        return jsonify({"token": token, "redirect_url": redirect_url})
    except Exception as e:
        logger.error(f"Export for submission failed: {e}", exc_info=True)
        return jsonify({"error": str(e)}), 500


@accession_routes.route("/accession_details/<accession_id>", methods=["GET"])
def get_accession_details(accession_id):
    """Get details for a specific ship accession."""
    try:
        if accession_id.startswith("SSA"):
            modal_data = create_accession_modal_data(accession_id)
        if accession_id.startswith("SSB"):
            modal_data = create_ship_accession_modal_data(accession_id)
        else:
            modal_data = {"error": "Invalid accession ID"}
    except Exception as e:
        logger.error(f"Error in get_accession_details: {str(e)}")
        modal_data = {"error": str(e)}

    return jsonify(modal_data)


@error_handlers.app_errorhandler(500)
def handle_500(e):
    """Handle internal server errors."""
    logger.error(f"Internal Server Error: {str(e)}")
    return jsonify(
        {
            "error": "Internal Server Error",
            "message": "The server encountered an error. Please try again later.",
        }
    ), 500


@error_handlers.app_errorhandler(429)
def handle_429(e):
    """Handle too many requests errors."""
    return jsonify(
        {
            "error": "Too Many Requests",
            "message": "Please wait before making more requests.",
        }
    ), 429

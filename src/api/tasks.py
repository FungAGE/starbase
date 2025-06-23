"""
API endpoints for task management.
Provides endpoints for submitting and checking status of async tasks.
"""

from flask import Blueprint, request, jsonify
from src.tasks.task_manager import TaskManager
from src.config.logging import get_logger

logger = get_logger(__name__)

tasks_bp = Blueprint("tasks", __name__, url_prefix="/api/tasks")


@tasks_bp.route("/status/<task_id>", methods=["GET"])
def get_task_status(task_id):
    """Get the status of a task."""
    try:
        status = TaskManager.get_task_status(task_id)
        return jsonify(status)
    except Exception as e:
        logger.error(f"Error getting task status for {task_id}: {e}")
        return jsonify({"status": "ERROR", "error": str(e)}), 500


@tasks_bp.route("/cleanup", methods=["POST"])
def cleanup_old_tasks():
    """Clean up old task records."""
    try:
        max_age = request.json.get("max_age", 3600) if request.json else 3600
        TaskManager.cleanup_old_tasks(max_age)
        return jsonify({"status": "success", "message": "Old tasks cleaned up"})
    except Exception as e:
        logger.error(f"Error cleaning up tasks: {e}")
        return jsonify({"status": "error", "error": str(e)}), 500


@tasks_bp.route("/list", methods=["GET"])
def list_tasks():
    """List all current tasks (for debugging)."""
    try:
        from src.tasks.task_manager import task_status_store

        return jsonify(
            {"tasks": list(task_status_store.keys()), "count": len(task_status_store)}
        )
    except Exception as e:
        logger.error(f"Error listing tasks: {e}")
        return jsonify({"status": "error", "error": str(e)}), 500

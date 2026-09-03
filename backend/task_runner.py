"""Synchronous task execution for backend routers.

Celery tasks are run inline via .apply().get(timeout=...) — the backend
process is the worker, so local eager execution is equivalent to a round
trip through the broker. Plain functions (e.g. in dev/tests without
celery) are called directly.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Tuple

DEFAULT_TIMEOUT = 360


def run_task_sync(
    fn: Callable[..., Any],
    args: Tuple = (),
    kwargs: Optional[Dict[str, Any]] = None,
    *,
    timeout: float = DEFAULT_TIMEOUT,
) -> Any:
    """Run fn synchronously and return its result."""
    kwargs = kwargs or {}
    if hasattr(fn, "apply"):
        return fn.apply(args=args, kwargs=kwargs).get(timeout=timeout)
    return fn(*args, **kwargs)

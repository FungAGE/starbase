"""Minimal Flask app used only to host Flask-Caching on the FastAPI backend."""

from __future__ import annotations

import os
from contextlib import contextmanager

from flask import Flask, has_app_context

from src.config.cache import cache

flask_app = Flask(__name__)
flask_app.config.update(
    CACHE_TYPE=os.getenv("BACKEND_CACHE_TYPE", "SimpleCache"),
    CACHE_DEFAULT_TIMEOUT=int(os.getenv("BACKEND_CACHE_DEFAULT_TIMEOUT", "300")),
)
cache.init_app(flask_app)


@contextmanager
def cache_app_context():
    """Push Flask app context when running outside the Dash/Flask frontend."""
    if has_app_context():
        yield
        return
    with flask_app.app_context():
        yield

"""Shared FastAPI dependencies (auth, Flask cache context for sql_manager)."""

from __future__ import annotations

import os
from collections.abc import Generator

from fastapi import Depends, Header, HTTPException, status
from flask import Flask

from src.config.cache import cache


def _init_flask_for_cache() -> Flask:
    """sql_manager + smart_cache expect Flask-Caching initialized on a Flask app."""
    flask_app = Flask(__name__)
    flask_app.config.update(
        CACHE_TYPE=os.getenv("BACKEND_CACHE_TYPE", "SimpleCache"),
        CACHE_DEFAULT_TIMEOUT=int(os.getenv("BACKEND_CACHE_DEFAULT_TIMEOUT", "300")),
    )
    cache.init_app(flask_app)
    return flask_app


flask_app = _init_flask_for_cache()


def flask_cache_context() -> Generator[None, None, None]:
    with flask_app.app_context():
        yield


def verify_backend_api_key(
    x_api_key: str | None = Header(default=None, alias="X-API-Key"),
    authorization: str | None = Header(default=None),
) -> None:
    expected = os.getenv("BACKEND_API_KEY")
    if not expected:
        return
    token = x_api_key
    if not token and authorization and authorization.lower().startswith("bearer "):
        token = authorization[7:].strip()
    if not token or token != expected:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or missing API key",
        )


RequireApiKey = Depends(verify_backend_api_key)
CacheContext = Depends(flask_cache_context)

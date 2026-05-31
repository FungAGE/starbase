"""FastAPI entrypoint for starbase compute backend."""

from __future__ import annotations

import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path

# Repo root must be on path (same layout as Flask app)
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from fastapi import Depends, FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.dependencies import verify_backend_api_key
from backend.routers import blast, data
from src.config.logging import get_logger
from src.database.migrations import create_database_indexes, run_alembic_migrations

logger = get_logger(__name__)


@asynccontextmanager
async def lifespan(_app: FastAPI):
    """Run DB migrations and create indexes on startup."""
    try:
        run_alembic_migrations()
        create_database_indexes()
    except Exception as exc:
        logger.error("Backend startup DB init failed: %s", exc)
        raise
    yield


def create_app() -> FastAPI:
    app = FastAPI(
        title="starbase backend",
        description="SQL + BLAST/HMMER compute API for starbase",
        lifespan=lifespan,
    )
    origins = os.getenv("BACKEND_CORS_ORIGINS", "").strip()
    if origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=[o.strip() for o in origins.split(",") if o.strip()],
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )

    @app.get("/health")
    def health() -> dict[str, str]:
        return {"status": "ok"}

    @app.get("/api/health", dependencies=[Depends(verify_backend_api_key)])
    def api_health() -> dict[str, str]:
        return {"status": "ok"}

    app.include_router(data.router)
    app.include_router(blast.router)
    return app


app = create_app()

"""FastAPI entrypoint for starbase compute backend."""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Repo root must be on path (same layout as Flask app)
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from fastapi import Depends, FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.dependencies import verify_backend_api_key
from backend.routers import admin, blast, curation, data, starfish
from src.config.logging import get_logger

logger = get_logger(__name__)


def create_app() -> FastAPI:
    app = FastAPI(
        title="starbase backend",
        description="SQL + BLAST/HMMER compute API for starbase",
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
    app.include_router(admin.router)
    app.include_router(curation.router)
    app.include_router(starfish.router)
    return app


app = create_app()

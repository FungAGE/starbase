"""Shared FastAPI dependencies (auth)."""

from __future__ import annotations

import os
from typing import Optional

from fastapi import Depends, Header, HTTPException, status


def verify_backend_api_key(
    x_api_key: Optional[str] = Header(default=None, alias="X-API-Key"),
    authorization: Optional[str] = Header(default=None),
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

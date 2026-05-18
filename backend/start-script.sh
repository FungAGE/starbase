#!/bin/bash
# Run FastAPI backend (ASGI). Build image from repo root: docker build -f backend/Dockerfile .
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT"

DEV_MODE="${BACKEND_DEV_MODE:-false}"
if [[ -f .env ]]; then
	_val="$(grep -E '^DEV_MODE=' .env 2>/dev/null | cut -d'=' -f2- || true)"
	if [[ -n "${_val}" ]]; then
		DEV_MODE="$_val"
	fi
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
	--dev)
		DEV_MODE=true
		shift
		;;
	*)
		echo "Unknown option: $1" >&2
		echo "Usage: $0 [--dev]" >&2
		exit 1
		;;
	esac
done

PORT="${BACKEND_PORT:-8001}"

if [[ "$DEV_MODE" == "true" ]]; then
	export ENVIRONMENT="${ENVIRONMENT:-development}"
	exec uvicorn backend.main:app \
		--host=0.0.0.0 \
		--port="$PORT" \
		--reload \
		--log-level=debug \
		--proxy-headers \
		--forwarded-allow-ips='*' \
		--timeout-keep-alive=5 \
		--timeout-graceful-shutdown=60
else
	export ENVIRONMENT="${ENVIRONMENT:-production}"
	export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
	export SQLALCHEMY_WARN_20="${SQLALCHEMY_WARN_20:-1}"
	export SQLALCHEMY_SILENCE_UBER_WARNING="${SQLALCHEMY_SILENCE_UBER_WARNING:-1}"
	WORKERS="${UVICORN_WORKERS:-4}"
	exec uvicorn backend.main:app \
		--host=0.0.0.0 \
		--port="$PORT" \
		--workers="$WORKERS" \
		--log-level=warning \
		--proxy-headers \
		--forwarded-allow-ips='*' \
		--timeout-keep-alive=5 \
		--timeout-graceful-shutdown=30 \
		--limit-max-requests=1000 \
		--no-access-log
fi

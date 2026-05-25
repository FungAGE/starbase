#!/usr/bin/env bash
# Start the compute backend in Docker and print how to run the frontend on the host.
#
# Phase 1 local test (same machine):
#   ./dev/scripts/setup-split-dev.sh backend   # terminal 1
#   ./dev/scripts/setup-split-dev.sh frontend  # terminal 2
#
# Phase 2 (frontend on laptop, backend on institute box via Tailscale):
#   On institute machine:  ./dev/scripts/setup-split-dev.sh backend
#   On laptop:             BACKEND_API_URL=http://100.x.y.z:8001 ./dev/scripts/setup-split-dev.sh frontend
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

MODE="${1:-}"

ensure_env() {
	if [[ -f .env ]]; then
		return
	fi
	echo "Creating .env from dev/config/env.template..."
	cp dev/config/env.template .env
	# Sensible split-dev defaults
	if ! grep -q '^BACKEND_API_KEY=change-me' .env 2>/dev/null; then
		:
	else
		_secret="dev-$(openssl rand -hex 16 2>/dev/null || date +%s)"
		sed -i "s/^BACKEND_API_KEY=.*/BACKEND_API_KEY=${_secret}/" .env
	fi
	sed -i 's/^BACKEND_API_URL=$/BACKEND_API_URL=http:\/\/localhost:8001/' .env
	sed -i "s|^DATA_DIR=$|DATA_DIR=${REPO_ROOT}/src/database/db|" .env
	echo "Edit .env if needed (BACKEND_API_KEY, DATA_DIR)."
}

load_env() {
	if [[ -f .env ]]; then
		set -a
		# shellcheck disable=SC1091
		source .env
		set +a
	fi
}

start_backend() {
	ensure_env
	load_env

	if [[ ! -d "${DATA_DIR:-$REPO_ROOT/src/database/db}" ]]; then
		echo "ERROR: DATA_DIR does not exist. Set DATA_DIR in .env to your database directory."
		exit 1
	fi

	echo "Starting compute backend (Redis + FastAPI + Celery)..."
	docker compose -f backend/docker-compose.dev.yaml up -d --build

	echo ""
	echo "Waiting for backend health..."
	for _ in $(seq 1 30); do
		if curl -sf "http://localhost:${BACKEND_PORT:-8001}/health" >/dev/null 2>&1; then
			echo "Backend is healthy at http://localhost:${BACKEND_PORT:-8001}"
			break
		fi
		sleep 2
	done

	echo ""
	echo "Backend logs:  docker compose -f backend/docker-compose.dev.yaml logs -f backend"
	echo "Stop backend:  docker compose -f backend/docker-compose.dev.yaml down"
}

start_frontend() {
	ensure_env
	load_env

	export BACKEND_API_URL="${BACKEND_API_URL:-http://localhost:8001}"
	export BACKEND_API_KEY="${BACKEND_API_KEY:-}"
	export DEV_MODE="${DEV_MODE:-true}"
	export ENVIRONMENT="${ENVIRONMENT:-development}"

	if [[ -z "$BACKEND_API_KEY" ]]; then
		echo "WARNING: BACKEND_API_KEY is empty — set it in .env to match the backend container."
	fi

	if ! curl -sf "${BACKEND_API_URL}/health" >/dev/null 2>&1; then
		echo "WARNING: Backend not reachable at ${BACKEND_API_URL}/health"
		echo "         Start it first: $0 backend"
	fi

	echo "Starting frontend on http://localhost:8000"
	echo "  BACKEND_API_URL=${BACKEND_API_URL}"
	echo ""

	if command -v conda >/dev/null 2>&1 && conda env list | grep -q '^starbase '; then
		exec conda run --no-capture-output -n starbase ./start-script.sh --dev
	fi

	exec ./start-script.sh --dev
}

case "$MODE" in
backend)
	start_backend
	;;
frontend)
	start_frontend
	;;
"")
	echo "Usage: $0 {backend|frontend}"
	echo ""
	echo "  backend   Build and run compute stack in Docker (port 8001)"
	echo "  frontend  Run Dash UI on host, talking to BACKEND_API_URL"
	exit 1
	;;
*)
	echo "Unknown mode: $MODE (use backend or frontend)" >&2
	exit 1
	;;
esac

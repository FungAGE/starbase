#!/bin/bash

set -o pipefail

# Load runtime configuration before starting Tailscale or the application.
if [ -n "${ENV_FILE:-}" ]; then
    _env_file="$ENV_FILE"
elif [ -f "$HOME/src/database/db/.env" ]; then
    _env_file="$HOME/src/database/db/.env"
elif [ -f .env ]; then
    _env_file=.env
else
    _env_file=""
fi

if [ -n "$_env_file" ]; then
    echo "Loading runtime configuration from $_env_file"
    set -a
    # shellcheck disable=SC1090
    source "$_env_file"
    set +a
fi

# Tailscale only when TS_AUTHKEY is provided
if [ -n "$TS_AUTHKEY" ]; then
    echo "Starting tailscaled (userspace networking)..."
    tailscaled \
        --state=mem: \
        --socket="$HOME/.tailscale/tailscaled.sock" \
        --tun=userspace-networking \
        --socks5-server=localhost:1055 \
        --outbound-http-proxy-listen=localhost:1055 \
        >"$HOME/.tailscale/tailscaled.log" 2>&1 &

    for _ in $(seq 1 30); do
        [ -S "$HOME/.tailscale/tailscaled.sock" ] && break
        sleep 1
    done

    if [ ! -S "$HOME/.tailscale/tailscaled.sock" ]; then
        echo "ERROR: tailscaled did not create its control socket." >&2
        tail -100 "$HOME/.tailscale/tailscaled.log" >&2
        exit 1
    fi

    if ! tailscale --socket="$HOME/.tailscale/tailscaled.sock" up \
        --auth-key="$TS_AUTHKEY" \
        --hostname="${TS_HOSTNAME:-starbase-frontend}" \
        --accept-dns=false \
        --timeout=30s; then
        echo "ERROR: Tailscale authentication failed." >&2
        tail -100 "$HOME/.tailscale/tailscaled.log" >&2
        exit 1
    fi

    # Scoped to backend_client.py only -- NOT HTTP_PROXY/HTTPS_PROXY globally,
    # since the tailnet netstack only routes tailnet addresses (no exit node
    # configured), so a global proxy would break NCBI/IPStack/SMTP calls to
    # the public internet.
    export TAILSCALE_PROXY="http://localhost:1055"
    echo "Tailscale up. Backend calls will route via tailnet proxy on :1055."

    # Early backend reachability probe
    if [ -n "${BACKEND_API_URL:-}" ]; then
        if _health=$(curl -sf --max-time 15 -x "http://localhost:1055" "${BACKEND_API_URL}/health" 2>&1); then
            echo "Backend reachable at ${BACKEND_API_URL}: ${_health}"
        else
            echo "WARNING: backend NOT reachable at ${BACKEND_API_URL}/health via tailnet." >&2
            _ts_sock="$HOME/.tailscale/tailscaled.sock"
            _backend_host="${BACKEND_API_URL#*://}"; _backend_host="${_backend_host%%/*}"
            _backend_peer="${_backend_host%%:*}"
            echo "         Pod tailnet status (look for the backend peer + online/offline):" >&2
            tailscale --socket="$_ts_sock" status 2>&1 | sed 's/^/           /' >&2
            echo "         Pinging ${_backend_peer} (2 tries; 'direct'/'relay' = path works):" >&2
            tailscale --socket="$_ts_sock" ping --c 2 "$_backend_peer" 2>&1 | sed 's/^/           /' >&2
            echo "         netcheck (DERP relay reachability from this pod):" >&2
            timeout 45 tailscale --socket="$_ts_sock" netcheck 2>&1 | tail -25 | sed 's/^/           /' >&2
            echo "         tailscaled.log (last 20 lines):" >&2
            tail -20 "$HOME/.tailscale/tailscaled.log" 2>&1 | sed 's/^/           /' >&2
            echo "         Peer offline/stale -> check backend machine: sudo tailscale status" >&2
            echo "         Peer online + ping OK -> port 8001 closed: docker compose ps on backend" >&2
        fi
    fi
fi


# Parse command line arguments
while [ $# -gt 0 ]; do
    case $1 in
        --dev)
            DEV_MODE=true
            shift # Remove --dev from processing
            ;;
        *)
            # Unknown option
            echo "Unknown option: $1"
            echo "Usage: $0 [--dev]"
            exit 1
            ;;
    esac
done

# Set environment variables based on mode
export PYTHONPATH=$(pwd)
if [ "$DEV_MODE" = "true" ]; then
    export ENVIRONMENT="development"
else
    export ENVIRONMENT="production"
    export PYTHONUNBUFFERED=1
    export SQLALCHEMY_WARN_20=1
    export SQLALCHEMY_SILENCE_UBER_WARNING=1
fi

# Check if --dev flag was provided
if [ "$DEV_MODE" = "true" ]; then
    # Development mode with reload
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --reload \
        --log-level=debug \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --access-log \
        app:server
else
    # Production mode with workers
    # we are provided 4cpu and 8gb of memory within a single pod
    # we can use 4 workers to maximize the utilization of the resources

    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --workers=4 \
        --log-level=warning \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=30 \
        --limit-max-requests=1000 \
        --no-access-log \
        app:server
fi

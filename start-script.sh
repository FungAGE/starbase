#!/bin/bash

# Load .env so uvicorn workers inherit ADMIN_TOKEN, DEV_MODE, etc.
if [ -f .env ]; then
    set -a
    # shellcheck disable=SC1091
    source .env
    set +a
fi

if [ -z "$DEV_MODE" ]; then
    DEV_MODE=false
fi

# Tailscale only when TS_AUTHKEY is provided
if [ -n "$TS_AUTHKEY" ]; then
    echo "Starting tailscaled (userspace networking)..."
    tailscaled \
        --state="$HOME/.tailscale/tailscaled.state" \
        --socket="$HOME/.tailscale/tailscaled.sock" \
        --tun=userspace-networking \
        --socks5-server=localhost:1055 \
        --outbound-http-proxy-listen=localhost:1055 \
        >"$HOME/.tailscale/tailscaled.log" 2>&1 &

    for _ in $(seq 1 30); do
        [ -S "$HOME/.tailscale/tailscaled.sock" ] && break
        sleep 1
    done

    tailscale --socket="$HOME/.tailscale/tailscaled.sock" up \
        --authkey="$TS_AUTHKEY" \
        --hostname="${TS_HOSTNAME:-starbase-frontend}" \
        --accept-dns=false \
        --timeout=30s

    # Scoped to backend_client.py only -- NOT HTTP_PROXY/HTTPS_PROXY globally,
    # since the tailnet netstack only routes tailnet addresses (no exit node
    # configured), so a global proxy would break NCBI/IPStack/SMTP calls to
    # the public internet.
    export TAILSCALE_PROXY="http://localhost:1055"
    echo "Tailscale up. Backend calls will route via tailnet proxy on :1055."
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

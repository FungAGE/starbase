#!/bin/sh

# Check if the Tailscale auth key is provided
if [ -z "$TS_AUTH_KEY" ]; then
  echo "Error: Tailscale auth key (TS_AUTH_KEY) not set."
  exit 1
fi

if [ -n "$TS_AUTH_KEY" ]; then
  sudo tailscaled --tun=userspace-networking &
  sudo tailscale up --authkey="${TS_AUTH_KEY}" --hostname="starbase-app" --accept-routes &

else
  echo "Error: No Tailscale auth key provided. Please set TS_AUTH_KEY."
  exit 1
fi
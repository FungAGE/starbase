#!/bin/bash

# Simple wrapper script to run the development setup
# This makes it easy to start development without remembering the full path

echo "🚀 Starting STARBASE Local Development Setup..."
echo ""

# Check if the dev directory exists
if [ ! -d "dev" ]; then
    echo "❌ Error: 'dev' directory not found!"
    echo "   Please ensure you're in the correct repository root."
    exit 1
fi

bash ./dev/scripts/setup-local-pod.sh

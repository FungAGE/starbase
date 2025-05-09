#!/bin/bash
set -e

# Start minikube if not running
if ! minikube status &>/dev/null; then
  echo "Starting Minikube..."
  ./kubernetes/minikube-start.sh
else
  echo "Minikube is already running"
fi

# Create persistence directories
mkdir -p ./src/database/db
mkdir -p /tmp/starbase-data

# Sync data between docker compose volumes and minikube volumes
echo "Syncing database files to minikube persistent storage..."
if [ -d "./src/database/db" ] && [ "$(ls -A ./src/database/db)" ]; then
  sudo cp -r ./src/database/db/* /tmp/starbase-data/ || true
  echo "Database files synced to minikube storage"
fi

# Set environment to use minikube's docker daemon
eval $(minikube docker-env)

# Build with docker compose using the minikube version
echo "Building with docker compose using minikube's Docker daemon..."
docker compose -f kubernetes/minikube-compose.yaml build

# Run the app using docker compose with minikube
echo "Running the app with docker compose in minikube environment..."
docker compose -f kubernetes/minikube-compose.yaml up -d app

# Provide connection information
echo "Your app is now running in the minikube environment"
echo "Access it at: http://localhost:8000"
echo ""
echo "To stop the app, run:"
echo "docker compose -f kubernetes/minikube-compose.yaml down"
echo ""
echo "To deploy to Kubernetes instead of using docker compose:"
echo "./kubernetes/build-and-deploy.sh" 
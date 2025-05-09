#!/bin/bash
set -e

# Create directory for database persistence
mkdir -p ./src/database/db
mkdir -p /tmp/starbase-data

# Copy database files to minikube's persistent volume location
echo "Syncing database files..."
sudo cp -r ./src/database/db/* /tmp/starbase-data/ || true

# Ensure we're using Minikube's Docker daemon
eval $(minikube docker-env)

# Build the Docker image
echo "Building Docker image..."
docker build -t starbase:latest .

# Apply Kubernetes manifests
echo "Deploying to Minikube..."
kubectl apply -f kubernetes/persistentvolume.yaml
kubectl apply -f kubernetes/deployment.yaml
kubectl apply -f kubernetes/service.yaml
kubectl apply -f kubernetes/ingress.yaml

# Add hosts entry for local development
echo "Add the following line to your /etc/hosts file if not already present:"
echo "$(minikube ip) starbase.local"

# Wait for deployment to be available
echo "Waiting for deployment to be ready..."
kubectl rollout status deployment/starbase

echo "Deployment complete! Access your app at http://starbase.local"
echo "For debugging, you can also run: kubectl port-forward svc/starbase 8000:80" 
#!/bin/bash

# Simple K3s Deployment Script for Serve-like Testing
# This mimics how your app would be deployed on the actual Serve platform

set -e

echo "🚀 Deploying Starbase to K3s (Serve-like environment)"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Build the Docker image and load it into k3s
echo -e "${BLUE}Building and loading Docker image...${NC}"
docker build -t starbase:local .
k3s ctr images import <(docker save starbase:local)

# Deploy to k3s
echo -e "${BLUE}Deploying to k3s...${NC}"
kubectl apply -f k3s-serve-test/simple-deploy.yaml

# Wait for deployment
echo -e "${BLUE}Waiting for deployment to be ready...${NC}"
kubectl wait --for=condition=available --timeout=300s deployment/starbase-app -n starbase-test

# Add to /etc/hosts if not already there
if ! grep -q "starbase.local" /etc/hosts; then
    echo -e "${YELLOW}Adding starbase.local to /etc/hosts (requires sudo)...${NC}"
    echo "127.0.0.1 starbase.local" | sudo tee -a /etc/hosts
fi

# Get the status
echo -e "${GREEN}✅ Deployment complete!${NC}"
echo
echo "📊 Status:"
kubectl get pods -n starbase-test
echo
echo "🌐 Access your app:"
echo "  - Via Ingress: http://starbase.local"
echo "  - Via Port Forward: kubectl port-forward -n starbase-test svc/starbase-service 8080:80"
echo
echo "🔍 Useful commands:"
echo "  - View logs: kubectl logs -f -n starbase-test deployment/starbase-app"
echo "  - Get pod details: kubectl describe pod -n starbase-test -l app=starbase"
echo "  - Delete deployment: kubectl delete namespace starbase-test" 
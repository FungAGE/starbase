#!/bin/bash
set -e

# Start Minikube with reasonable resources
minikube start --memory=4096 --cpus=2 --driver=docker --force

# Enable ingress addon for routing
minikube addons enable ingress

# Set docker env to use Minikube's Docker daemon
eval $(minikube docker-env)

# Display cluster info
kubectl cluster-info

echo "Minikube Kubernetes cluster is now running!"
echo "Use 'minikube dashboard' to access the Kubernetes dashboard" 
# Starbase Kubernetes Development Environment

This directory contains all necessary files for testing the Starbase application in a local Kubernetes environment before deploying to production.

## Prerequisites

- Docker
- Docker Compose
- Minikube (https://minikube.sigs.k8s.io/docs/start/)
- kubectl (https://kubernetes.io/docs/tasks/tools/)

install minikube:
```bash
curl -LO https://github.com/kubernetes/minikube/releases/latest/download/minikube-linux-amd64
sudo install minikube-linux-amd64 /usr/local/bin/minikube && rm minikube-linux-amd64
```

install kubectl:
```bash
curl -LO "https://dl.k8s.io/release/$(curl -L -s https://dl.k8s.io/release/stable.txt)/bin/linux/amd64/kubectl"
echo "$(cat kubectl.sha256)  kubectl" | sha256sum --check
sudo install -o root -g root -m 0755 kubectl /usr/local/bin/kubectl
kubectl version --client
```
## Getting Started

### Option 1: Use docker compose with Minikube (Easiest)

This option lets you use your familiar docker compose workflow but with Minikube's Docker daemon:

```bash
chmod +x kubernetes/compose-to-kube.sh
./kubernetes/compose-to-kube.sh
```

This will:
1. Start Minikube if not already running
2. Sync your database files to the Minikube environment
3. Build your Docker images using Minikube's Docker daemon
4. Start your containers with docker compose using the Minikube daemon

To stop the services:
```bash
docker compose -f kubernetes/minikube-compose.yaml down
```

### Option 2: Full Kubernetes Deployment

If you want to test with actual Kubernetes resources (deployments, services, etc.):

1. Start the Minikube cluster:

```bash
chmod +x kubernetes/minikube-start.sh
./kubernetes/minikube-start.sh
```

2. Build the Docker image and deploy to Minikube:

```bash
chmod +x kubernetes/build-and-deploy.sh
./kubernetes/build-and-deploy.sh
```

3. Add the host entry shown by the script to your `/etc/hosts` file:

```
# Example (your IP may differ)
192.168.49.2 starbase.local
```

4. Access the application at http://starbase.local

## Managing Database Persistence

Both options use the same persistent storage located at `/tmp/starbase-data` on your host machine, which is mounted into both docker compose and Kubernetes deployments.

## Useful Commands

### View application logs
```bash
# For docker compose
docker compose -f kubernetes/minikube-compose.yaml logs -f app

# For Kubernetes
kubectl logs -l app=starbase
```

### Get a shell inside the container
```bash
# For docker compose
docker compose -f kubernetes/minikube-compose.yaml exec app /bin/bash

# For Kubernetes
kubectl exec -it $(kubectl get pod -l app=starbase -o jsonpath='{.items[0].metadata.name}') -- /bin/bash
```

### Delete the deployment
```bash
kubectl delete -f kubernetes/deployment.yaml -f kubernetes/service.yaml -f kubernetes/ingress.yaml -f kubernetes/persistentvolume.yaml
```

### Stop Minikube
```bash
minikube stop
```

## Testing Updates

After making changes to your application:

1. Rebuild using docker compose:
```bash
docker compose -f kubernetes/minikube-compose.yaml build
docker compose -f kubernetes/minikube-compose.yaml up -d
```

2. Or for Kubernetes:
```bash
./kubernetes/build-and-deploy.sh
``` 
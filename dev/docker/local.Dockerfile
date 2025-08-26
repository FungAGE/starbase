# Local Single-Pod Dockerfile for STARBASE
# This Dockerfile simulates a Kubernetes pod environment for local testing

FROM python:3.9
LABEL org.opencontainers.image.authors="adrian.e.forsythe@gmail.com"
LABEL org.opencontainers.image.description="STARBASE local single-pod deployment for testing"
LABEL org.opencontainers.image.source="https://github.com/genecore/starbase"
LABEL pod.type="single-pod-simulation"
LABEL deployment.type="local-testing"

# Build arguments (passed from docker-compose)
ARG IPSTACK_API_KEY
ARG MAINTENANCE_TOKEN

# Environment variables (mimicking Kubernetes pod environment)
ENV USER=starbase
ENV HOME=/home/$USER
ENV IPSTACK_API_KEY=$IPSTACK_API_KEY
ENV MAINTENANCE_TOKEN=$MAINTENANCE_TOKEN

# Kubernetes-like environment variables
ENV KUBERNETES_SERVICE_HOST=127.0.0.1
ENV KUBERNETES_SERVICE_PORT=443
ENV POD_NAME=starbase-local-pod
ENV POD_NAMESPACE=default
ENV POD_IP=127.0.0.1
ENV NODE_NAME=local-node

# Security context (mimicking Kubernetes security context)
ENV RUN_AS_USER=1000
ENV RUN_AS_GROUP=1000
ENV FS_GROUP=1000

# Add user with specific UID/GID (mimicking Kubernetes security context)
RUN groupadd -g ${FS_GROUP} starbase && \
    useradd -m -u ${RUN_AS_USER} -g ${RUN_AS_GROUP} $USER

# Set working directory
WORKDIR $HOME/

# System dependencies (including tools for Kubernetes-like monitoring)
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
        curl \
        iptables \
        wget \
        redis-server \
        build-essential \
        procps \
        htop \
        net-tools \
        iputils-ping \
        dnsutils \
        jq \
        && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Mambaforge (same as original Dockerfile)
RUN wget https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p $HOME/miniconda && \
    rm miniforge.sh && \
    echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> $HOME/.bashrc && \
    echo "conda activate starbase" >> $HOME/.bashrc

ENV PATH=$HOME/miniconda/bin:$PATH

# Copy environment file and create conda environment using mamba
COPY environment.yaml .
RUN mamba env create -f environment.yaml && \
    mamba clean -afy

# Set conda environment to activate by default
ENV PATH=$HOME/miniconda/envs/starbase/bin:$HOME/miniconda/bin:$PATH
ENV CONDA_DEFAULT_ENV=starbase

# Install Node.js, npm, and blasterjs
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g biojs-vis-blasterjs

# Create directory structure (mimicking Kubernetes volume mounts)
RUN mkdir -p \
    $HOME/src/database/db \
    $HOME/src/database/cache \
    $HOME/src/database/logs && \
    chown -R $USER:$USER $HOME && \
    chmod -R 755 $HOME/src/database/db && \
    chmod -R 755 $HOME/src/database/cache && \
    chmod -R 755 $HOME/src/database/logs

# Add enhanced healthcheck (mimicking Kubernetes liveness/readiness probes)
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/api/cache/status || exit 1

# Copy application code (changes most frequently, so do this last)
COPY ./ ./
RUN chmod +x dev/scripts/start-script.sh && \
    # Ensure all directories and files are owned by starbase user
    chown -R $USER:$USER $HOME/src && \
    chown -R $USER:$USER $HOME/src/database/cache && \
    chown -R $USER:$USER $HOME/src/database/logs

# Create a startup script that mimics Kubernetes init container behavior
RUN cat > $HOME/k8s-init.sh << 'EOF'
#!/bin/bash
echo "=== Kubernetes Pod Simulation Init ==="
echo "Pod Name: $POD_NAME"
echo "Pod Namespace: $POD_NAMESPACE"
echo "Pod IP: $POD_IP"
echo "Node Name: $NODE_NAME"
echo "Kubernetes Service Host: $KUBERNETES_SERVICE_HOST"
echo "Kubernetes Service Port: $KUBERNETES_SERVICE_PORT"
echo "======================================"

# Create directories if they don't exist (mimicking init container behavior)
mkdir -p /home/starbase/src/database/db
mkdir -p /home/starbase/src/database/cache
mkdir -p /home/starbase/src/database/logs

# Set proper permissions (mimicking Kubernetes volume mount behavior)
chown -R starbase:starbase /home/starbase/src/database/db
chown -R starbase:starbase /home/starbase/src/database/cache
chown -R starbase:starbase /home/starbase/src/database/logs

echo "Pod initialization completed successfully"
EOF

RUN chmod +x $HOME/k8s-init.sh && chown $USER:$USER $HOME/k8s-init.sh

# Create enhanced start script for single-pod deployment
RUN cat > $HOME/start-local-pod.sh << 'EOF'
#!/bin/bash

# Run Kubernetes-like initialization
/home/starbase/k8s-init.sh

# Set environment variables for single-pod deployment
export PYTHONPATH=/home/starbase:/home/starbase/src
export ENVIRONMENT=${ENVIRONMENT:-development}
export DEV_MODE=${DEV_MODE:-true}

# Check if we're in development mode
if [ "$DEV_MODE" = "true" ]; then
    echo "Starting STARBASE in single-pod development mode..."
    # Use single worker for local testing (mimicking single pod)
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --workers=1 \
        --reload \
        --log-level=info \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --access-log \
        app:server
else
    echo "Starting STARBASE in single-pod production mode..."
    # Single worker for pod simulation
    uvicorn --host=0.0.0.0 \
        --port=8000 \
        --workers=1 \
        --log-level=warning \
        --interface wsgi \
        --proxy-headers \
        --forwarded-allow-ips='*' \
        --timeout-keep-alive=5 \
        --timeout-graceful-shutdown=60 \
        --limit-max-requests=1000 \
        --no-access-log \
        app:server
fi
EOF

RUN chmod +x $HOME/start-local-pod.sh && chown $USER:$USER $HOME/start-local-pod.sh

# Switch to user (mimicking Kubernetes security context)
USER $USER

# Expose port (same as original)
EXPOSE 8000

# Set labels for Kubernetes-like identification
LABEL pod.name="starbase-local-pod"
LABEL pod.namespace="default"
LABEL app.name="starbase"
LABEL app.version="local"
LABEL deployment.type="single-pod"

# Use the enhanced start script
ENTRYPOINT ["./start-local-pod.sh"]


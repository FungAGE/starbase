# STARBASE Local Single-Pod Testing Environment

This document describes how to set up and use a local testing environment that simulates a single Kubernetes pod deployment for STARBASE.

## Overview

The local single-pod environment mimics the deployment service's Kubernetes pod behavior, allowing you to test your application in conditions similar to production without requiring a full Kubernetes cluster.

### What This Environment Provides

- **Single Pod Simulation**: Runs STARBASE in a single container that mimics a Kubernetes pod
- **Resource Constraints**: Applies memory and CPU limits similar to Kubernetes resource limits
- **Environment Variables**: Sets Kubernetes-like environment variables for testing
- **Volume Mounts**: Simulates persistent volume claims with local directories
- **Health Checks**: Implements liveness and readiness probe equivalents
- **Security Context**: Applies user/group permissions similar to Kubernetes security contexts
- **Network Simulation**: Creates a custom network to simulate pod networking

## Prerequisites

- Docker Engine (version 20.0 or higher)
- Docker Compose (version 2.0 or higher)
- At least 4GB of available RAM
- At least 2 CPU cores

## Quick Start

### 1. Automatic Setup (Recommended)

Use the provided setup script for automated configuration:

```bash
# Make setup script executable (if not already)
chmod +x setup-local-pod.sh

# Run the complete setup
./setup-local-pod.sh
```

This script will:
- Detect your local IP address
- Create the environment configuration file
- Build the Docker image
- Create required directories
- Start the single-pod environment
- Verify the application is healthy

### 2. Manual Setup

If you prefer manual setup or need to customize the configuration:

#### Step 1: Create Environment Configuration

```bash
# Copy the environment template
cp env.template .env

# Edit .env file and update the following variables:
# - AUTH_DOMAIN: Set to your local IP (e.g., 192.168.1.10.nip.io)
# - HOST_IP: Set to your local IP address
# - IPSTACK_API_KEY: (optional) Your IP geolocation API key
# - MAINTENANCE_TOKEN: (optional) Your maintenance token
```

#### Step 2: Create Required Directories

```bash
mkdir -p data cache logs src/database/db
```

#### Step 3: Build the Docker Image

```bash
docker build -f local.Dockerfile -t starbase:local-pod .
```

#### Step 4: Start the Environment

```bash
docker-compose -f docker-compose.local.yaml up -d
```

## Configuration

### Environment Variables

The `.env` file contains all configuration options. Key variables include:

| Variable | Description | Default |
|----------|-------------|---------|
| `AUTH_DOMAIN` | Domain for the application | `127.0.0.1.nip.io` |
| `HOST_IP` | Local IP address | `127.0.0.1` |
| `EXTERNAL_PORT` | External port for the application | `8080` |
| `ENVIRONMENT` | Application environment | `development` |
| `DEV_MODE` | Enable development features | `true` |
| `MEMORY_LIMIT` | Container memory limit | `2G` |
| `CPU_LIMIT` | Container CPU limit | `1.0` |
| `POD_NAME` | Simulated pod name | `starbase-local-pod` |
| `POD_NAMESPACE` | Simulated pod namespace | `default` |

### Resource Limits

The environment simulates Kubernetes resource constraints:

- **Memory Limit**: 2GB (configurable via `MEMORY_LIMIT`)
- **CPU Limit**: 1 CPU core (configurable via `CPU_LIMIT`)
- **Memory Request**: 1GB (configurable via `MEMORY_REQUEST`)
- **CPU Request**: 0.5 CPU cores (configurable via `CPU_REQUEST`)

### Volume Mounts

The following directories are mounted to simulate Kubernetes persistent volumes:

- `./data` → `/home/starbase/data` (application data)
- `./cache` → `/home/starbase/src/database/db/cache` (cache storage)
- `./logs` → `/home/starbase/logs` (application logs)
- `./src` → `/home/starbase/src` (source code for development)

## Usage

### Accessing the Application

Once the environment is running, you can access STARBASE at:

- **Local**: http://localhost:8080
- **Network**: http://YOUR_LOCAL_IP:8080
- **Domain**: http://studio.YOUR_LOCAL_IP.nip.io:8080

### Common Commands

#### View Container Status
```bash
docker ps | grep starbase
```

#### View Application Logs
```bash
docker logs starbase-local-pod -f
```

#### Access Container Shell
```bash
docker exec -it starbase-local-pod /bin/bash
```

#### Restart the Application
```bash
docker-compose -f docker-compose.local.yaml restart
```

#### Stop the Environment
```bash
docker-compose -f docker-compose.local.yaml down
```

#### Rebuild and Restart
```bash
docker-compose -f docker-compose.local.yaml down
docker build -f local.Dockerfile -t starbase:local-pod .
docker-compose -f docker-compose.local.yaml up -d
```

### Health Monitoring

The environment includes health checks that mimic Kubernetes liveness and readiness probes:

- **Health Check Endpoint**: http://localhost:8080/api/cache/status
- **Check Interval**: 30 seconds
- **Timeout**: 10 seconds
- **Retries**: 3

#### Check Application Health
```bash
curl -f http://localhost:8080/api/cache/status
```

#### View Health Status
```bash
docker inspect starbase-local-pod | jq '.[0].State.Health'
```

## Troubleshooting

### Common Issues

#### Port Already in Use
If port 8080 is already in use, change the `EXTERNAL_PORT` in your `.env` file:
```bash
EXTERNAL_PORT=8081
```

#### Permission Denied Errors
Ensure the user has proper permissions for Docker:
```bash
sudo usermod -aG docker $USER
newgrp docker
```

#### Container Health Check Failures
Check the application logs for startup errors:
```bash
docker logs starbase-local-pod
```

#### DNS Resolution Issues
If using `.nip.io` domains doesn't work, try using the IP address directly or add entries to `/etc/hosts`.

### Debugging

#### Enable Debug Logging
Set the following in your `.env` file:
```bash
LOG_LEVEL=DEBUG
DEV_MODE=true
```

#### Access Container Filesystem
```bash
docker exec -it starbase-local-pod /bin/bash
```

#### View Container Resource Usage
```bash
docker stats starbase-local-pod
```

#### Inspect Container Configuration
```bash
docker inspect starbase-local-pod
```

## Development Workflow

### Making Code Changes

When developing with this environment:

1. **Source Code**: Changes to `./src/` are automatically mounted and visible in the container
2. **Configuration**: Changes to configuration files require a container restart
3. **Dependencies**: Changes to `environment.yaml` or `pyproject.toml` require rebuilding the image

### Live Reloading

The development mode enables live reloading. Changes to Python files will automatically restart the application.

### Testing

Run tests inside the container:
```bash
docker exec -it starbase-local-pod pytest
```

Or run the test profile:
```bash
docker-compose -f docker-compose.local.yaml --profile testing up
```

## Comparison with Kubernetes

This local environment simulates many aspects of a Kubernetes pod deployment:

### Similarities
- Single container/pod architecture
- Resource limits and requests
- Environment variable injection
- Volume mounts
- Health checks (liveness/readiness probes)
- Security context (user/group permissions)
- Network isolation
- Labels and annotations

### Differences
- No actual Kubernetes API server
- No service discovery (using static IPs/domains)
- No automatic pod scheduling/rescheduling
- No horizontal pod autoscaling
- No ingress controller (using direct port mapping)
- No secrets management (using environment variables)

## Advanced Configuration

### Custom Resource Limits

Modify resource limits in `docker-compose.local.yaml`:
```yaml
deploy:
  resources:
    limits:
      memory: 4G
      cpus: '2.0'
    reservations:
      memory: 2G
      cpus: '1.0'
```

### Adding Additional Services

To add services like Redis or PostgreSQL, uncomment and modify the relevant sections in `docker-compose.local.yaml`.

### Custom Network Configuration

Modify the network configuration to simulate different network topologies:
```yaml
networks:
  default:
    driver: bridge
    ipam:
      config:
        - subnet: 172.30.0.0/16
```

## Monitoring and Observability

### Container Metrics
```bash
# View real-time resource usage
docker stats starbase-local-pod

# View container processes
docker exec starbase-local-pod ps aux
```

### Application Metrics
Access the metrics endpoint (if available):
```bash
curl http://localhost:8080/metrics
```

### Log Aggregation
Configure log forwarding in `docker-compose.local.yaml`:
```yaml
logging:
  driver: "fluentd"
  options:
    fluentd-address: localhost:24224
    tag: starbase.local
```

## Security Considerations

The local environment implements security practices similar to Kubernetes:

- **Non-root user**: Application runs as user `starbase` (UID 1000)
- **Read-only filesystem**: Some directories are mounted read-only
- **No new privileges**: Security option prevents privilege escalation
- **Resource limits**: Prevents resource exhaustion

## Migration to Production

When moving from local testing to production Kubernetes:

1. **Environment Variables**: Use Kubernetes ConfigMaps and Secrets
2. **Volumes**: Replace local volumes with Persistent Volume Claims
3. **Networking**: Configure Kubernetes Services and Ingress
4. **Health Checks**: Convert Docker health checks to Kubernetes probes
5. **Resource Limits**: Apply ResourceQuotas and LimitRanges
6. **Security**: Implement Pod Security Standards

## Support

For issues specific to the local single-pod environment:

1. Check the troubleshooting section above
2. Review container logs: `docker logs starbase-local-pod`
3. Verify environment configuration in `.env` file
4. Ensure all prerequisites are met

For application-specific issues, refer to the main STARBASE documentation.


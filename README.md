# `starbase`: Database and Toolkit for Exploration of _Starship_ Elements in Fungi

<img src=assets/logos/favicon.png width=200>

<!-- badges: start -->

![Starbase status](https://img.shields.io/website?url=https%3A%2F%2Fstarbase.serve.scilifelab.se%2Fapp%2Fstarbase)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/FungAGE/starbase/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FungAGE/starbase?branch=main)

<!-- badges: end -->

## Overview

`starbase` is a web-based application that provides various analytical and visualization tools for exploring large eukaryotic transposable elements known as _Starships_.

Access `starbase` [here](https://starbase.serve.scilifelab.se/).

## What is a _Starship_?

_Starships_ are extremely large (~20-700 kb) DNA transposons that can constitute up to 2% of fungal genomes. These novel class II DNA transposons are endemic to _Pezizomycotina_ and can significantly impact the architecture of fungal genomes. They replicate within the host genome via tyrosine recombinases (captain genes) and can carry significant genetic 'cargo', such as:

- Genes for metal resistance in _Paecilomyces_
- Cheese-making genes in _Penicillium_
- Formaldehyde resistance genes in _Aspergillus nidulans_ and _Penicillium chrysogenum_

## Usage

- [Explore existing _Starship_ sequences](https://starbase.serve.scilifelab.se/wiki): Browse and select individual or collections of _Starship_ sequences for download. View detailed information about _Starship_ Families and their general characteristics.
- [Search for _Starships_](https://starbase.serve.scilifelab.se/blast): Utilize the BLAST/HMMER search functions to analyze sequences against the existing database. Runs a classification workflow to compare identify if the sequence is within an existing family/navis/haplotype.
- [Submit a Sequence](https://starbase.serve.scilifelab.se/submit): Use the submission tab to upload new _Starship_ sequences for curation and inclusion in the database.

### Features Under Development

- `starfish` webserver
- [_Starship_ `Browser/Comparison](https://starbase.serve.scilifelab.se/pgv): Visualize and compare _Starships_ and their gene annotations.

## Contributions

Contributions to the development of `starbase` are welcome! If you have suggestions or improvements, please feel free to submit an issue or pull request.

### Building and running from the Docker image

- **Note:** Currently, `starbase` will not run without database files (`src/database/db/`), which need to be mounted into the container.

Versions of the database are stored in a [Zenodo repository](https://doi.org/10.5281/zenodo.17533381).

## Local Development

### Mimicking [Serve](https://github.com/ScilifelabDataCentre/serve?tab=readme-ov-file#deploy-serve-for-local-development-with-docker-compose) Production Environment (Recommended)
The SciLifeLab Serve platform uses a single-pod kubernetes environment for each app.

### Quick Start

```bash
# Option 1: Simple setup (recommended)
./setup-dev.sh

# Option 2: Direct setup
bash ./dev/scripts/setup-local-pod.sh
```

This will:
- Detect your local IP address
- Create the environment configuration
- Build the Docker image
- Create required directories
- Start the single-pod environment
- Verify the application is healthy

#### Manual Setup

If you prefer manual setup or need to customize the configuration:

```bash
# 1. Create environment configuration
cp dev/config/env.template .env

# 2. Update .env with your local IP (optional - will be auto-detected)
# Edit AUTH_DOMAIN, HOST_IP, POD_IP, KUBERNETES_SERVICE_HOST

# 3. Create required directories
mkdir -p src/database/cache src/database/logs src/database/db

# 4. Build the Docker image
docker build -f dev/docker/local.Dockerfile -t starbase:local-pod .

# 5. Start the environment
docker compose -f dev/docker/docker-compose.local.yaml up -d
```

#### Access Your Application

Once running, access your Dash app at:
- **Local**: http://localhost:8080
- **Network**: http://YOUR_LOCAL_IP:8080
- **Domain**: http://studio.YOUR_LOCAL_IP.nip.io:8080

#### Development Workflow

```bash
# View application logs
docker logs starbase-local-pod -f

# Access container shell
docker exec -it starbase-local-pod /bin/bash

# Restart the application
docker compose -f dev/docker/docker-compose.local.yaml restart

# Stop the environment
docker compose -f dev/docker/docker-compose.local.yaml down

# Rebuild and restart (after dependency changes)
docker compose -f dev/docker/docker-compose.local.yaml down
docker build -f dev/docker/local.Dockerfile -t starbase:local-pod .
docker compose -f dev/docker/docker-compose.local.yaml up -d
```

#### Features

This environment provides:
- **Single Pod Simulation**: Mimics a Kubernetes pod deployment
- **Resource Constraints**: 2GB memory, 1 CPU core limits
- **Volume Mounts**: Database files properly mounted and persistent
- **Health Checks**: Application health monitoring
- **Development Mode**: Hot reloading enabled for code changes
- **Security Context**: Runs as non-root user (UID 1000)
- **Network Isolation**: Custom bridge network

### Legacy Development Methods

#### Docker Compose (Original)

```bash
# Build and run the application
docker compose up app --build

# Run with verbose logging and hot reload
docker compose up app --build && docker exec starbase_app ./start-script.sh --dev
```

#### Direct Docker

```bash
# Find the most recent version under the "Packages" tab
docker pull ghcr.io/fungage/starbase:[tag]

# Run Starbase locally with Docker
docker build -t starbase .
docker run -it --rm -p 8000:8000 starbase ./start-script.sh
```

You can reach the app in your local browser by visiting `localhost:8000` (legacy) or `localhost:8080` (single-pod environment).

### Debug Mode

Enable debug logging in your `.env` file:
```bash
LOG_LEVEL=DEBUG
DEV_MODE=true
```

## Unit tests

Run Pytest in a Docker container:

```bash
# Run tests (will build image if needed)
docker compose --profile testing up unit-tests

# Or run tests against existing image
docker compose run --rm unit-tests
```

## Tear down

```bash
# Stop all services and remove volumes
docker compose down --volumes

# Remove built images (optional)
docker rmi starbase:latest
```
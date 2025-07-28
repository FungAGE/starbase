# `starbase`: Database and Toolkit for Exploration of _Starship_ Elements in Fungi

<img src=assets/logos/favicon.png width=200>

<!-- badges: start -->

![Starbase status](https://img.shields.io/website?url=https%3A%2F%2Fstarbase.serve.scilifelab.se%2Fapp%2Fstarbase)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/FungAGE/starbase/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FungAGE/starbase?branch=main)

<!-- badges: end -->

# `starbase`: A Database and Toolkit for Exploring Eukaryotic Transposable Elements in Fungi

## Overview

`starbase` is a web-based application that provides various analytical and visualization tools for exploring large eukaryotic transposable elements known as _Starships_.

Access `starbase` [here](https://starbase.serve.scilifelab.se/).

## What is a _Starship_?

_Starships_ are extremely large (~20-700 kb) DNA transposons that can constitute up to 2% of fungal genomes. These novel class II DNA transposons are endemic to _Pezizomycotina_ and can significantly impact the architecture of fungal genomes. They replicate within the host genome via tyrosine recombinases (captain genes) and can carry significant genetic 'cargo', such as:

- Genes for metal resistance in _Paecilomyces_
- Cheese-making genes in _Penicillium_
- Formaldehyde resistance genes in _Aspergillus nidulans_ and _Penicillium chrysogenum_

## Usage

- [Download _Starship_ sequences](https://starbase.serve.scilifelab.se/download): Browse and select individual or collections of _Starship_ sequences for download.
- [Access the Wiki](https://starbase.serve.scilifelab.se/wiki): View detailed information about _Starship_ Families and their general characteristics.
- [Search for _Starships_](https://starbase.serve.scilifelab.se/blast): Utilize the BLAST/HMMER search functions to analyze sequences against the existing database.
- [Submit a Sequence](https://starbase.serve.scilifelab.se/submit): Use the submission tab to upload new _Starship_ sequences for curation and inclusion in the database.

### Features Under Development

- `starfish` webserver
- [_Starship_ Browser/Comparison](https://starbase.serve.scilifelab.se/pgv): Visualize and compare _Starships_ and their gene annotations.

## Contributions

Contributions to the development of `starbase` are welcome! If you have suggestions or improvements, please feel free to submit an issue or pull request.

### Building and running from the Docker image

- **Note:** Currently, `starbase` will not run without database files (`src/database/db/`), which need to be mounted into the container. We will provide a link to an archive of the database files, when available.

## Local development

* Find the most recent version under the "Packages" tab: `docker pull ghcr.io/fungage/starbase:[tag]`

* Run Starbase locally with Docker:
```bash
docker build -t starbase .
docker run -it --rm -p 8000:8000 starbase ./start-script.sh
```

* **Recommended**: Run Starbase locally using Docker Compose:
```bash
# Build and run the application
docker compose up app --build

# Run with verbose logging and hot reload
docker compose up app --build && docker exec starbase_app ./start-script.sh --dev
```

You can reach the app in your local browser by visiting `localhost:8000`.

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
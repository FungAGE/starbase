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

1. Find the most recent version under the "Packages" tab in this repository: `docker pull ghcr.io/fungage/starbase:[tag]`
2. Run `docker build -t starbase .` then `docker run --rm -p 8000:8000 starbase`
3. Launch the app by visiting `localhost:8000` in your browser.

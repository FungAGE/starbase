# `starbase`: Database and Toolkit for Exploration of _Starship_ Elements in Fungi

<img src=assets/logos/favicon.png width=200>

<!-- badges: start -->

![Starbase status](https://img.shields.io/website?url=https%3A%2F%2Fstarbase.serve.scilifelab.se%2Fapp%2Fstarbase)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/FungAGE/starbase/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FungAGE/starbase?branch=main)

<!-- badges: end -->

# `starbase`: A Database and Toolkit for Exploring Eukaryotic Transposable Elements in Fungi

## Overview

`starbase` is a web-based application that provides various analytical and visualization tools for exploring large eukaryotic transposable elements known as Starships.

Access `starbase` [here](https://starbase.serve.scilifelab.se/).

## What is a Starship?

Starships are extremely large (~20-700 kb) DNA transposons that can constitute up to 2% of fungal genomes. These novel class II DNA transposons are endemic to _Pezizomycotina_ and can significantly impact the architecture of fungal genomes. They replicate within the host genome via tyrosine recombinases (captain genes) and can carry significant genetic 'cargo', such as:

- Genes for metal resistance in _Paecilomyces_
- Cheese-making genes in _Penicillium_
- Formaldehyde resistance genes in _Aspergillus nidulans_ and _Penicillium chrysogenum_

## Usage

- Access the Wiki: View detailed information about Starship Families and their general characteristics.
- Submit a Sequence: Use the submission tab to upload new Starship sequences for curation and inclusion in the database.
- Search for Starships: Utilize the BLAST/HMMER search functions to analyze sequences against the existing database.

### Features Under Development

- `starfish` webserver
- Synteny/Genome Browser

## Contributions

Contributions to the development of `starbase` are welcome! If you have suggestions or improvements, please feel free to submit an issue or pull request.

## Build from the docker image

1. Find the most recent version under the "Packages" tab: `docker pull ghcr.io/fungage/starbase:[tag]`
1. Run `docker build -t starbase` then `docker run -p 7000:80 starbase`
1. Launch the app by visiting `localhost:7000` in your browser.

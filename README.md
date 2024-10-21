# `starbase`: Database and Toolkit for Exploration of _Starship_ Elements in Fungi

<img src=assets/logos/favicon.png width=200>

<!-- badges: start -->

![Starbase status](https://img.shields.io/website?url=https%3A%2F%2Fstarbase.serve.scilifelab.se%2Fapp%2Fstarbase)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test coverage](https://codecov.io/gh/FungAGE/starbase/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FungAGE/starbase?branch=main)

<!-- badges: end -->

# `starbase`: A Database and Toolkit for Exploring Eukaryotic Transposable Elements in Fungi
## Overview

Starbase is a web-based application designed for exploring large eukaryotic transposable elements known as Starships. These novel class II DNA transposons are endemic to *Pezizomycotina* and can significantly impact the genomic architecture of fungi. This application provides various tools for analyzing and visualizing these genetic elements.

## Usage
- Access the Wiki: Click on the "Catalogue/Wiki of Starship Metadata" button to view detailed information about the metadata associated with various Starship sequences.
- Submit a Sequence: Use the submission button to upload new Starship sequences for analysis.
- Perform Searches: Utilize the BLAST/HMMER search features to analyze submitted sequences against existing databases.

## Contribution
Contributions to this project are welcome! If you have suggestions or improvements, please feel free to submit an issue or pull request.

## Features
### Current Functionalities

The application currently supports the following features:
- Wiki: Catalogue/Wiki of Starship Metadata
- Submit: Submission of new Starship sequences
- BLAST/HMMER Searches: Perform BLAST or HMMER searches on the submitted sequences

### Future Development
The following functionalities are under active development and will be available soon:
- Synteny/Genome Browser

## What is a Starship?

Starships are extremely large (~20-700 kb) DNA transposons that can constitute up to 2% of fungal genomes. They replicate within the host genome via tyrosine recombinases (captain genes) and can carry significant genetic 'cargo', such as:

- Genes for metal resistance in *Paecilomyces*
- Cheese-making genes in *Penicillium*
- Formaldehyde resistance genes in *Aspergillus nidulans* and *Penicillium chrysogenum*

## Getting Started
### Current Web-App
Access `starbase` [here](https://starbase.serve.scilifelab.se/).

### Build from the docker image
1. Find the most recent version under the "Packages" tab: `docker pull ghcr.io/fungage/starbase:[tag]`
1. Run `docker build -t starbase` then `docker run -p 7000:80 starbase`
2. Launch the app by visiting `localhost:7000` in your browser.

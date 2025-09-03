# Starbase Database Cleanup & Validation Pipeline

This pipeline analyzes, validates, and repairs data integrity across Starbase tables, focusing on sequence identity, referential integrity, and taxonomy correctness. A single runner orchestrates tasks, and flags allow targeted operations.

## Overview

Major goals:
- Validate/repair relationships among `ships`, `accessions`, and `joined_ships`
- Detect/consolidate duplicate sequences (including reverse complements)
- Validate/fix `joined_ships.ship_id` using the `sequence_reference` table
- Track discovered issues in a dedicated `cleanup_issues` table
- Fill or reconcile missing taxonomy for `genomes` using multiple strategies

Primary entrypoint: `src/database/cleanup/run_comprehensive_cleanup.py`

Core utilities: `src/database/cleanup/utils/database_cleanup.py`

## Quick Start

Analyze only (no writes):
```bash
python src/database/cleanup/run_comprehensive_cleanup.py --analyze-relationships
python src/database/cleanup/run_comprehensive_cleanup.py --check-relationships
```

Full pipeline (dry run):
```bash
python src/database/cleanup/run_comprehensive_cleanup.py --report /tmp/cleanup_report.txt --skip-accession-cleanup
```

Apply general fixes after analysis:
```bash
python src/database/cleanup/run_comprehensive_cleanup.py --apply-fixes --apply
```

## What the Pipeline Does

### 1) Accession cleanup (optional)
- Skipped with `--skip-accession-cleanup`.
- Reverse-complement checks and nested sequence audits to reduce accession duplication.

### 2) Integrity checks
- Genome table: missing fields, orphans
- Taxonomy table: missing levels, orphans, duplicates
- Genomic features: coordinate sanity, orphans, invalid strands
- Foreign keys: broken ship/accession/joined/captain/genome-taxonomy links
- Schema violations: null/empty required fields

### 3) Ships ↔ Accessions ↔ JoinedShips consistency
- Checks for orphaned rows, missing sequence links, `ships` missing from `joined_ships`, inconsistent `joined_ships`.
- Fixes (with `--apply-fixes`): remove broken references, delete truly orphaned accessions, create missing `joined_ships` for `ships`, fill minimal defaults.

### 4) Duplicate sequences consolidation
- De-duplicates `ships` by canonical MD5 (min of forward vs reverse complement) and updates dependent references in `joined_ships`.

### 5) Ship ID validation and repair via sequence_reference
- Uses FASTA-derived `sequence_reference` to validate `joined_ships.ship_id` by sequence identity (MD5 + rev-comp awareness).
- Flags mismatches and adds missing `ship_id` where unambiguous.

### 6) Issue tracking table
- `cleanup_issues` is auto-created and populated after analysis.
- Columns: `id, issue_type, category, table_name, record_id, details(JSON), status, source, created_at, updated_at`.
- Unique index prevents duplicates.

## CLI Reference (selected)

General:
- `--report <path>` save text summary
- `--apply` perform writes (omit for dry run)
- `--skip-accession-cleanup` skip slow accession step

Relationship analysis/fixes:
- `--analyze-relationships` detailed stats and recommendations
- `--check-relationships` counts of core issues
- `--apply-fixes` run general fixes (ships/accessions/joined_ships + ship_id fixes)

Duplicate handling:
- `--identify-duplicates`
- `--cleanup-duplicates [--apply]`
- `--consolidate-duplicate-ships [--apply]`

Ship ID validation/repair:
- `--validate-ship-relationships`
- `--fix-ship-relationships [--apply]`

Taxonomy enrichment (genome.taxonomy_id):
- From joined_ships tax_id mode:
  - `--fix-missing-genome-taxonomy-from-joined [--apply]`
- From NCBI E-utilities (needs OME→scientific name CSV):
  - `--fix-missing-genome-taxonomy-via-ncbi --ome-name-map <csv> [--ncbi-email <email>] [--ncbi-api-key <key>] [--apply]`
- From NCBI taxdump (offline):
  - `--fix-missing-genome-taxonomy-via-taxdump --names-dmp <path> [--include-synonyms] [--apply]`
- From curated OME map TSV (preferred if available):
  - `--fix-missing-genome-taxonomy-via-map --ome-tax-map <tsv> [--apply]`
  - Parser prefers JSON block when present; builds `Genus species`; links to existing taxonomy by name or (genus,species), or creates minimal row.
  - Genus-only fallback is conservative (only if unique match and no species in map).
- Reconcile existing assignments to the map:
  - `--reconcile-genome-taxonomy-via-map --reconcile-ome-tax-map <tsv> [--apply]`

Mislabeling detection:
- `--analyze-mislabeling` heuristics for `joined_ships.ship_id` misuse
- `--fix-mislabeling [--apply]`

## Recording Issues to cleanup_issues

Runner automatically:
1) Ensures `cleanup_issues` exists
2) Records all discovered issues after analysis

You’ll see logs like: `Recorded issues summary: inserted=<N> skipped=<M>`

## Sequence Reference Table

Expected columns (built by a separate tool): `starshipID, sequence, md5_original, md5_revcomp, canonical_md5`.
Used as the source of truth to validate `joined_ships.ship_id` via sequence identity.

## Safety & Execution Model

- Dry run by default; add `--apply` to write
- DB writes wrapped in transactions; rollback on error
- Extensive logging for visibility and troubleshooting

## Typical Workflows

Minimal audit + issue logging:
```bash
python src/database/cleanup/run_comprehensive_cleanup.py --report /tmp/cleanup.txt --skip-accession-cleanup
```

Full pass with repairs:
```bash
python src/database/cleanup/run_comprehensive_cleanup.py --report /tmp/cleanup.txt --apply-fixes --apply
```

Taxonomy completion (fastest → richest):
```bash
# 1) From joined_ships
python src/database/cleanup/run_comprehensive_cleanup.py --fix-missing-genome-taxonomy-from-joined --apply

# 2) Curated OME map and reconciliation
python src/database/cleanup/run_comprehensive_cleanup.py --fix-missing-genome-taxonomy-via-map --ome-tax-map /path/to/ome_map.tsv --apply
python src/database/cleanup/run_comprehensive_cleanup.py --reconcile-genome-taxonomy-via-map --reconcile-ome-tax-map /path/to/ome_map.tsv --apply

# 3) Optionally supplement via taxdump or NCBI
python src/database/cleanup/run_comprehensive_cleanup.py --fix-missing-genome-taxonomy-via-taxdump --names-dmp /path/to/names.dmp --apply
# or
python src/database/cleanup/run_comprehensive_cleanup.py --fix-missing-genome-taxonomy-via-ncbi --ome-name-map /path/to/ome_to_name.csv --ncbi-email you@org --apply
```

## Notes & Best Practices

- Always dry-run first to understand impact
- Prefer curated OME map reconciliation over genus-only inference to avoid mislabeling
- After sequence-based ship_id fixes, re-run duplicate sequence checks
- `cleanup_issues` insertions are idempotent (unique index)

## Where Things Live

- Runner: `src/database/cleanup/run_comprehensive_cleanup.py`
- Core logic: `src/database/cleanup/utils/database_cleanup.py`
- Issue table: `create_cleanup_issues_table`, `record_cleanup_issues` in utils
- Additional docs:
  - `README_database_cleaning.md`: cleaning tasks
  - `README_relationship_checks.md`: core relationship checks

## Troubleshooting

- Import/path errors: run from project root; runner adjusts `sys.path`
- SQLite DateTime errors: pipeline uses proper `datetime` objects
- Parameter binding errors: sequence/ship_id logic uses raw `sqlite3` where needed
- Nothing recorded to `cleanup_issues`: ensure report step ran; check logs for table creation

# Original plan for database cleaning tasks
## correction of `accession_tag`
### checking reverse compliment of sequence
- for each sequence:
1. generate md5sums for both normal and reverse compliment
2. generate list of accessions with identical md5sums for normal and reverse complimented sequences
3. check for existing sequences in the database that are nested, and therefore should be grouped under the same accession
6. if the list is longer than 1, consolidate the accessions based on the lowest number
### checking of nesting and updating accession version
- for each sequence:
1. look for existing seqeuences that are nested in other sequences
2. create a list of these nested sequences and the accessions of both seqeuences
3. consolidate nested sequences under the accession of the longer sequence (i.e. the sequence that they are nested within)
- once accessions have been corrected, proceed with other database cleaning tasks
## check genome table
- we need to add genome information where it is missing
## check taxonomic table
- check that taxonomic information is internally consistent
## check genomic features
- check that coordinates make sense and are linked to genomes
## after all previous steps have been completed
- update foreign keys and check consistency
## general schema check?
- can we check if database entries violate schema?

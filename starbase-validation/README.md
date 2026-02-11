# Starbase Validation

Standalone validation/ETL app that sits between the Django curation app (MySQL) and the Starbase web app (SQLite). It extracts data, validates against quality rules, enriches metadata, and loads versioned SQLite databases only when quality gates pass.

## Architecture

```
Django (MySQL) → Validation App → Versioned SQLite
                      ↓ (if issues)
                 Report to curators
```

## Install

```bash
pip install -e /path/to/starbase-validation
```

## Run pipeline

```bash
# Dry run (validate only, no SQLite output)
starbase-validate --mysql-config=config.yaml --dry-run

# Full run
starbase-validate --mysql-config=config.yaml --output=starbase_v1.0.0.db
```

## Layout

- **starbase_validation/config/** – Quality gates and enrichment source config (YAML).
- **starbase_validation/validators/** – Taxonomy, genome, classification, sequence, relationship checks.
- **starbase_validation/extractors/** – Read from Django MySQL.
- **starbase_validation/transformers/** – Enrich (e.g. NCBI) and normalize.
- **starbase_validation/loaders/** – Build SQLite and manage versions.
- **starbase_validation/reporters/** – Validation reports and issue export to Django.
- **starbase_validation/pipeline.py** – ETL orchestration.

## Quality gates

See `config/quality_rules.yaml` for blocking vs warning rules and thresholds. Blocking issues stop publication; warnings are reported but allow publication.

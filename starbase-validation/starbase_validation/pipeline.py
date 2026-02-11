"""Main ETL orchestration: extract → validate → enrich → transform → gate → load → report."""

import argparse
import logging
from pathlib import Path

import yaml

from starbase_validation.config import get_quality_rules_path
from starbase_validation.extractors import extract_from_mysql
from starbase_validation.loaders import DatabaseVersionManager, load_to_sqlite
from starbase_validation.reporters import generate_report, report_to_django
from starbase_validation.transformers import enrich_taxonomy, normalize
from starbase_validation.validators import validate_all

logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Raised when quality gates fail (blocking issues)."""

    pass


class ValidationPipeline:
    """Orchestrates extract, validate, enrich, transform, gate, load, and report."""

    def __init__(self, mysql_config, quality_rules=None, quality_rules_path=None):
        self.mysql_config = mysql_config
        if quality_rules is not None:
            self.quality_rules = quality_rules
        else:
            path = quality_rules_path or get_quality_rules_path()
            with open(path) as f:
                self.quality_rules = yaml.safe_load(f)
        self.version_manager = DatabaseVersionManager()

    def run(self, output_path=None, output_dir=None, dry_run=False):
        """Run the full pipeline. On blocking issues, report to curators and raise ValidationError."""
        logger.info("Extracting data from Django MySQL...")
        records = extract_from_mysql(self.mysql_config)

        logger.info("Validating data quality...")
        validation_results = validate_all(records, self.quality_rules)

        if validation_results.has_blocking_issues():
            report_to_django(
                validation_results.blocking_issues,
                self.mysql_config,
            )
            raise ValidationError(
                f"{len(validation_results.blocking_issues)} blocking issues found"
            )

        logger.info("Enriching metadata from external sources...")
        enriched = enrich_taxonomy(records)

        logger.info("Transforming data...")
        transformed = normalize(enriched)

        logger.info("Running final validation...")
        final_results = validate_all(transformed, self.quality_rules)
        if not final_results.is_acceptable():
            report_to_django(final_results.blocking_issues, self.mysql_config)
            raise ValidationError("Final validation failed")

        if dry_run:
            logger.info("[DRY RUN] Would create database")
            generate_report(validation_results, output_path=output_path)
            return "dry_run"

        version = self.version_manager.calculate_version(transformed)
        if output_path and output_path.endswith(".db"):
            out_path = output_path
        elif output_dir or (output_path and not output_path.endswith(".db")):
            out_path = str(Path(output_dir or output_path) / f"starbase_v{version}.db")
        else:
            out_path = f"starbase_v{version}.db"

        logger.info("Loading into SQLite: %s", out_path)
        load_to_sqlite(transformed, out_path)
        new_version = self.version_manager.calculate_version(transformed)
        logger.info("Created database version %s", new_version)

        generate_report(validation_results, version=new_version, output_path=output_path)
        return new_version


def main():
    """CLI entry point for starbase-validate."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    parser = argparse.ArgumentParser(description="Starbase validation pipeline")
    parser.add_argument(
        "--mysql-config",
        required=True,
        help="Path to YAML file with MySQL connection (or Django DATABASES-style dict)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output SQLite path (or directory for versioned filename)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate only; do not create SQLite",
    )
    args = parser.parse_args()

    with open(args.mysql_config) as f:
        mysql_config = yaml.safe_load(f)

    pipeline = ValidationPipeline(mysql_config=mysql_config)
    try:
        output_path = args.output
        output_dir = None
        if output_path and not output_path.endswith(".db"):
            output_dir = output_path
            output_path = None
        version = pipeline.run(
            output_path=output_path,
            output_dir=output_dir,
            dry_run=args.dry_run,
        )
        if args.dry_run:
            print("Dry run completed successfully.")
        else:
            print(f"Successfully created version {version}")
    except ValidationError as e:
        logger.error("Validation failed: %s", e)
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Genome metadata validation (assembly, taxonomy link)."""


def validate_genome(records, rules):
    """Check genome metadata completeness. Return ValidationResults-compatible list."""
    from starbase_validation.validators.results import ValidationResults

    results = ValidationResults()
    # TODO: Extract from database_cleanup genome validation
    # Rules: max_missing_assembly, require_taxonomy_link
    return results

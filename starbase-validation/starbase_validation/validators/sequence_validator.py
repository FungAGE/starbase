"""Sequence validation (duplicates, reverse complement checks)."""


def validate_sequences(records, rules):
    """Check for duplicate accessions and reverse complement consistency. Return ValidationResults-compatible list."""
    from starbase_validation.validators.results import ValidationResults

    results = ValidationResults()
    # TODO: Extract from cleanup_accessions duplicate detection, reverse complement checks
    # Rules: max_duplicate_accessions, require_reverse_complement_check
    return results

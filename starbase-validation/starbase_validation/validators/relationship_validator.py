"""Foreign key and relationship validation."""


def validate_relationships(records):
    """Check for orphaned foreign keys and referential integrity. Return ValidationResults-compatible list."""
    from starbase_validation.validators.results import ValidationResults

    results = ValidationResults()
    # TODO: Orphaned FK checks
    # blocking_issues: orphaned_foreign_keys
    return results

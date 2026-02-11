"""Taxonomy completeness and hierarchy validation."""


def validate_taxonomy(records, rules):
    """Check taxonomy completeness (species, strain, genus) and hierarchy. Return ValidationResults-compatible list."""
    from starbase_validation.validators.results import ValidationResults

    results = ValidationResults()
    # TODO: Extract from taxonomy_consolidator.check_taxonomy_completeness, validate_hierarchy
    # Rules: max_missing_species, max_missing_strain, require_genus
    return results

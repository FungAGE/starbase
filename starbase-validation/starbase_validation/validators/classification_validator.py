"""Classification validation (family, navis, haplotype)."""


def validate_classification(records, rules):
    """Check family/navis/haplotype classification. Return ValidationResults-compatible list."""
    from starbase_validation.validators.results import ValidationResults

    results = ValidationResults()
    # TODO: Extract from fill_family_ids classification logic
    # Rules: max_missing_family (0% = blocking)
    return results

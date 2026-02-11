"""Validators for taxonomy, genome, classification, sequences, and relationships."""

from starbase_validation.validators.classification_validator import validate_classification
from starbase_validation.validators.genome_validator import validate_genome
from starbase_validation.validators.relationship_validator import validate_relationships
from starbase_validation.validators.sequence_validator import validate_sequences
from starbase_validation.validators.taxonomy_validator import validate_taxonomy


def validate_all(records, quality_rules):
    """Run all validators and aggregate results."""
    from starbase_validation.validators.results import ValidationResults

    gates = quality_rules.get("quality_gates", quality_rules)
    results = ValidationResults()
    results.add(validate_taxonomy(records, gates.get("taxonomy", {})))
    results.add(validate_genome(records, gates.get("genome", {})))
    results.add(validate_classification(records, gates.get("classification", {})))
    results.add(validate_sequences(records, gates.get("sequences", {})))
    results.add(validate_relationships(records))
    results.evaluate_thresholds(quality_rules)
    return results

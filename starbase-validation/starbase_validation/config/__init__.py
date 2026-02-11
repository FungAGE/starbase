"""Configuration for quality rules and enrichment sources."""

import os

_CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))


def get_quality_rules_path():
    """Return path to quality_rules.yaml."""
    return os.path.join(_CONFIG_DIR, "quality_rules.yaml")


def get_enrichment_sources_path():
    """Return path to enrichment_sources.yaml."""
    return os.path.join(_CONFIG_DIR, "enrichment_sources.yaml")

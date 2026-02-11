"""Enrich taxonomy metadata from external sources (e.g. NCBI)."""


def enrich_taxonomy(records, enrichment_config=None):
    """Fill missing taxonomy from NCBI or other sources. Returns enriched records (in-place or new structure)."""
    # TODO: Use enrichment_sources.yaml (NCBI), rate limiting, optional API key.
    return records


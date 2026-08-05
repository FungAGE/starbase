"""Shared whitelist/config for admin grid tables.

Single source of truth for src/pages/admin.py (UI) and
src/database/admin_manager_impl.py (backend-side write validation) --
the whitelist must be enforced on whichever side actually executes SQL,
so both must see the same EDITABLE_COLS/TABLE_CONFIG.
"""

EDITABLE_COLS = {
    "joined_ships": {"curated_status", "evidence", "source", "tax_id", "accession_id"},
    "submissions": {
        "needs_review",
        "classification_family",
        "classification_navis",
        "classification_haplotype",
        "classification_confidence",
        "comment",
    },
    "taxonomy": {
        "name",
        "taxID",
        "genus",
        "species",
        "strain",
        "superkingdom",
        "clade",
        "kingdom",
        "subkingdom",
        "phylum",
        "subphylum",
        "class",
        "subclass",
        "order",
        "suborder",
        "family",
        "section",
        "species_group",
        "subgenus",
    },
    "papers": {
        "Title",
        "Author",
        "PublicationYear",
        "DOI",
        "Url",
        "shortCitation",
        "starshipMentioned",
        "typePaper",
    },
    "family_names": {
        "familyName",
        "notes",
        "type_element_reference",
        "clade",
        "longFamilyID",
        "oldFamilyID",
        "otherFamilyID",
    },
    "navis_names": {"navis_name", "previous_navis_name", "activity"},
    "haplotype_names": {"haplotype_name", "previous_haplotype_name", "activity"},
    "accessions": {"version_tag"},
    "ship_accessions": {"ship_accession_display", "ship_version_tag"},
    "genomes": {
        "ome",
        "version",
        "genomeSource",
        "citation",
        "biosample",
        "acquisition_date",
        "assembly_accession",
        "taxonomy_id",
    },
    "ships": {"md5", "rev_comp_md5"},
    "starship_features": {
        "contigID",
        "elementBegin",
        "elementEnd",
        "elementLength",
        "strand",
    },
}

# One DB (starbase.sqlite) for every table now -- no per-table "db" key needed.
TABLE_CONFIG = {
    "joined_ships": {"sql_table": "joined_ships", "pk": "id"},
    "submissions": {"sql_table": "submissions", "pk": "id"},
    "taxonomy": {"sql_table": "taxonomy", "pk": "id"},
    "papers": {"sql_table": "papers", "pk": "id"},
    "family_names": {"sql_table": "family_names", "pk": "id"},
    "navis_names": {"sql_table": "navis_names", "pk": "id"},
    "haplotype_names": {"sql_table": "haplotype_names", "pk": "id"},
    "accessions": {"sql_table": "accessions", "pk": "id"},
    "ship_accessions": {"sql_table": "ship_accessions", "pk": "id"},
    "genomes": {"sql_table": "genomes", "pk": "id"},
    "ships": {"sql_table": "ships", "pk": "id"},
    "starship_features": {"sql_table": "starship_features", "pk": "id"},
}

# SQLite reserved / quoted column names in UPDATE statements
SQL_QUOTED_COLS = {"class", "order", "source"}

# Columns that should be stored as integers (SQLite booleans)
BOOLEAN_COLS = {"needs_review"}

# Grid IDs keyed by table_key -- used for batch save/discard
GRID_IDS = {
    "joined_ships": "admin-joined-ships-grid",
    "submissions": "admin-submissions-grid",
    "taxonomy": "admin-taxonomy-grid",
    "papers": "admin-papers-grid",
    "family_names": "admin-family-names-grid",
    "navis_names": "admin-navis-names-grid",
    "haplotype_names": "admin-haplotype-names-grid",
    "accessions": "admin-accessions-grid",
    "ship_accessions": "admin-ship-accessions-grid",
    "genomes": "admin-genomes-grid",
}


def sql_col_ref(col_id):
    if col_id in SQL_QUOTED_COLS:
        return f"`{col_id}`"
    return col_id

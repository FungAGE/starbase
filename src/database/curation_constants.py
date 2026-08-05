"""Shared constants for the curation/annotation-review workflow.

Single source of truth for src/pages/curation.py (UI) and
curation_manager_impl.py (backend) -- matches the admin_table_config.py
pattern used for the admin grid whitelist.
"""

# Flag values ported from MAS4starships' starship.Annotation.flag choices.
# No 6; 8/9 exist as undocumented legacy sentinels in the source app and are
# intentionally not carried over.
FLAG_LABELS = {
    0: "GREEN",
    1: "YELLOW",
    2: "RED",
    3: "REVIEW_NAME",
    4: "N_A",
    5: "ORANGE",
    7: "UNANNOTATED",
}

FLAG_COLORS = {
    0: "green",
    1: "yellow",
    2: "red",
    3: "pink",
    4: "gray",
    5: "orange",
    7: "gray",
}

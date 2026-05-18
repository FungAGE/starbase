import os

# Get the project root directory (where the app runs from)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

DATA_DIR = os.environ.get("DATA_DIR", "")

# Create required directories
REQUIRED_DIRS = [
    os.path.join(DATA_DIR),
    os.path.join(PROJECT_ROOT, "src", "database", "cache"),
    os.path.join(DATA_DIR, "ships", "fna", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "fna", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "faa", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "fna", "hmm"),
    os.path.join(DATA_DIR, "captain", "tyr", "faa", "hmm"),
    os.path.join(DATA_DIR, "ships", "gbks"),
]

for directory in REQUIRED_DIRS:
    os.makedirs(directory, exist_ok=True)

# Database paths
DB_PATHS = {
    "starbase": os.path.join(DATA_DIR, "starbase.sqlite"),
    "submissions": os.path.join(DATA_DIR, "submissions.sqlite"),
    "telemetry": os.path.join(DATA_DIR, "telemetry.sqlite"),
}

# BLAST database paths
BLAST_DB_PATHS = {
    "ship": {
        "all": {
            "nucl": os.path.join(DATA_DIR, "ships", "fna", "blastdb", "ships_all.fa")
        },
        "curated": {
            "nucl": os.path.join(
                DATA_DIR, "ships", "fna", "blastdb", "ships_curated.fa"
            )
        },
    },
    "gene": {
        "tyr": {
            "nucl": os.path.join(
                DATA_DIR, "captain", "tyr", "fna", "blastdb", "captains.fna"
            ),
            "prot": os.path.join(
                DATA_DIR, "captain", "tyr", "faa", "blastdb", "captains.faa"
            ),
            "hmm": {
                "nucl": os.path.join(
                    DATA_DIR, "captain", "tyr", "fna", "hmm", "combined.hmm"
                ),
                "prot": os.path.join(
                    DATA_DIR, "captain", "tyr", "faa", "hmm", "combined.hmm"
                ),
            },
        },
    },
}

# GBK paths
GBK_PATH = os.path.join(DATA_DIR, "ships", "gbks")

# Phylogeny paths
PHYLOGENY_PATHS = {
    "tree": os.path.join(
        DATA_DIR,
        "captain",
        "tyr",
        "faa",
        "tree",
        "funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile",
    ),
    "msa": os.path.join(
        DATA_DIR,
        "captain",
        "tyr",
        "faa",
        "alignments",
        "funTyr50_cap25_crp3_p1-512_activeFilt.clipkit",
    ),
    "clades": os.path.join(
        DATA_DIR, "captain", "tyr", "faa", "tree", "superfam-clades.tsv"
    ),
}

# URL routes
HOME_URL = os.getenv("HOME_URL", "/")
WIKI_URL = os.getenv("WIKI_URL", "/wiki")
BLAST_URL = os.getenv("BLAST_URL", "/blast")
ABOUT_URL = os.getenv("ABOUT_URL", "/about")
SYNTENY_URL = os.getenv("SYNTENY_URL", "/synteny")
SUBMIT_URL = os.getenv("SUBMIT_URL", "/submit")
METRICS_URL = os.getenv("METRICS_URL", "/metrics")

# API Keys
IPSTACK_API_KEY = os.environ.get("IPSTACK_API_KEY")
MAINTENANCE_TOKEN = os.environ.get("MAINTENANCE_TOKEN")

# Backend API (compute split)
# BACKEND_API_URL: base URL of the FastAPI backend, e.g. http://100.x.y.z:8001 (Tailscale)
#                  or http://backend:8001 (local docker-compose).
#                  When set, sql_manager routes all DB/BLAST calls over HTTP.
BACKEND_API_URL = os.environ.get("BACKEND_API_URL", "")
BACKEND_API_KEY = os.environ.get("BACKEND_API_KEY", "")

# GenBank files path
GBK_PATH = os.path.join(DATA_DIR, "ships", "gbks")

# Cache settings
CACHE_TIMEOUT = (
    None if os.getenv("CACHE_TIMEOUT") is None else int(os.getenv("CACHE_TIMEOUT"))
)
CACHE_DIR = os.path.join(PROJECT_ROOT, "src", "database", "cache")

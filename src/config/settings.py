import os

from dotenv import load_dotenv

# Get the project root directory (where the app runs from)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
_DEFAULT_DB_DIR = os.path.join(PROJECT_ROOT, "src", "database", "db")

load_dotenv(os.path.join(PROJECT_ROOT, ".env"))

# Development mode
IS_DEV = os.getenv("DEV_MODE", "false").lower() == "true"

# Backend API (compute split) — read early; drives path layout below.
# BACKEND_API_URL: e.g. http://100.x.y.z:8001 (Tailscale) or http://backend:8001 (compose).
BACKEND_API_URL = os.environ.get("BACKEND_API_URL", "")
BACKEND_API_KEY = os.environ.get("BACKEND_API_KEY", "")

# DATA_DIR: SQLite + BLAST data on the compute backend (institute machine).
# FRONTEND_DATA_DIR: lightweight local DBs (submissions, telemetry) on the Serve pod / dev host.
DATA_DIR = os.environ.get("DATA_DIR", "")
FRONTEND_DATA_DIR = os.environ.get("FRONTEND_DATA_DIR") or _DEFAULT_DB_DIR
_local_db_root = FRONTEND_DATA_DIR if BACKEND_API_URL else (DATA_DIR or _DEFAULT_DB_DIR)
_compute_data_root = DATA_DIR or _DEFAULT_DB_DIR

# Create directories: full tree on backend/monolith; only local DB dir on split frontend.
if BACKEND_API_URL:
    os.makedirs(_local_db_root, exist_ok=True)
    os.makedirs(os.path.join(PROJECT_ROOT, "src", "database", "cache"), exist_ok=True)
elif DATA_DIR:
    REQUIRED_DIRS = [
        DATA_DIR,
        os.path.join(PROJECT_ROOT, "src", "database", "cache"),
        os.path.join(DATA_DIR, "ships", "fna", "blastdb"),
        os.path.join(DATA_DIR, "captain", "tyr", "fna", "blastdb"),
        os.path.join(DATA_DIR, "captain", "tyr", "faa", "blastdb"),
        os.path.join(DATA_DIR, "captain", "tyr", "fna", "hmm"),
        os.path.join(DATA_DIR, "captain", "tyr", "faa", "hmm"),
        os.path.join(DATA_DIR, "ships", "gbks"),
    ]
    for _dir in REQUIRED_DIRS:
        os.makedirs(_dir, exist_ok=True)

# Database paths
# "submissions" now lives in the same file as "starbase" -- Submission is on
# the same SQLAlchemy Base as everything else, it was only ever split into a
# separate physical file via a separate engine/session. One DB, one
# transaction (the promote workflow used to write two files non-atomically).
# telemetry stays separate: raw CREATE TABLE, not on Base.metadata, high
# write volume, no need to join against curator data.
_starbase_db_path = os.path.join(_compute_data_root, "starbase.sqlite")
DB_PATHS = {
    "starbase": _starbase_db_path,
    "submissions": _starbase_db_path,
    "telemetry": os.path.join(_local_db_root, "telemetry.sqlite"),
}

# BLAST database paths (backend / monolith only)
BLAST_DB_PATHS = {
    "ship": {
        "all": {
            "nucl": os.path.join(
                _compute_data_root, "ships", "fna", "blastdb", "ships_all.fa"
            )
        },
        "curated": {
            "nucl": os.path.join(
                _compute_data_root, "ships", "fna", "blastdb", "ships_curated.fa"
            )
        },
    },
    "gene": {
        "tyr": {
            "nucl": os.path.join(
                _compute_data_root, "captain", "tyr", "fna", "blastdb", "captains.fna"
            ),
            "prot": os.path.join(
                _compute_data_root, "captain", "tyr", "faa", "blastdb", "captains.faa"
            ),
            "hmm": {
                "nucl": os.path.join(
                    _compute_data_root,
                    "captain",
                    "tyr",
                    "fna",
                    "hmm",
                    "combined.hmm",
                ),
                "prot": os.path.join(
                    _compute_data_root,
                    "captain",
                    "tyr",
                    "faa",
                    "hmm",
                    "combined.hmm",
                ),
            },
        },
    },
}

GBK_PATH = os.path.join(_compute_data_root, "ships", "gbks")

# Phylogeny paths
PHYLOGENY_PATHS = {
    "tree": os.path.join(
        _compute_data_root,
        "captain",
        "tyr",
        "faa",
        "tree",
        "funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.treefile",
    ),
    "msa": os.path.join(
        _compute_data_root,
        "captain",
        "tyr",
        "faa",
        "alignments",
        "funTyr50_cap25_crp3_p1-512_activeFilt.clipkit",
    ),
    "clades": os.path.join(
        _compute_data_root,
        "captain",
        "tyr",
        "faa",
        "tree",
        "superfam-clades.tsv",
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
ADMIN_URL = os.getenv("ADMIN_URL", "/admin")
ADMIN_TOKEN = os.environ.get("ADMIN_TOKEN")

# Define valid pages
PAGES = {
    HOME_URL,
    WIKI_URL,
    BLAST_URL,
    ABOUT_URL,
    SYNTENY_URL,
    SUBMIT_URL,
    METRICS_URL,
}

PAGE_MAPPING = ",".join(PAGES)


# API Keys
IPSTACK_API_KEY = os.environ.get("IPSTACK_API_KEY")
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")
SECRET_KEY = os.environ.get("SECRET_KEY")

# Cache settings
CACHE_TIMEOUT = (
    None if os.getenv("CACHE_TIMEOUT") is None else int(os.getenv("CACHE_TIMEOUT"))
)
CACHE_DIR = os.path.join(PROJECT_ROOT, "src", "database", "cache")

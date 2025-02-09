import os
from pathlib import Path

# Get the project root directory (where the app runs from)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

# Try potential database directories in order of preference
potential_db_dirs = [
    os.path.join(PROJECT_ROOT,"src","database", "db"),
    os.path.join(PROJECT_ROOT,"database", "db")
]

# Use the first valid directory path, or default to the last option
DATA_DIR = next((path for path in potential_db_dirs if path is not None), potential_db_dirs[-1])

# Create required directories
REQUIRED_DIRS = [
    os.path.join(DATA_DIR),
    os.path.join(DATA_DIR, "cache"),
    os.path.join(DATA_DIR, "ships", "fna", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "fna", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "faa", "blastdb"),
    os.path.join(DATA_DIR, "captain", "tyr", "fna", "hmm"),
    os.path.join(DATA_DIR, "captain", "tyr", "faa", "hmm"),
]

for directory in REQUIRED_DIRS:
    os.makedirs(directory, exist_ok=True)

# Database paths
DB_PATHS = {
    'starbase': os.path.join(DATA_DIR, 'starbase.sqlite'),
    'submissions': os.path.join(DATA_DIR, 'submissions.sqlite'),
    'telemetry': os.path.join(DATA_DIR, 'telemetry.sqlite')
}

# BLAST database paths
BLAST_DB_PATHS = {
    "ship": {"nucl": os.path.join(DATA_DIR, "ships", "fna", "blastdb", "ships.fa")},
    "gene": {
        "tyr": {
            "nucl": os.path.join(DATA_DIR, "captain", "tyr", "fna", "blastdb", "captains.fna"),
            "prot": os.path.join(DATA_DIR, "captain", "tyr", "faa", "blastdb", "captains.faa"),
            "hmm": {
                "nucl": os.path.join(DATA_DIR, "captain", "tyr", "fna", "hmm", "combined.hmm"),
                "prot": os.path.join(DATA_DIR, "captain", "tyr", "faa", "hmm", "combined.hmm"),
            },
        },
    },
}

# URL routes
HOME_URL = os.getenv("HOME_URL", "/")
WIKI_URL = os.getenv("WIKI_URL", "/wiki")
BLAST_URL = os.getenv("BLAST_URL", "/blast")
ABOUT_URL = os.getenv("ABOUT_URL", "/about")
PGV_URL = os.getenv("PGV_URL", "/pgv")
SUBMIT_URL = os.getenv("SUBMIT_URL", "/submit")
DL_URL = os.getenv("DL_URL", "/download")
METRICS_URL = os.getenv("METRICS_URL", "/metrics")

# API Keys
IPSTACK_API_KEY = os.environ.get('IPSTACK_API_KEY')
MAINTENANCE_TOKEN = os.environ.get('MAINTENANCE_TOKEN')

# Cache settings
CACHE_TIMEOUT = int(os.getenv('CACHE_TIMEOUT', 86400))
CACHE_DIR = os.getenv('CACHE_DIR', '/tmp/starbase_cache')
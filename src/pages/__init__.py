import os
from pathlib import Path
from dotenv import load_dotenv

env_path = Path('.') / '.env'
if env_path.exists():
    load_dotenv(dotenv_path=env_path)

HOME_URL = os.getenv("HOME_URL", "/")
WIKI_URL = os.getenv("WIKI_URL", "/wiki")
BLAST_URL = os.getenv("BLAST_URL", "/blast")
ABOUT_URL = os.getenv("ABOUT_URL", "/about")
PGV_URL = os.getenv("PGV_URL", "/pgv")
SUBMIT_URL = os.getenv("SUBMIT_URL", "/submit")
DL_URL = os.getenv("DL_URL", "/download")
METRICS_URL = os.getenv("METRICS_URL", "/metrics")
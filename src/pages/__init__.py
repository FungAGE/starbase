from dotenv import load_dotenv
import os

load_dotenv("env")

HOME_URL = os.getenv("HOME_URL")
WIKI_URL = os.getenv("WIKI_URL")
# EXPLORE_URL = os.getenv("EXPLORE_URL")
BLAST_URL = os.getenv("BLAST_URL")
ABOUT_URL = os.getenv("ABOUT_URL")
PGV_URL = os.getenv("PGV_URL")
SUBMIT_URL = os.getenv("SUBMIT_URL")
DL_URL = os.getenv("DL_URL")

import pandas as pd
from sqlalchemy import create_engine
from src.components.config import MOUNTED_DIRECTORY_PATH

engine = None

if MOUNTED_DIRECTORY_PATH:
    # Initialize the SQLite engine using the mounted directory path
    engine = create_engine(f"sqlite:///{MOUNTED_DIRECTORY_PATH}/starbase.sqlite")
    # print(f"Engine created for database at: {MOUNTED_DIRECTORY_PATH}")
else:
    print("Error: Mounted directory not found.")

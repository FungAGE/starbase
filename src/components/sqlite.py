import pandas as pd
from sqlalchemy import create_engine
from app import MOUNTED_DIRECTORY_PATH

if MOUNTED_DIRECTORY_PATH:
    engine = create_engine(f"sqlite:///{MOUNTED_DIRECTORY_PATH}/starbase.sqlite")
    query = "SELECT name FROM sqlite_master WHERE type='table'"
    sql_tbls = pd.read_sql_query(query, engine)

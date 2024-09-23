import pandas as pd
from sqlalchemy import create_engine

engine = create_engine("sqlite:///database_folder/starbase.sqlite")

query = "SELECT name FROM sqlite_master WHERE type='table'"
sql_tbls = pd.read_sql_query(query, engine)

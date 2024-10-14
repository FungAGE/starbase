import pandas as pd
from sqlalchemy import create_engine

engine = None

# Initialize the SQLite engine
engine = create_engine("sqlite:///src/data/db/starbase.sqlite")

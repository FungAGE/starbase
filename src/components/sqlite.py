from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from models import Base

# Connect to the SQLite database file
DATABASE_URL = "sqlite:///src/data/db/starbase.sqlite"
# Initialize the SQLite engine
engine = create_engine(DATABASE_URL)

# Bind the models (Base) to the engine
Base.metadata.bind = engine

# Create all tables in the database (if they don't exist)
Base.metadata.create_all(engine)

# Create a session factory
Session = sessionmaker(bind=engine)
session = Session()

from src.config.database import Base, starbase_engine, submissions_engine
from src.database.models.schema import *

def init_databases():
    """Initialize all database tables"""
    # Create all tables in starbase database
    Base.metadata.create_all(starbase_engine)
    
    # Create submissions tables
    Base.metadata.create_all(submissions_engine)

if __name__ == "__main__":
    init_databases() 
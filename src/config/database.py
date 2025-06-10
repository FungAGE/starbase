from src.config.settings import DB_PATHS
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import os

# Get the environment
ENV = os.getenv("ENVIRONMENT", "development")
IS_DEV = ENV.lower() == "development"


def create_db_engine(db_path):
    """Create a database engine with the given path"""
    if not os.path.exists(db_path):
        db_dir = os.path.dirname(db_path)
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

    # Set echo based on environment - only log SQL in development
    return create_engine(f"sqlite:///{db_path}", echo=IS_DEV)


# Create engines
starbase_engine = create_db_engine(DB_PATHS["starbase"])
submissions_engine = create_db_engine(DB_PATHS["submissions"])
telemetry_engine = create_db_engine(DB_PATHS["telemetry"])

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionsSession = sessionmaker(bind=submissions_engine)
TelemetrySession = sessionmaker(bind=telemetry_engine)

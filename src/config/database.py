from src.config.settings import DB_PATHS
from sqlalchemy import create_engine
from sqlalchemy.pool import QueuePool
from sqlalchemy.orm import sessionmaker
import os

# Get the environment
ENV = os.getenv("ENVIRONMENT", "development")
IS_DEV = ENV.lower() == "development"

DATABASE_URLS = {
    "starbase": f"sqlite:///{DB_PATHS['starbase']}",
    "submissions": f"sqlite:///{DB_PATHS['submissions']}",
    "telemetry": f"sqlite:///{DB_PATHS['telemetry']}",
}

# Create engines with proper connection pooling
engines = {
    name: create_engine(
        url,
        poolclass=QueuePool,
        pool_size=10,
        max_overflow=20,
        pool_timeout=30,
        pool_recycle=1800,
        echo=IS_DEV,  # Only log SQL in development mode
    )
    for name, url in DATABASE_URLS.items()
}

# Extract individual engines for backward compatibility
starbase_engine = engines["starbase"]
submissions_engine = engines["submissions"]
telemetry_engine = engines["telemetry"]

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionsSession = sessionmaker(bind=submissions_engine)
TelemetrySession = sessionmaker(bind=telemetry_engine)

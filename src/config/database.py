from src.config.settings import DB_PATHS, IS_DEV
from sqlalchemy import create_engine
from sqlalchemy.pool import QueuePool
from sqlalchemy.orm import sessionmaker

DATABASE_URLS = {
    "starbase": f"sqlite:///{DB_PATHS['starbase']}",
    "submissions": f"sqlite:///{DB_PATHS['submissions']}",
    "telemetry": f"sqlite:///{DB_PATHS['telemetry']}",
}

# Create engines with proper connection pooling
# For SQLite, we need special handling due to threading restrictions
engines = {}
for name, url in DATABASE_URLS.items():
    if url.startswith("sqlite"):
        engines[name] = create_engine(
            url,
            poolclass=QueuePool,
            pool_size=10,
            max_overflow=20,
            pool_timeout=30,
            pool_recycle=1800,
            echo=IS_DEV,
            connect_args={"check_same_thread": False},
        )
    else:
        # MySQL/PostgreSQL configuration
        engines[name] = create_engine(
            url,
            poolclass=QueuePool,
            pool_size=10,
            max_overflow=20,
            pool_timeout=30,
            pool_recycle=1800,
            echo=IS_DEV,
        )

# Extract individual engines for backward compatibility
starbase_engine = engines["starbase"]
submissions_engine = engines["submissions"]
telemetry_engine = engines["telemetry"]

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionsSession = sessionmaker(bind=submissions_engine)
TelemetrySession = sessionmaker(bind=telemetry_engine)

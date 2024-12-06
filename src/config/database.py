from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
import os

# Create base class for declarative models
Base = declarative_base()

# Database configuration
DB_PATHS = {
    'starbase': os.getenv('STARBASE_PATH', 'src/database/db/starbase.sqlite'),
    'submissions': os.getenv('SUBMISSIONS_PATH', 'src/database/db/submissions.sqlite'),
    'telemetry': os.getenv('TELEMETRY_PATH', 'src/database/db/telemetry.sqlite')
}

def create_db_engine(db_path):
    """Create a database engine with the given path"""
    if not os.path.exists(db_path):
        db_dir = os.path.dirname(db_path)
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)
    return create_engine(f'sqlite:///{db_path}', echo=False)

# Create engines
starbase_engine = create_db_engine(DB_PATHS['starbase'])
submissions_engine = create_db_engine(DB_PATHS['submissions'])
telemetry_engine = create_db_engine(DB_PATHS['telemetry'])

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionsSession = sessionmaker(bind=submissions_engine)
TelemetrySession = sessionmaker(bind=telemetry_engine)
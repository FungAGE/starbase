from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
import os

# Get the project root directory (where the app runs from)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

# Default database directory
DEFAULT_DB_DIR = os.path.join(PROJECT_ROOT, "src", "database", "db")

# Create base class for declarative models
Base = declarative_base()

# Database configuration
DB_PATHS = {
    'starbase': os.getenv('STARBASE_PATH', os.path.join(os.environ.get('STARBASE_DB_DIR', DEFAULT_DB_DIR), 'starbase.sqlite')),
    'submissions': os.getenv('SUBMISSIONS_PATH', os.path.join(os.environ.get('STARBASE_DB_DIR', DEFAULT_DB_DIR), 'submissions.sqlite')),
    'telemetry': os.getenv('TELEMETRY_PATH', os.path.join(os.environ.get('STARBASE_DB_DIR', DEFAULT_DB_DIR), 'telemetry.sqlite'))
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
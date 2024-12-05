from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
import os

# Create base class for declarative models
Base = declarative_base()

# Database configuration
DB_PATH = os.getenv('DB_PATH', 'src/database/db/starbase.sqlite')
SUBMISSIONS_PATH = os.getenv('SUBMISSIONS_PATH', 'src/database/db/submissions.sqlite')

# Create engines
starbase_engine = create_engine(f'sqlite:///{DB_PATH}', echo=True)
submissions_engine = create_engine(f'sqlite:///{SUBMISSIONS_PATH}', echo=True)

# Create session factories
StarbaseSession = sessionmaker(bind=starbase_engine)
SubmissionsSession = sessionmaker(bind=submissions_engine)
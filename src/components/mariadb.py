import os
import logging
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from pymysql.err import OperationalError
from models import Base
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

# Load the environment variables from the .env file
load_dotenv(dotenv_path=".env", override=True)

db_user = os.getenv("DB_USER")
db_password = os.getenv("DB_PASSWORD")
db_host = os.getenv("DB_HOST")
db_port = os.getenv("DB_PORT", "3307")
db_name = os.getenv("DB_NAME")

connection_str = (
    f"mysql+pymysql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}"
)

sql_connected = False  # Default to False before trying to connect

try:
    # Ensure the old engine is deleted from memory
    if "engine" in globals():
        del engine  # Remove reference to the existing engine

    # Attempt to connect to the SQL database
    engine = create_engine(
        connection_str,
        pool_pre_ping=True,
        pool_size=5,
        max_overflow=10,
        pool_recycle=1800,
        pool_timeout=30,
    )
    sql_connected = True
    logger.info("Successfully connected to the SQL database.")

    # Create session
    Session = sessionmaker(bind=engine)
    session = Session()
    logger.info("Session factory created and session started.")

except OperationalError as e:
    logger.error("Could not connect to SQL server: %s", e)
except Exception as e:
    logger.exception(
        "An unexpected error occurred while trying to connect to the SQL server."
    )

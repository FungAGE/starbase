import os
import logging
from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from models import Base

logger = logging.getLogger(__name__)

db_user = os.getenv("DB_USER")
db_password = os.getenv("DB_PASSWORD")
db_host = os.getenv("DB_HOST")
db_port = os.getenv("DB_PORT", "3307")
db_name = os.getenv("DB_NAME")

if db_user and db_password and db_host and db_name:
    logger.debug("Successfully fetched DB credentials from environment variables.")
else:
    logger.error(
        "Failed to fetch one or more DB credentials from environment variables."
    )

connection_str = (
    f"mysql+pymysql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}"
)
logger.debug(f"Constructed connection string: {connection_str}")

try:
    engine = create_engine(connection_str)
    logger.info(
        f"Successfully created the engine for database: {db_name} at {db_host}:{db_port}"
    )
except Exception as e:
    logger.exception("Error creating the engine.")
    raise e

metadata = MetaData()

Base = declarative_base()

Base.metadata.bind = engine
logger.debug("Bound the Base metadata to the engine.")

try:
    Session = sessionmaker(bind=engine)
    session = Session()
    logger.info("Session factory created and session started.")
except Exception as e:
    logger.exception("Failed to create session or start session.")
    raise e

try:
    Base.metadata.create_all(engine)
    logger.info("All tables created (if they did not already exist).")
except Exception as e:
    logger.exception("Error creating tables in the database.")
    raise e

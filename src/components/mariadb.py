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
# override=True will overwrite existing variables
load_dotenv(dotenv_path=".env", override=True)

db_user = os.getenv("DB_USER")
db_password = os.getenv("DB_PASSWORD")
db_host = os.getenv("DB_HOST")
db_port = os.getenv("DB_PORT", "3307")
db_name = os.getenv("DB_NAME")


connection_str = "mysql+pymysql://%s:%s@%s:%s/%s" % (
    db_user,
    db_password,
    db_host,
    db_port,
    db_name,
)

try:
    # Ensure the old engine is deleted from memory
    if "engine" in globals():
        del engine  # Remove reference to the existing engine
        # engine.dispose()  # Clears out any stale connections
except Exception as e:
    logger.debug("No previous engine to dispose of, or disposal failed.")

# Attempt to connect to the SQL database
try:
    engine = create_engine(
        connection_str,
        pool_pre_ping=True,
        pool_size=5,
        max_overflow=10,
        pool_recycle=1800,
        pool_timeout=30,
        # echo=True,
    )

    sql_connected = True
except OperationalError as e:
    print("Could not connect to SQL server:", e)
    sql_connected = False


# metadata = MetaData()

# Base = declarative_base()

# Base.metadata.bind = engine
# logger.debug("Bound the Base metadata to the engine.")

try:
    Session = sessionmaker(bind=engine)
    session = Session()
    logger.info("Session factory created and session started.")
except Exception as e:
    logger.exception("Failed to create session or start session.")
    raise e

# try:
#     Base.metadata.create_all(engine)
#     logger.info("All tables created (if they did not already exist).")
# except Exception as e:
#     logger.exception("Error creating tables in the database.")
#     raise e

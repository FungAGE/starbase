from sqlalchemy import create_engine, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from models import Base
import os

# Fetch credentials from environment variables
db_user = os.getenv("DB_USER")
db_password = os.getenv("DB_PASSWORD")
db_host = os.getenv("DB_HOST")
db_port = os.getenv("DB_PORT", "3307")
db_name = os.getenv("DB_NAME")


connection_str = (
    f"mysql+pymysql://{db_user}:{db_password}@{db_host}:{db_port}/{db_name}"
)
engine = create_engine(connection_str)

# Create metadata object to handle table creation
metadata = MetaData()

Base = declarative_base()

# Bind the models (Base) to the engine
Base.metadata.bind = engine

# Create a session factory
Session = sessionmaker(bind=engine)
session = Session()

# Create all tables in the database (if they don't exist)
Base.metadata.create_all(engine)

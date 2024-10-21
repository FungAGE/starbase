import sqlite3
from datetime import datetime

import warnings

import datetime
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String, Text
from sqlalchemy.exc import NoSuchTableError

warnings.filterwarnings("ignore")

db_url = "sqlite:///telemetry.sqlite"


def init_db(db_url):
    engine = create_engine(db_url)
    metadata = MetaData()
    requests_table = Table(
        "request_logs",
        metadata,
        Column("id", Integer, primary_key=True, autoincrement=True),
        Column("endpoint", Text, nullable=False),
        Column("timestamp", Text, nullable=False),
    )

    try:
        requests_table = Table("request_logs", metadata, autoload_with=engine)
    except NoSuchTableError:
        metadata.create_all(engine)

    return engine


engine = init_db(db_url)


# Connect to your telemetry database
def add_request():
    with engine.connect() as connection:
        connection.execute(
            """
          CREATE TABLE IF NOT EXISTS request_logs (
            id INTEGER PRIMARY KEY,
            ip TEXT,
            endpoint TEXT,
            timestamp TEXT
              )
          """
        )


def log_request(ip, endpoint):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with engine.connect() as connection:
        connection.execute(
            """INSERT INTO request_logs (ip, endpoint, timestamp) VALUES (?, ?, ?)""",
            (ip, endpoint, timestamp),
        )

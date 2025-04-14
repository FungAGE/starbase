from src.config.database import (
    Base,
    starbase_engine,
    submissions_engine,
    telemetry_engine,
)


def init_databases():
    """Initialize all database tables"""
    Base.metadata.create_all(starbase_engine)

    Base.metadata.create_all(submissions_engine)

    Base.metadata.create_all(telemetry_engine)


if __name__ == "__main__":
    init_databases()

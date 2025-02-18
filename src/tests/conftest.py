import pytest
from src.config.database import StarbaseSession

@pytest.fixture(scope="function")
def db_session():
    """Provide a database session for tests"""
    session = StarbaseSession
    yield session
    session.rollback()
    session.close()
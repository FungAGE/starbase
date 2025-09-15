"""
Updated conftest.py with mixed test data strategy integration.
"""

import pytest
import pandas as pd
from src.database.models.schema import Base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from flask import jsonify

from flask import Flask
from src.config.cache import cache
from src.config.limiter import limiter
from src.api import register_routes

from werkzeug.middleware.proxy_fix import ProxyFix
from flask_compress import Compress
from src.config.cache import cleanup_old_cache

from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium import webdriver

# Import our new test data management
from tests.utils.test_data_manager import TestDataManager
from tests.test_config import test_config


@pytest.fixture(scope="session")
def test_client():
    """Create a test client for the Flask server."""
    server = Flask(__name__)
    server.wsgi_app = ProxyFix(server.wsgi_app, x_for=1, x_proto=1)
    Compress(server)

    server.config.update(
        TESTING=True,
        MAX_CONTENT_LENGTH=10 * 1024 * 1024,  # 10MB limit
        CACHE_TYPE="SimpleCache",
        CACHE_DEFAULT_TIMEOUT=300,
        SEND_FILE_MAX_AGE_DEFAULT=0,
        COMPRESS_MIMETYPES=["text/html", "text/css", "application/javascript"],
        COMPRESS_LEVEL=6,
        COMPRESS_ALGORITHM=["gzip", "br"],
    )

    cache.init_app(server)
    cleanup_old_cache()
    limiter.init_app(server)

    register_routes(server, limiter)
    limiter.enabled = True
    limiter.limit("10 per hour")(server.view_functions["blast.check_blast_limit"])

    with server.test_client() as testing_client:
        with server.app_context():
            yield testing_client


@pytest.fixture(scope="function")
def db_session():
    """Provide a database session for tests."""
    # Implementation depends on your database setup
    # This is a placeholder
    pass


# Test Data Manager Fixtures
@pytest.fixture
def test_data_manager():
    """Fixture providing a test data manager with mock data."""
    return TestDataManager(use_real_data=False)


@pytest.fixture
def real_data_manager():
    """Fixture providing a test data manager with real data from database."""
    return TestDataManager(use_real_data=True)


# Mock Data Fixtures
@pytest.fixture
def mock_ships_df(test_data_manager):
    """Fixture providing mock ships data."""
    return test_data_manager.get_ships_data(with_sequence=True)


@pytest.fixture
def mock_meta_df(test_data_manager):
    """Fixture providing mock metadata."""
    return test_data_manager.get_meta_data()


@pytest.fixture
def mock_sequence(test_data_manager):
    """Fixture providing a test sequence."""
    return test_data_manager.get_test_sequence(length=1000)


# Real Data Fixtures (conditional on database availability)
@pytest.fixture
def real_ships_df(real_data_manager):
    """Fixture providing real ships data from database."""
    try:
        return real_data_manager.get_ships_data(with_sequence=True, limit=10)
    except Exception:
        pytest.skip("Real ships data not available")


@pytest.fixture
def real_meta_df(real_data_manager):
    """Fixture providing real metadata from database."""
    try:
        return real_data_manager.get_meta_data(limit=10)
    except Exception:
        pytest.skip("Real metadata not available")


# Configuration Fixtures
@pytest.fixture
def config():
    """Fixture providing test configuration."""
    return test_config


@pytest.fixture
def use_real_data_flag(config):
    """Fixture providing real data usage flag."""
    return config.should_use_real_data()


@pytest.fixture
def database_available(config):
    """Fixture indicating if database is available."""
    return config.is_database_available()


# Legacy fixtures for backward compatibility
@pytest.fixture
def sample_meta_dataframe():
    """Legacy fixture - use mock_meta_df instead."""
    return pd.DataFrame({
        "accession_tag": ["SBS000001", "SBS000002", "SBS000003"],
        "familyName": ["Voyager", "Voyager", "Voyager"],
        "name": ["Test Organism 1", "Test Organism 2", "Test Organism 3"],
        "curated_status": ["curated", "curated", "uncurated"],
        "elementLength": [1000, 1000, 1000],
        "upDR": [100, 100, 100],
        "downDR": [100, 100, 100],
        "contigID": ["contig1", "contig2", "contig3"],
        "captainID": ["captain1", "captain2", "captain3"],
        "elementBegin": [100, 200, 300],
        "elementEnd": [1100, 1200, 1300],
        "type_element_reference": ["ref1", "ref2", "ref3"],
        "navis_name": ["navis1", "navis2", "navis3"],
        "haplotype_name": ["haplo1", "haplo2", "haplo3"],
        "ome": ["ome1", "ome2", "ome3"],
        "version": ["1.0", "1.0", "1.0"],
        "genomeSource": ["source1", "source2", "source3"],
        "citation": ["citation1", "citation2", "citation3"],
        "assembly_accession": ["assembly1", "assembly2", "assembly3"]
    })


@pytest.fixture
def sample_papers_dataframe():
    """Legacy fixture for papers data."""
    return pd.DataFrame({
        "accession_tag": ["SBS000001", "SBS000002", "SBS000003"],
        "title": ["Paper 1", "Paper 2", "Paper 3"],
        "authors": ["Author 1", "Author 2", "Author 3"],
        "journal": ["Journal 1", "Journal 2", "Journal 3"],
        "year": [2020, 2021, 2022],
        "doi": ["10.1000/1", "10.1000/2", "10.1000/3"]
    })


# Utility functions for tests
def skip_if_no_real_data(config, test_name: str = None):
    """Skip test if real data is not available."""
    if not config.should_use_real_data(test_name):
        pytest.skip("Real data not enabled")
    
    if not config.is_database_available():
        pytest.skip("Database not available")


def skip_if_no_mock_data():
    """Skip test if mock data is not available (should rarely happen)."""
    # Mock data should always be available, but this is here for completeness
    pass


# Pytest configuration
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "real_data: mark test as requiring real data from database"
    )
    config.addinivalue_line(
        "markers", "mock_data: mark test as using mock data only"
    )
    config.addinivalue_line(
        "markers", "blast_pipeline: mark test as part of BLAST pipeline"
    )
    config.addinivalue_line(
        "markers", "classification_pipeline: mark test as part of classification pipeline"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )
    config.addinivalue_line(
        "markers", "unit: mark test as unit test"
    )


# Pytest collection hooks
def pytest_collection_modifyitems(config, items):
    """Modify test collection based on configuration."""
    if not test_config.use_real_data:
        # Skip real data tests if not enabled
        skip_real_data = pytest.mark.skip(reason="Real data tests not enabled")
        for item in items:
            if "real_data" in item.keywords:
                item.add_marker(skip_real_data)
    
    if not test_config.is_database_available():
        # Skip tests that require database
        skip_database = pytest.mark.skip(reason="Database not available")
        for item in items:
            if "real_data" in item.keywords or "integration" in item.keywords:
                item.add_marker(skip_database)

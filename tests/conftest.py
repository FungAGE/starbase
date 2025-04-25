import pytest
import pandas as pd
from src.database.models import Base
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
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    try:
        yield session
    finally:
        session.rollback()
        session.close()


@pytest.fixture
def chrome_options():
    """Chrome options for the webdriver"""
    options = webdriver.ChromeOptions()
    options.add_argument("--headless")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    return options


@pytest.fixture
def driver(chrome_options):
    """Create a Chrome webdriver"""
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)
    yield driver
    driver.quit()


@pytest.fixture(scope="function")
def mock_accession_data():
    """Fixture to provide mock accession data as a DataFrame."""
    return pd.DataFrame(
        {
            "accession_tag": ["ABC123", "DEF456"],
            "starshipID": ["123", "456"],
            "curated_status": ["curated", "uncurated"],
            "familyName": ["Family1", "Family2"],
            "starship_navis": ["Navis1", "Navis2"],
            "starship_haplotype": ["Haplo1", "Haplo2"],
            "order": ["Order1", "Order2"],
            "family": ["Family1", "Family2"],
            "name": ["Species1", "Species2"],
            "strain": ["Strain1", "Strain2"],
            "taxID": [12345, 67890],
            "assembly_accession": ["GCA_000001", "GCA_000002"],
            "genomeSource": ["Source1", "Source2"],
            "contigID": ["Contig1", "Contig2"],
            "elementBegin": [100, 200],
            "elementEnd": [500, 600],
            "size": [400, 400],
        }
    )


@pytest.fixture(scope="function")
def simple_accession_data():
    yield jsonify({"content": "content", "title": "title"})

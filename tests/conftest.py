import pytest
import pandas as pd
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

from src.database.models.schema import Base
from src.database.sql_manager import fetch_ships, fetch_meta_data
from src.utils.seq_utils import clean_sequence
from src.utils.classification_utils import generate_md5_hash

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
    options.add_argument("--disable-gpu")
    options.add_argument("--window-size=1920,1080")
    options.add_argument("--disable-extensions")
    options.add_argument("--remote-debugging-port=9222")
    return options


@pytest.fixture
def driver(chrome_options):
    """Create a Chrome webdriver"""
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)
    yield driver
    driver.quit()



@pytest.fixture
def test_meta_data():
    try:
        meta_df = fetch_meta_data(accession_tags=["SBS000001", "SBS000002"])
        if not meta_df.empty:
            return meta_df
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return pd.DataFrame(
        {
            "accession_tag": ["ABC123", "DEF456"],
            "starshipID": ["123", "456"],
            "curated_status": ["curated", "uncurated"],
            "familyName": ["Family1", "Family2"],
            "navis_name": ["Navis1", "Navis2"],
            "haplotype_name": ["Haplo1", "Haplo2"],
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




@pytest.fixture
def test_ships_df():
    try:
        ships_df = fetch_ships(accession_tags=["SBS000001", "SBS000002"], with_sequence=True)
        if not ships_df.empty:
            return ships_df
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return pd.DataFrame(
        {
            "accession_tag": ["SBS000001", "SBS000002"],
            "sequence": [
                "ATGCATGCATGC",  # Simple sequence for exact match
                "ATGCATGCATGCATGC",  # Longer sequence for contained match
                "ATGCATGCATTT",  # Similar sequence for similarity match
            ],
            "md5": [
                
                generate_md5_hash(clean_sequence("ATGCATGCATGC")),
                generate_md5_hash(clean_sequence("ATGCATGCATGCATGC")), 
                generate_md5_hash(clean_sequence("ATGCATGCATTT")),
            ],
            "rev_comp_md5": [
                generate_md5_hash(clean_sequence("GCATGCATGCAT")),
                generate_md5_hash(clean_sequence("GCATGCATGCATGCAT")),
                generate_md5_hash(clean_sequence("AAATGCATGCAT")),
            ]
        }
    )


@pytest.fixture
def test_sequence():
    # looking for a real sequence from the database
    try:
        ship_df = fetch_ships(accession_tags=["SBS000002"], with_sequence=True)
        if not ship_df.empty:
            sequence = ship_df.iloc[0]["sequence"]
            return sequence.upper()  # Convert to uppercase to match stored MD5
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return "ATGCATGCATGC"  # Mock sequence

@pytest.fixture
def test_sequence_revcomp(test_sequence):
    try:
        # return the reverse complement of the real sequence
        complement = str.maketrans("ATGC", "TACG")
        revcomp = test_sequence.translate(complement)[::-1]
        return revcomp
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return "GCATGCATGCAT"  # Mock sequence

@pytest.fixture
def test_contained_sequence(test_sequence):
    try:
        # return a contained subsequence of the real sequence
        return test_sequence[50:-50]  # Take from position 50 to 50 positions from the end
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return "ATGCATGC"  # Mock sequence

@pytest.fixture
def test_similar_sequence(test_sequence):
    # introduce a small mutation to create a similar sequence
    try:
        return test_sequence[:10] + "A" + test_sequence[11:]  # Change one base
    except Exception:
        pass
    
    # Fallback to mock data if real data isn't available
    return "ATGCATGCAA"  # Mock sequence

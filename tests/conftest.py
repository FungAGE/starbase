import pytest
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

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
def single_accession_meta_data():
    return {
        "curated_status": "curated",
        "starshipID": "Prometheus_Ssp_2",
        "ship_id": 331,
        "joined_ship_id": 599,
        "ship_accession_tag": "SSB0000313",
        "version_tag": 1,
        "taxID": 2044277,
        "strain": "SCSIO F190",
        "order": "Onygenales",
        "family": None,
        "name": "Spiromastix sp. SCSIO F190",
        "elementLength": None,
        "upDR": None,
        "downDR": None,
        "contigID": None,
        "captainID": None,
        "elementBegin": None,
        "elementEnd": None,
        "familyName": "Prometheus",
        "type_element_reference": "Urquhart et al., unpublished",
        "navis_name": "Prometheus",
        "haplotype_name": "2",
        "ome": None,
        "version": None,
        "genomeSource": None,
        "citation": None,
        "assembly_accession": None,
        "md5": "23cde5cdde7fa335c9cad900c1275ddd",
        "rev_comp_md5": "cdcbf2c7875f729700851711baa9222f",
        "sequence_length": 108263,
    }


@pytest.fixture
def multiple_accession_meta_data():
    return {
        "curated_status": [None, "curated"],
        "starshipID": ["fususs2_s06540", "08-36-03-25_s00073"],
        "ship_id": [630, 754],
        "joined_ship_id": [418, 754],
        "ship_accession_tag": ["SSB0000566", "SSB0000676"],
        "version_tag": [1, 1],
        "taxID": [None, 746128],
        "strain": [None, "Aspergillus fumigatus"],
        "order": ["Hypocreales", "Eurotiales"],
        "family": ["Nectriaceae", "Aspergillaceae"],
        "name": ["Fusarium ussurianum", "Aspergillus fumigatus"],
        "elementLength": [None, None],
        "upDR": [None, None],
        "downDR": [None, None],
        "contigID": [None, None],
        "captainID": [None, None],
        "elementBegin": [None, None],
        "elementEnd": [None, None],
        "familyName": ["Arwing", "Prometheus"],
        "type_element_reference": [
            "Gluck-Thaler et al., 2022",
            "Urquhart et al., unpublished",
        ],
        "navis_name": [None, "navis06"],
        "haplotype_name": [None, "var03"],
        "ome": ["fususs2", "08-36-03-25"],
        "version": [None, None],
        "genomeSource": ["ncbi", None],
        "citation": [None, "10.1128/mBio.00536-15"],
        "assembly_accession": ["GCA_017656605.1", "10.5281/zenodo.5775265"],
        "md5": ["743b3591a2a317d9780ca4eced108b60", "d598a244e33d2131699489fdaac7f701"],
        "rev_comp_md5": [
            "a76457e834456d250e42cc33edcf81ca",
            "aaa598945706fbe41f4312397492a637",
        ],
        "sequence_length": [18485, 76116],
    }


@pytest.fixture
def single_accession_multiple_ships_meta_data():
    return {
        "curated_status": [
            None,
            "curated",
            "curated",
            "auto_created",
            "auto_created",
            "uncurated",
            "uncurated",
            "uncurated",
            "uncurated",
        ],
        "starshipID": [
            "aspfum6_s01218",
            "CEA10-lr_s00107",
            "CM2733_s04518",
            "aspfum137_e02579",
            "CEA10_s00026",
            "fillin_00420",
            "fillin_00678",
            "fillin_00711",
            "fillin_00765",
        ],
        "ship_id": [189, 2556, 2134, 1529, 2098, 3131, 3324, 3355, 3391],
        "joined_ship_id": [130, 768, 769, 8057, 8529, 9376, 9569, 9600, 9636],
        "ship_accession_tag": [
            "SSB0000182",
            "SSB0002279",
            "SSB0001884",
            "SSB0001361",
            "SSB0001848",
            "SSB0002716",
            "SSB0002909",
            "SSB0002940",
            "SSB0002976",
        ],
        "version_tag": [1, 1, 1, 1, 1, 1, 1, 1, 1],
        "taxID": [None, 746128, 746128, None, 9606, None, None, None, None],
        "strain": ["A1163", None, None, None, "CEA10", None, None, None, None],
        "order": [
            "Eurotiales",
            "Eurotiales",
            "Eurotiales",
            None,
            "Eurotiales",
            None,
            None,
            None,
            None,
        ],
        "family": [
            "Aspergillaceae",
            "Aspergillaceae",
            "Aspergillaceae",
            None,
            "Aspergillaceae",
            None,
            None,
            None,
            None,
        ],
        "name": [
            "Aspergillus fumigatus",
            "Aspergillus fumigatus",
            "Aspergillus fumigatus",
            None,
            "Aspergillus fumigatus CEA10",
            None,
            None,
            None,
            None,
        ],
        "elementLength": [None, None, None, None, None, None, None, None, None],
        "upDR": [None, None, None, None, None, None, None, None, None],
        "downDR": [None, None, None, None, None, None, None, None, None],
        "contigID": [None, None, None, None, None, None, None, None, None],
        "captainID": [None, None, None, None, None, None, None, None, None],
        "elementBegin": [None, None, None, None, None, None, None, None, None],
        "elementEnd": [None, None, None, None, None, None, None, None, None],
        "familyName": [
            "Prometheus",
            "Prometheus",
            "Prometheus",
            None,
            None,
            None,
            None,
            None,
            None,
        ],
        "type_element_reference": [
            "Urquhart et al., unpublished",
            "Urquhart et al., unpublished",
            "Urquhart et al., unpublished",
            None,
            None,
            None,
            None,
            None,
            None,
        ],
        "navis_name": [None, "navis06", "navis06", None, None, None, None, None, None],
        "haplotype_name": [None, "var41", "var41", None, None, None, None, None, None],
        "ome": [
            "aspfum6",
            "CEA10",
            "CM2733",
            "aspfum137",
            "CEA10",
            None,
            None,
            None,
            None,
        ],
        "version": ["04-08-2014 0:00", None, None, None, None, None, None, None, None],
        "genomeSource": ["ncbi", None, None, None, None, None, None, None, None],
        "citation": [
            "26-06-2007",
            "10.1038/s41467-022-32924-7",
            "10.3390/genes9070363",
            "10.1038/s41564-021-00993-x",
            "10.1038/s41467-022-32924-7",
            None,
            None,
            None,
            None,
        ],
        "assembly_accession": [
            "Kaspfum6",
            "GCA_051225625.1",
            "10.5281/zenodo.5775265",
            "GCA_020501425.1",
            "GCA_051225625.1",
            None,
            None,
            None,
            None,
        ],
        "md5": [
            "946152598b782419263f931fdc45a122",
            "630b72bd313d36fabbb5a5dbcf415294",
            "c9db7e93d7312d58292748f190c0c755",
            "4dcdf281a059e4c5fff9303d0c729448",
            "4da37e0160b91c2ed64f9c1fb4e3abf1",
            "653c75e2e04a6082bba2b40097c85013",
            "548042213935d8b3e81911fc8a4c2588",
            "935fcdf406abeef437621e49e22a36b2",
            "b9c75132f0f6cb2747eae2152b203b80",
        ],
        "rev_comp_md5": [
            "5c923204abff36a5d07dde528ec09d3a",
            "c2e2176b5c36f1b85bf2de0060628e23",
            "4545f5f0908dc79ec22c4054572f04f1",
            "907a3b26480c1b8d958879a960ba8d47",
            "93584cb7b38ef0a3c182df19c676cf84",
            "df08149df7535e84e319a24dcdec85e3",
            "b127083b8eb439e67b95509d54af4bc0",
            "f2bcf5bfaf437ab6592c572c930512f0",
            "8d328163b67dd708dcd2b2ad6a54d98e",
        ],
        "sequence_length": [
            92517,
            93334,
            93249,
            92294,
            92365,
            92453,
            92518,
            92525,
            92525,
        ],
    }


@pytest.fixture
def test_haplotype_ships_df():
    return pd.DataFrame(
        {
            "accession_tag": ["SSA002851", "SSA002904", "SSA002596"],
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
            ],
        }
    )


@pytest.fixture
def test_captains_df():
    return pd.DataFrame(
        {
            "captainID": ["Phoenix", "Enterprise", "Voyager"],
            "sequence": [
                "MALWMRLLPLLALLALWGPDPAAA",  # Phoenix sequence
                "MALWMRLLPLLALLALWGPDPAAA",  # Same sequence for clustering
                "MALWMRLLPLLALLALWGPDPBBB",  # Different sequence
            ],
            "navis_name": ["Phoenix", "Phoenix", "Galactica"],
            "accession_tag": [
                "SSA002596",
                "SSA002904",
                "SSA002851",
            ],  # Multiple accession tags
        }
    )

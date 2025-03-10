import pytest
from dash.testing.application_runners import import_app
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium import webdriver

@pytest.fixture
def chrome_options():
    """Chrome options for the webdriver"""
    options = webdriver.ChromeOptions()
    options.add_argument('--headless')  # Run in headless mode (no GUI)
    options.add_argument('--no-sandbox')
    options.add_argument('--disable-dev-shm-usage')
    return options

@pytest.fixture
def driver(chrome_options):
    """Create a Chrome webdriver"""
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)
    yield driver
    driver.quit()

@pytest.fixture(scope="function")
def db_session():
    """Provide a database session for tests"""
    # Import the session class through the app to ensure proper initialization
    app = import_app("app")
    from src.config.database import StarbaseSession
    
    session = StarbaseSession
    yield session
    session.rollback()
    session.close()
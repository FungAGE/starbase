import pytest
from dash.testing.application_runners import import_app
from dash.testing.browser import Browser
import time
import multiprocessing
import requests
import os
import sys

# Add project root to path if needed
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if root_dir not in sys.path:
    sys.path.append(root_dir)


def heavy_request(url, duration=1):
    """Simulate a heavy request that takes resources"""
    start = time.time()
    while time.time() - start < duration:
        requests.get(url)
        time.sleep(0.1)

def test_server_under_load(dash_duo: Browser):
    """Test how the application behaves under heavy load"""
    # Import the app using dash's import_app utility
    app = import_app("app")
    
    # Start the dash server
    dash_duo.start_server(app)
    base_url = dash_duo.server_url
    
    # Create multiple processes to simulate load
    processes = []
    for _ in range(10):  # Create 10 concurrent processes
        p = multiprocessing.Process(
            target=heavy_request,
            args=(f"{base_url}/pgv",)
        )
        processes.append(p)
    
    try:
        # Start load generation
        for p in processes:
            p.start()
            
        # Try to load the page while under load
        dash_duo.wait_for_element("#pgv-table", timeout=10)
        
        # Check if we got any 500 errors
        logs = dash_duo.get_logs()
        errors = [log for log in logs if "500" in str(log)]
        
        # Assert we handled the errors gracefully
        for error in errors:
            assert "Internal Server Error" not in dash_duo.find_element("#pgv-table").text
            
    finally:
        # Clean up processes
        for p in processes:
            p.terminate()
            p.join()

def test_slow_database(dash_duo: Browser, monkeypatch):
    """Test how the application handles slow database responses"""
    # Import the app using dash's import_app utility
    app = import_app("app")
    
    # Patch the database function to be slow
    def slow_fetch(*args, **kwargs):
        time.sleep(5)  # Simulate slow DB
        return original_fetch(*args, **kwargs)
    
    # Get the original function before patching
    from src.database.sql_manager import fetch_ship_table as original_fetch
    monkeypatch.setattr("src.database.sql_manager.fetch_ship_table", slow_fetch)
    
    # Start the dash server
    dash_duo.start_server(app)
    
    # Try to load the page
    dash_duo.wait_for_element("#pgv-table", timeout=10)
    
    # Check if we handled the slow response gracefully
    assert "Error" not in dash_duo.find_element("#pgv-table").text
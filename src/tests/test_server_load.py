import pytest
from dash.testing.application_runners import import_app
from dash.testing.browser import Browser
import time
import multiprocessing
import requests
import os
import sys
from pathlib import Path
import logging
import base64  # Add this import for test_blast_submission

logger = logging.getLogger(__name__)

# Add project root to path if needed
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if root_dir not in sys.path:
    sys.path.append(root_dir)

sys.path.append(str(Path(__file__).parent.parent.parent))  # Add project root to Python path

def heavy_request(url, duration=1):
    """Simulate a heavy request that takes resources"""
    start = time.time()
    errors = []
    while time.time() - start < duration:
        try:
            response = requests.get(url, timeout=5)
            if response.status_code >= 500:
                errors.append(f"Got {response.status_code} error")
        except requests.exceptions.RequestException as e:
            errors.append(f"Request failed: {str(e)}")
        time.sleep(0.1)
    return errors

@pytest.mark.skip_browser
def test_server_under_load(dash_duo: Browser):
    """Test how the application behaves under heavy load"""
    # Import the app using dash's import_app utility
    app = import_app("app")
    
    # Start the server without threading parameter
    dash_duo.start_server(
        app,
        port=8050,
        host="127.0.0.1",
        debug=False,
        use_reloader=False,  # Keep this
        # Remove threaded=True as it's causing the conflict
    )
    
    base_url = dash_duo.server_url
    
    try:
        # Navigate to a specific page first
        dash_duo.driver.get(f"{base_url}/pgv")
        
        # Create multiple processes to simulate load
        processes = []
        manager = multiprocessing.Manager()
        error_list = manager.list()
        
        for _ in range(5):
            p = multiprocessing.Process(
                target=lambda: error_list.extend(heavy_request(f"{base_url}/pgv")),
            )
            processes.append(p)
        
        # Start load generation
        for p in processes:
            p.start()
            
        # Try to load the page while under load
        dash_duo.wait_for_element("#pgv-table", timeout=10)
        
        # Check for error handling
        error_elements = dash_duo.find_elements(".error-message")
        if error_elements:
            for elem in error_elements:
                assert "Internal Server Error" not in elem.text
                assert "try again" in elem.text.lower()
            
    finally:
        # Clean up processes
        for p in processes:
            if p.is_alive():
                p.terminate()
            p.join()
        
        # Make sure to quit the browser
        if hasattr(dash_duo, 'driver'):
            dash_duo.driver.quit()

@pytest.mark.skip_browser
def test_slow_database(dash_duo: Browser, monkeypatch):
    """Test how the application handles slow database responses"""
    app = import_app("app")
    
    def slow_fetch(*args, **kwargs):
        time.sleep(2)  # Reduce sleep time for testing
        return original_fetch(*args, **kwargs)
    
    from src.database.sql_manager import fetch_ship_table as original_fetch
    monkeypatch.setattr("src.database.sql_manager.fetch_ship_table", slow_fetch)
    
    dash_duo.start_server(
        app,
        port=8051,
        host="127.0.0.1",
        debug=False,
        use_reloader=False
    )
    
    # Try to load the page
    dash_duo.wait_for_element("#pgv-table", timeout=15)
    
    # Check for loading indicator
    loading = dash_duo.find_element(".loading-indicator")
    assert loading is not None

@pytest.mark.skip_browser
def test_timeout_handling(dash_duo: Browser):
    """Test that timeouts are handled gracefully"""
    app = import_app("app")
    
    dash_duo.start_server(
        app,
        port=8052,
        host="127.0.0.1",
        debug=False,
        use_reloader=False
    )
    
    # Navigate to blast page
    dash_duo.driver.get(f"{dash_duo.server_url}/blast")
    
    # Submit a large sequence
    large_sequence = "A" * 10000  # Reduced size for testing
    dash_duo.find_element("#query-text").send_keys(large_sequence)
    dash_duo.find_element("#submit-button").click()
    
    # Wait for error message
    error_elem = dash_duo.wait_for_element(".error-message", timeout=10)
    assert error_elem is not None
    assert "try again" in error_elem.text.lower()

def test_blast_submission(dash_duo: Browser):
    """Test submitting a BLAST search with a test file"""
    # Import the app
    app = import_app("app")
    
    # Start the server
    dash_duo.start_server(
        app,
        port=8052,  # Use a different port
        host="localhost"
    )
    
    # Read and encode the test file
    test_file_path = "tmp/aspfum3_s01168.fna"
    with open(test_file_path, "rb") as file:
        content = file.read()
        encoded = base64.b64encode(content).decode()
    
    # Create the file content string that Dash expects
    file_contents = f"data:application/octet-stream;base64,{encoded}"
    
    # Use JavaScript to set the file upload content
    dash_duo.driver.execute_script(
        """
        const fileUpload = document.getElementById('blast-fasta-upload');
        const event = new Event('drop');
        fileUpload.dispatchEvent(event);
        """
    )
    
    # Set the file contents using JavaScript
    dash_duo.driver.execute_script(
        f"""
        const contents = '{file_contents}';
        const filename = 'aspfum3_s01168.fna';
        window.dash_clientside.clientside.update_upload_status(contents, filename);
        """
    )
    
    # Click the submit button
    dash_duo.find_element("#submit-button").click()
    
    # Wait for results
    dash_duo.wait_for_element("#blast-table", timeout=30)
    
    # Verify results appeared
    results_table = dash_duo.find_element("#blast-table")
    assert results_table is not None
    
    # Check for specific content in results
    table_text = results_table.text
    assert "BLAST Results" in table_text
    
    # Optional: Check for specific matches if expected
    # assert "Expected Gene Name" in table_text
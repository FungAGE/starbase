import pytest
from dash.testing.application_runners import import_app
from dash.testing.browser import Browser
import time
import multiprocessing
import requests
import os
import sys
from pathlib import Path

# Add project root to path if needed
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if root_dir not in sys.path:
    sys.path.append(root_dir)

sys.path.append(str(Path(__file__).parent.parent.parent))  # Add project root to Python path

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
    
    # Start the server with specific host
    dash_duo.start_server(
        app,
        port=8050,
        host="localhost"  # Explicitly set host to localhost
    )
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
    
    # Start the server with specific host
    dash_duo.start_server(
        app,
        port=8051,  # Use a different port for second test
        host="localhost"  # Explicitly set host to localhost
    )
    
    # Try to load the page
    dash_duo.wait_for_element("#pgv-table", timeout=10)
    
    # Check if we handled the slow response gracefully
    assert "Error" not in dash_duo.find_element("#pgv-table").text
    assert "Error" not in dash_duo.find_element("#pgv-table").text

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
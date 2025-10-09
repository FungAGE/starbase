#!/usr/bin/env python3
"""
Test script for the database cleanup functionality.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.database.cleanup_accessions import main

if __name__ == "__main__":
    print("Testing database cleanup script...")
    print("Running in dry-run mode to analyze the database...")
    
    # Run the cleanup script in dry-run mode
    main(dry_run=True, output_report="cleanup_report.txt")
    
    print("Test completed. Check cleanup_report.txt for results.")
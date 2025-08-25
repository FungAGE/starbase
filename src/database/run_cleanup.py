#!/usr/bin/env python3
"""
Script to run optimized accession cleanup.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.database.cleanup_accessions import main

if __name__ == "__main__":
    print("Running optimized accession cleanup...")
    print("This will analyze sequences and identify duplicates/nested sequences with progress reporting.")
    print("Running in dry-run mode by default.")
    
    # Run the optimized accession cleanup
    main(dry_run=True, output_report="optimized_cleanup_report.txt")
    
    print("Optimized accession cleanup completed. Check optimized_cleanup_report.txt for results.")
    print("To apply changes, run: python run_optimized_cleanup.py --apply")
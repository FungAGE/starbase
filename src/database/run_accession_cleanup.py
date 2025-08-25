#!/usr/bin/env python3
"""
Script to run complete accession cleanup including reverse complement detection.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.database.cleanup_accessions import main

if __name__ == "__main__":
    print("Running complete accession cleanup...")
    print("This will analyze sequences and identify:")
    print("1. Self-complementary sequences (normal = reverse complement)")
    print("2. Reverse complement pairs")
    print("3. Nested sequences")
    print("4. Exact duplicates")
    print("\nRunning in dry-run mode by default.")
    
    # Run the complete accession cleanup
    main(dry_run=True, output_report="complete_cleanup_report.txt")
    
    print("\nComplete accession cleanup finished. Check complete_cleanup_report.txt for results.")
    print("To apply changes, run: python run_accession_cleanup.py --apply")

#!/usr/bin/env python3
"""
Runner script for comprehensive database cleanup.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from src.database.database_cleanup import main

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run comprehensive database cleanup")
    parser.add_argument("--apply", action="store_true", 
                       help="Apply changes to database (default is dry run)")
    parser.add_argument("--report", type=str, 
                       help="Path to save cleanup report")
    
    args = parser.parse_args()
    
    main(dry_run=not args.apply, output_report=args.report)

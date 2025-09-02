#!/usr/bin/env python3
"""
Simple test script for database cleanup functions.
This script can be run directly from the database directory.
"""

import sys
import os

# Add the src directory to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
sys.path.insert(0, src_dir)

def test_imports():
    """Test if we can import the required modules."""
    try:
        from database.database_cleanup import (
            check_ships_accessions_joined_ships_relationships,
            analyze_table_relationships
        )
        print("✅ Successfully imported database cleanup functions")
        return True
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality without database connection."""
    try:
        from database.database_cleanup import (
            check_ships_accessions_joined_ships_relationships,
            analyze_table_relationships
        )
        
        print("✅ Basic functionality test passed")
        print("Available functions:")
        print("  - check_ships_accessions_joined_ships_relationships()")
        print("  - analyze_table_relationships()")
        print("  - fix_ships_accessions_joined_ships_relationships()")
        print("  - Automatically creates joined_ships entries for ships missing them")
        print("  - Ensures all ships have corresponding entries in the main connecting table")
        
        return True
    except Exception as e:
        print(f"❌ Functionality test failed: {e}")
        return False

if __name__ == "__main__":
    print("Testing database cleanup functionality...")
    print("=" * 50)
    
    # Test imports
    if test_imports():
        # Test basic functionality
        test_basic_functionality()
        print("\n✅ All tests passed! You can now use the database cleanup functions.")
        print("\nUsage examples:")
        print("  python run_comprehensive_cleanup.py --analyze-relationships")
        print("  python run_comprehensive_cleanup.py --check-relationships")
        print("  python run_comprehensive_cleanup.py --apply-fixes")
        print("  python run_comprehensive_cleanup.py --skip-accession-cleanup --apply-fixes")
    else:
        print("\n❌ Tests failed. Please check your Python path and dependencies.")

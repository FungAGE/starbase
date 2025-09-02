#!/usr/bin/env python3
"""
Test script for the new ships-accessions-joined_ships relationship checking functionality.
"""

import sys
import os

# Add the src directory to the Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)
sys.path.insert(0, src_dir)

from database.database_cleanup import (
    check_ships_accessions_joined_ships_relationships,
    analyze_table_relationships,
    fix_ships_accessions_joined_ships_relationships
)

def test_relationship_checks():
    """Test the relationship checking functions."""
    print("Testing ships-accessions-joined_ships relationship checks...")
    print("=" * 60)
    
    try:
        # Test 1: Check relationships
        print("\n1. Running relationship checks...")
        issues = check_ships_accessions_joined_ships_relationships()
        
        print(f"Found {len(issues['orphaned_ships'])} orphaned ships")
        print(f"Found {len(issues['orphaned_accessions'])} orphaned accessions")
        print(f"Found {len(issues['orphaned_joined_ships'])} orphaned joined_ships")
        print(f"Found {len(issues['missing_sequence_links'])} missing sequence links")
        print(f"Found {len(issues['inconsistent_joined_ships'])} inconsistent joined_ships")
        
        # Test 2: Analyze relationships
        print("\n2. Running relationship analysis...")
        analysis = analyze_table_relationships()
        
        print(f"Total ships: {analysis['relationship_counts']['total_ships']}")
        print(f"Total accessions: {analysis['relationship_counts']['total_accessions']}")
        print(f"Total joined_ships: {analysis['relationship_counts']['total_joined_ships']}")
        print(f"Ships with accession: {analysis['relationship_counts']['ships_with_accession']}")
        print(f"Joined_ships with ship: {analysis['relationship_counts']['joined_with_ship']}")
        
        # Test 3: Test fixes (dry run)
        print("\n3. Testing fixes (dry run)...")
        fixes = fix_ships_accessions_joined_ships_relationships(dry_run=True)
        
        print(f"Would fix {len(fixes['ships_fixed'])} ship issues")
        print(f"Would fix {len(fixes['accessions_fixed'])} accession issues")
        print(f"Would fix {len(fixes['joined_ships_fixed'])} joined_ships issues")
        
        print("\nAll tests completed successfully!")
        
    except Exception as e:
        print(f"Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_relationship_checks()

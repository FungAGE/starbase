#!/usr/bin/env python3
"""
Test script to verify nested sequence detection logic.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.database.cleanup_accessions import find_nested_sequences, fetch_all_sequences
import pandas as pd

def test_nested_detection():
    """Test the nested sequence detection with a small sample."""
    
    print("Testing nested sequence detection logic...")
    
    # Create a small test dataset
    test_data = {
        'accession_id': [1, 2, 3, 4, 5],
        'accession_tag': ['SSA000001', 'SSA000002', 'SSA000003', 'SSA000004', 'SSA000005'],
        'version_tag': [None, None, None, None, None],
        'sequence': [
            'ATCGATCGATCG',  # Longest sequence
            'ATCGATCG',      # Nested in sequence 1
            'GCTAGCTA',      # Different sequence
            'ATCG',          # Nested in sequences 1 and 2
            'ATCGATCGATCG'   # Identical to sequence 1 (should be filtered out)
        ],
        'md5': ['hash1', 'hash2', 'hash3', 'hash4', 'hash5']
    }
    
    test_df = pd.DataFrame(test_data)
    
    print("Test sequences:")
    for _, row in test_df.iterrows():
        print(f"  {row['accession_tag']}: {row['sequence']} (length: {len(row['sequence'])})")
    
    # Test the nested detection
    nested_pairs = find_nested_sequences(test_df)
    
    print(f"\nFound {len(nested_pairs)} nested sequence pairs:")
    for nested, containing in nested_pairs:
        print(f"  {nested} is nested in {containing}")
    
    # Verify the results
    expected_pairs = [
        ('SSA000004', 'SSA000001'),  # ATCG nested in ATCGATCGATCG
        ('SSA000004', 'SSA000002'),  # ATCG nested in ATCGATCG
        ('SSA000002', 'SSA000001'),  # ATCGATCG nested in ATCGATCGATCG
    ]
    
    print(f"\nExpected pairs: {len(expected_pairs)}")
    for nested, containing in expected_pairs:
        print(f"  {nested} is nested in {containing}")
    
    # Check if results match expectations
    actual_pairs = set(nested_pairs)
    expected_pairs_set = set(expected_pairs)
    
    if actual_pairs == expected_pairs_set:
        print("\n✅ Test PASSED: All expected nested sequences found correctly!")
    else:
        print("\n❌ Test FAILED:")
        print(f"  Expected: {expected_pairs_set}")
        print(f"  Actual: {actual_pairs}")
        print(f"  Missing: {expected_pairs_set - actual_pairs}")
        print(f"  Extra: {actual_pairs - expected_pairs_set}")

if __name__ == "__main__":
    test_nested_detection()

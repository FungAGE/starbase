#!/usr/bin/env python3
"""
Test script to verify reverse complement detection logic.
"""

import sys
import os

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from src.database.cleanup_accessions import find_reverse_complement_pairs, generate_sequence_hashes
import pandas as pd

def test_reverse_complement_detection():
    """Test the reverse complement detection with a small sample."""
    
    print("Testing reverse complement detection logic...")
    
    # Create a small test dataset
    test_data = {
        'accession_id': [1, 2, 3, 4, 5, 6],
        'accession_tag': ['SSA000001', 'SSA000002', 'SSA000003', 'SSA000004', 'SSA000005', 'SSA000006'],
        'version_tag': [None, None, None, None, None, None],
        'sequence': [
            'ATCGATCG',      # Normal sequence
            'CGATCGAT',      # Reverse complement of SSA000001
            'GCTAGCTA',      # Different sequence
            'TAGCTAGC',      # Reverse complement of SSA000003
            'ATCGATCG',      # Identical to SSA000001 (should be filtered out as duplicate)
            'TTTTTTTT',      # Self-complementary sequence
        ],
        'md5': ['hash1', 'hash2', 'hash3', 'hash4', 'hash5', 'hash6']
    }
    
    test_df = pd.DataFrame(test_data)
    
    print("Test sequences:")
    for _, row in test_df.iterrows():
        normal_hash, rev_comp_hash = generate_sequence_hashes(row['sequence'])
        print(f"  {row['accession_tag']}: {row['sequence']}")
        print(f"    Normal hash: {normal_hash}")
        print(f"    Rev comp hash: {rev_comp_hash}")
        print(f"    Self-complementary: {normal_hash == rev_comp_hash}")
    
    # Test the reverse complement detection
    rev_comp_pairs = find_reverse_complement_pairs(test_df)
    
    print(f"\nFound {len(rev_comp_pairs)} reverse complement pairs:")
    for acc1, acc2 in rev_comp_pairs:
        print(f"  {acc1} <-> {acc2}")
    
    # Verify the results
    expected_pairs = [
        ('SSA000001', 'SSA000002'),  # ATCGATCG <-> CGATCGAT
        ('SSA000003', 'SSA000004'),  # GCTAGCTA <-> TAGCTAGC
        ('SSA000005', 'SSA000002'),  # SSA000005 is identical to SSA000001, so also <-> SSA000002
    ]
    
    print(f"\nExpected pairs: {len(expected_pairs)}")
    for acc1, acc2 in expected_pairs:
        print(f"  {acc1} <-> {acc2}")
    
    # Check if results match expectations
    actual_pairs = set()
    for acc1, acc2 in rev_comp_pairs:
        actual_pairs.add(tuple(sorted([acc1, acc2])))
    
    expected_pairs_set = set()
    for acc1, acc2 in expected_pairs:
        expected_pairs_set.add(tuple(sorted([acc1, acc2])))
    
    if actual_pairs == expected_pairs_set:
        print("\n✅ Test PASSED: All expected reverse complement pairs found correctly!")
    else:
        print("\n❌ Test FAILED:")
        print(f"  Expected: {expected_pairs_set}")
        print(f"  Actual: {actual_pairs}")
        print(f"  Missing: {expected_pairs_set - actual_pairs}")
        print(f"  Extra: {actual_pairs - expected_pairs_set}")

if __name__ == "__main__":
    test_reverse_complement_detection()

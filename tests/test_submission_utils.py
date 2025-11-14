#!/usr/bin/env python3
"""
Test script for submission_utils.

This demonstrates various use cases for the submission utility.
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.utils.submission_utils import (
    SubmissionProcessor,
    process_single,
    validate_submission,
    check_sequence_duplicate
)


def test_validation():
    """Test input validation."""
    print("=" * 80)
    print("TEST 1: Input Validation")
    print("=" * 80)
    
    # Valid submission
    valid_data = {
        'sequence': 'ATGCATGCATGCATGC' * 100,  # Valid sequence
        'starshipID': 'SS-test-1',
        'evidence': 'starfish',
        'source': 'test_data',
        'genus': 'Fusarium',
        'species': 'oxysporum'
    }
    
    is_valid, errors = validate_submission(valid_data)
    print(f"Valid submission: {is_valid}")
    if errors:
        print(f"Errors: {errors}")
    
    # Invalid submission - missing required field
    invalid_data = {
        'sequence': 'ATGC',
        'starshipID': 'SS-test-2',
        # Missing 'evidence' and 'source'
    }
    
    is_valid, errors = validate_submission(invalid_data)
    print(f"\nInvalid submission: {is_valid}")
    print(f"Errors: {errors}")
    
    print()


def test_duplicate_detection():
    """Test duplicate sequence detection."""
    print("=" * 80)
    print("TEST 2: Duplicate Detection")
    print("=" * 80)
    
    # Test with a new sequence
    test_sequence = 'ATGCATGCATGC' * 200
    
    duplicate_info = check_sequence_duplicate(
        test_sequence,
        genus='Fusarium',
        species='oxysporum'
    )
    
    print(f"Is duplicate: {duplicate_info.is_duplicate}")
    if duplicate_info.is_duplicate:
        print(f"Existing ship_id: {duplicate_info.existing_ship_id}")
        print(f"Existing accession: {duplicate_info.existing_accession}")
        print(f"Different taxon: {duplicate_info.different_taxon}")
        print(f"Match type: {duplicate_info.match_type}")
    
    print()


def test_single_submission():
    """Test processing a single submission."""
    print("=" * 80)
    print("TEST 3: Single Submission (Dry Run)")
    print("=" * 80)
    
    processor = SubmissionProcessor(dry_run=True)
    
    submission = {
        'sequence': 'ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT' * 100,
        'starshipID': 'SS-test-fusarium-1',
        'evidence': 'manual_curation',
        'source': 'test_dataset_2024',
        'genus': 'Fusarium',
        'species': 'oxysporum',
        'strain': 'f.sp. lycopersici',
        'ship_family': 'Voyager',
        'assembly_accession': 'GCA_000000001.1',
        'curated_status': 'curated',
        'notes': 'Test submission for validation'
    }
    
    result = processor.process_submission(submission)
    
    print(f"Success: {result['success']}")
    print(f"Ship ID: {result['ship_id']}")
    print(f"Accession: {result['accession']}")
    if result['errors']:
        print(f"Errors: {result['errors']}")
    if result['warnings']:
        print(f"Warnings: {result['warnings']}")
    
    print()


def test_batch_submission():
    """Test processing multiple submissions."""
    print("=" * 80)
    print("TEST 4: Batch Submission (Dry Run)")
    print("=" * 80)
    
    submissions = [
        {
            'sequence': 'ATGCATGCATGC' * 150,
            'starshipID': 'SS-batch-1',
            'evidence': 'starfish',
            'source': 'batch_test',
            'genus': 'Aspergillus',
            'species': 'fumigatus',
        },
        {
            'sequence': 'GCTAGCTAGCTA' * 150,
            'starshipID': 'SS-batch-2',
            'evidence': 'starfish',
            'source': 'batch_test',
            'genus': 'Aspergillus',
            'species': 'nidulans',
        },
        {
            # Invalid - missing required field
            'sequence': 'ATGC' * 100,
            'starshipID': 'SS-batch-3',
            'evidence': 'starfish',
            # Missing 'source'
        }
    ]
    
    processor = SubmissionProcessor(dry_run=True)
    results = processor.process_batch(submissions)
    
    print(f"\nProcessed {len(results)} submissions")
    for i, result in enumerate(results, 1):
        print(f"\nSubmission {i}:")
        print(f"  Success: {result['success']}")
        if result['errors']:
            print(f"  Errors: {result['errors']}")
    
    print()


def test_convenience_function():
    """Test the convenience function."""
    print("=" * 80)
    print("TEST 5: Convenience Function")
    print("=" * 80)
    
    # Note: This is just for demonstration - shows the simple API
    print("Example usage of process_single():")
    print("""
    result = process_single(
        sequence='ATGCATGC...',
        starshipID='SS-1.1',
        evidence='starfish',
        source='publication_2024',
        genus='Fusarium',
        species='oxysporum',
        ship_family='Voyager'
    )
    """)
    print()


if __name__ == '__main__':
    print("\n" + "=" * 80)
    print("SUBMISSION UTILS TEST SUITE")
    print("=" * 80 + "\n")
    
    try:
        test_validation()
        test_duplicate_detection()
        test_single_submission()
        test_batch_submission()
        test_convenience_function()
        
        print("=" * 80)
        print("ALL TESTS COMPLETED")
        print("=" * 80)
        
    except Exception as e:
        print(f"\n❌ Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()


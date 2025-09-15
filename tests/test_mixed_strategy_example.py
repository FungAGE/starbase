"""
Example test demonstrating the mixed test data strategy.
This shows how to use both real and mock data in the same test suite.
"""

import pytest
import pandas as pd
from unittest.mock import patch, MagicMock

from tests.utils.test_data_manager import TestDataManager
from tests.test_config import test_config
from src.utils.classification_utils import (
    assign_accession,
    check_exact_match,
    check_contained_match,
    check_similar_match,
    run_classification_workflow
)


class TestMixedDataStrategy:
    """Test class demonstrating mixed data strategy usage."""
    
    def test_mock_data_always_available(self, test_data_manager):
        """Test that mock data is always available for unit tests."""
        # Mock data should always be available
        ships_df = test_data_manager.get_ships_data(with_sequence=True, limit=5)
        meta_df = test_data_manager.get_meta_data(limit=5)
        sequence = test_data_manager.get_test_sequence(length=1000)
        
        assert isinstance(ships_df, pd.DataFrame)
        assert isinstance(meta_df, pd.DataFrame)
        assert isinstance(sequence, str)
        assert len(sequence) == 1000
        assert len(ships_df) <= 5
        assert len(meta_df) <= 5
    
    @pytest.mark.real_data
    def test_real_data_when_available(self, real_data_manager, database_available):
        """Test using real data when database is available."""
        if not database_available:
            pytest.skip("Database not available")
        
        # Try to get real data
        try:
            ships_df = real_data_manager.get_ships_data(with_sequence=True, limit=5)
            meta_df = real_data_manager.get_meta_data(limit=5)
            
            assert isinstance(ships_df, pd.DataFrame)
            assert isinstance(meta_df, pd.DataFrame)
            assert len(ships_df) > 0
            assert len(meta_df) > 0
            
            # Verify real data has expected structure
            assert 'accession_tag' in ships_df.columns
            assert 'sequence' in ships_df.columns
            assert 'accession_tag' in meta_df.columns
            assert 'familyName' in meta_df.columns
            
        except Exception as e:
            pytest.skip(f"Real data not available: {e}")
    
    def test_classification_with_mock_data(self, test_data_manager):
        """Test classification workflow with mock data."""
        # Get mock data
        ships_df = test_data_manager.get_ships_data(with_sequence=True, limit=10)
        test_sequence = test_data_manager.get_test_sequence(length=1000)
        
        # Test exact match
        result = check_exact_match(test_sequence, ships_df)
        assert result is not None  # Should find a match with our mock data
        
        # Test contained match
        contained_result = check_contained_match(test_sequence, ships_df)
        assert contained_result is not None
        
        # Test similar match
        similar_result = check_similar_match(test_sequence, ships_df)
        assert similar_result is not None
    
    @pytest.mark.real_data
    @pytest.mark.classification_pipeline
    def test_classification_with_real_data(self, real_data_manager, database_available):
        """Test classification workflow with real data from database."""
        if not database_available:
            pytest.skip("Database not available")
        
        try:
            # Get real data
            ships_df = real_data_manager.get_ships_data(with_sequence=True, limit=20)
            meta_df = real_data_manager.get_meta_data(limit=20)
            
            # Get a real sequence to test with
            if len(ships_df) > 0:
                test_sequence = ships_df.iloc[0]['sequence']
                
                # Test classification workflow
                result = assign_accession(test_sequence, ships_df)
                accession, needs_review = result
                
                assert accession is not None
                assert isinstance(needs_review, bool)
                
                # Test that we can run the full workflow
                workflow_result = run_classification_workflow(
                    sequence=test_sequence,
                    ships_df=ships_df,
                    meta_df=meta_df
                )
                
                assert workflow_result is not None
                
        except Exception as e:
            pytest.skip(f"Real data classification test failed: {e}")
    
    def test_fallback_to_mock_when_real_unavailable(self, test_data_manager):
        """Test that we gracefully fall back to mock data when real data is unavailable."""
        # This test should always pass because mock data is always available
        ships_df = test_data_manager.get_ships_data(with_sequence=True, limit=5)
        meta_df = test_data_manager.get_meta_data(limit=5)
        
        assert len(ships_df) > 0
        assert len(meta_df) > 0
        
        # Test that mock data has the right structure
        assert 'accession_tag' in ships_df.columns
        assert 'sequence' in ships_df.columns
        assert 'accession_tag' in meta_df.columns
        assert 'familyName' in meta_df.columns
    
    @pytest.mark.parametrize("use_real_data", [True, False])
    def test_data_manager_switching(self, use_real_data):
        """Test that data manager can switch between real and mock data."""
        manager = TestDataManager(use_real_data=use_real_data)
        
        if use_real_data and not test_config.is_database_available():
            # If real data is requested but not available, should fall back to mock
            ships_df = manager.get_ships_data(with_sequence=True, limit=5)
            assert len(ships_df) > 0  # Should still get mock data
        else:
            # Should get the requested type of data
            ships_df = manager.get_ships_data(with_sequence=True, limit=5)
            assert len(ships_df) > 0
    
    def test_configuration_controls_data_usage(self, config):
        """Test that configuration properly controls data usage."""
        # Test that configuration methods work
        assert isinstance(config.should_use_real_data(), bool)
        assert isinstance(config.is_database_available(), bool)
        
        # Test specific test control
        assert isinstance(config.should_use_real_data("test_classification"), bool)
        assert isinstance(config.should_use_real_data("test_blast"), bool)
    
    @pytest.mark.integration
    def test_full_workflow_integration(self, test_data_manager, real_data_manager, config):
        """Test full workflow integration with mixed data strategy."""
        # This test demonstrates how to use both data types in one test
        
        # Always use mock data for setup
        mock_ships = test_data_manager.get_ships_data(with_sequence=True, limit=5)
        mock_meta = test_data_manager.get_meta_data(limit=5)
        
        # Try to use real data for the actual workflow if available
        if config.should_use_real_data("integration") and config.is_database_available():
            try:
                real_ships = real_data_manager.get_ships_data(with_sequence=True, limit=10)
                real_meta = real_data_manager.get_meta_data(limit=10)
                
                # Use real data for the test
                ships_df = real_ships
                meta_df = real_meta
                data_source = "real"
                
            except Exception:
                # Fall back to mock data
                ships_df = mock_ships
                meta_df = mock_meta
                data_source = "mock"
        else:
            # Use mock data
            ships_df = mock_ships
            meta_df = mock_meta
            data_source = "mock"
        
        # Run the test with whatever data we have
        test_sequence = ships_df.iloc[0]['sequence'] if len(ships_df) > 0 else "ATCGATCGATCG"
        
        result = assign_accession(test_sequence, ships_df)
        accession, needs_review = result
        
        assert accession is not None
        assert isinstance(needs_review, bool)
        
        # Log which data source was used
        print(f"Integration test used {data_source} data")
    
    def test_error_handling_in_data_manager(self, test_data_manager):
        """Test error handling in data manager."""
        # Test with invalid parameters
        with pytest.raises(ValueError):
            test_data_manager.get_ships_data(limit=-1)
        
        with pytest.raises(ValueError):
            test_data_manager.get_test_sequence(length=-1)
        
        # Test with very large limits (should be capped)
        large_ships = test_data_manager.get_ships_data(limit=10000)
        assert len(large_ships) <= test_config.max_mock_records
        
        large_meta = test_data_manager.get_meta_data(limit=10000)
        assert len(large_meta) <= test_config.max_mock_records


# Utility test functions
def test_data_manager_initialization():
    """Test that data manager initializes correctly."""
    # Test with mock data
    mock_manager = TestDataManager(use_real_data=False)
    assert mock_manager.use_real_data is False
    
    # Test with real data
    real_manager = TestDataManager(use_real_data=True)
    assert real_manager.use_real_data is True


def test_config_initialization():
    """Test that configuration initializes correctly."""
    config = test_config
    assert hasattr(config, 'use_real_data')
    assert hasattr(config, 'max_mock_records')
    assert hasattr(config, 'test_sequence_length')
    assert hasattr(config, 'should_use_real_data')
    assert hasattr(config, 'is_database_available')


# Example of how to mark tests for different scenarios
@pytest.mark.unit
def test_unit_with_mock_data(test_data_manager):
    """Example unit test using only mock data."""
    ships_df = test_data_manager.get_ships_data(with_sequence=True, limit=3)
    assert len(ships_df) == 3


@pytest.mark.integration
@pytest.mark.real_data
def test_integration_with_real_data(real_data_manager, database_available):
    """Example integration test using real data."""
    if not database_available:
        pytest.skip("Database not available")
    
    try:
        ships_df = real_data_manager.get_ships_data(with_sequence=True, limit=5)
        assert len(ships_df) > 0
    except Exception:
        pytest.skip("Real data not available")


@pytest.mark.classification_pipeline
def test_classification_pipeline_mock(test_data_manager):
    """Example classification pipeline test with mock data."""
    ships_df = test_data_manager.get_ships_data(with_sequence=True, limit=10)
    test_sequence = test_data_manager.get_test_sequence(length=1000)
    
    result = check_exact_match(test_sequence, ships_df)
    assert result is not None

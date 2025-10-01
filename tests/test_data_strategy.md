# Test Data Strategy

This document outlines the mixed test data strategy implemented for the Starbase project, which combines real data from the database with mock data for comprehensive testing.

## Overview

The test data strategy addresses the challenge of testing with large sequence files while maintaining efficient Git workflows. Instead of tracking large test files in Git, we use a mixed approach:

- **Real Data**: Pulled from the database using existing functions (`fetch_meta_data`, `fetch_ships`) for tests involving the BLAST/classification pipeline
- **Mock Data**: Generated programmatically for unit tests and other testing scenarios

## Benefits

1. **No Large Files in Git**: Large sequence files are not tracked, keeping the repository lightweight
2. **Real Data Testing**: Critical pipeline tests use actual data from the database
3. **Fast Unit Tests**: Mock data enables fast, isolated unit tests
4. **Flexible Configuration**: Easy to switch between real and mock data based on test requirements
5. **Database Independence**: Tests can run without database access using mock data

## Architecture

### Core Components

1. **TestDataManager** (`tests/utils/test_data_manager.py`)
   - Central class for managing test data
   - Handles both real and mock data generation
   - Provides consistent interface for all test data needs

2. **Test Configuration** (`tests/test_config.py`)
   - Controls when to use real vs mock data
   - Configures data limits and parameters
   - Manages database availability checks

3. **Updated Conftest** (`tests/conftest_updated.py`)
   - Provides fixtures for both real and mock data
   - Integrates with the new test data strategy
   - Maintains backward compatibility

### Data Sources

#### Real Data
- **Source**: Database via `fetch_meta_data()` and `fetch_ships()`
- **Use Cases**: BLAST pipeline, classification workflow, integration tests
- **Requirements**: Database connection, network access
- **Fallback**: Automatically falls back to mock data if unavailable

#### Mock Data
- **Source**: Programmatically generated test data
- **Use Cases**: Unit tests, isolated component testing, CI/CD pipelines
- **Requirements**: None (always available)
- **Characteristics**: Deterministic, fast, lightweight

## Usage

### Basic Usage

```python
import pytest
from tests.utils.test_data_manager import TestDataManager

def test_with_mock_data():
    """Test using mock data (always available)."""
    manager = TestDataManager(use_real_data=False)
    ships_df = manager.get_ships_data(with_sequence=True, limit=5)
    meta_df = manager.get_meta_data(limit=5)
    
    assert len(ships_df) > 0
    assert len(meta_df) > 0

@pytest.mark.real_data
def test_with_real_data():
    """Test using real data from database."""
    manager = TestDataManager(use_real_data=True)
    try:
        ships_df = manager.get_ships_data(with_sequence=True, limit=10)
        meta_df = manager.get_meta_data(limit=10)
        
        # Test with real data
        assert len(ships_df) > 0
        assert len(meta_df) > 0
    except Exception:
        pytest.skip("Real data not available")
```

### Using Fixtures

```python
def test_with_fixtures(mock_ships_df, mock_meta_df):
    """Test using pytest fixtures."""
    assert len(mock_ships_df) > 0
    assert len(mock_meta_df) > 0

@pytest.mark.real_data
def test_with_real_fixtures(real_ships_df, real_meta_df):
    """Test using real data fixtures."""
    assert len(real_ships_df) > 0
    assert len(real_meta_df) > 0
```

### Configuration-Based Testing

```python
def test_with_config(config, test_data_manager, real_data_manager):
    """Test that adapts based on configuration."""
    if config.should_use_real_data("classification"):
        manager = real_data_manager
    else:
        manager = test_data_manager
    
    ships_df = manager.get_ships_data(with_sequence=True, limit=5)
    assert len(ships_df) > 0
```

## Test Markers

The strategy uses pytest markers to categorize tests:

- `@pytest.mark.real_data`: Test requires real data from database
- `@pytest.mark.mock_data`: Test uses only mock data
- `@pytest.mark.blast_pipeline`: Test is part of BLAST pipeline
- `@pytest.mark.classification_pipeline`: Test is part of classification pipeline
- `@pytest.mark.integration`: Integration test
- `@pytest.mark.unit`: Unit test

## Configuration

### Environment Variables

- `USE_REAL_DATA`: Set to "true" to enable real data usage
- `DATABASE_URL`: Database connection string
- `MAX_MOCK_RECORDS`: Maximum number of mock records to generate
- `TEST_SEQUENCE_LENGTH`: Length of generated test sequences

### Test-Specific Configuration

```python
# In test_config.py
REAL_DATA_TESTS = {
    "test_classification": True,
    "test_blast": True,
    "test_integration": True,
    "test_unit": False
}
```

## File Organization

```
tests/
├── conftest_updated.py          # Updated pytest configuration
├── test_config.py               # Test configuration
├── test_mixed_strategy_example.py  # Example tests
├── utils/
│   └── test_data_manager.py     # Core data management
├── test_data/
│   ├── large_sequences/         # Large files (not tracked)
│   │   ├── multi.fasta
│   │   └── SBS000357.1.fasta
│   └── mock/                    # Small mock files (tracked)
│       └── small_test_sequence.fasta
└── test_classification_mixed.py # Updated classification tests
```

## Migration Guide

### From Old Test Files

1. **Replace direct file access**:
   ```python
   # Old way
   with open("tests/test_data/large_file.fasta") as f:
       content = f.read()
   
   # New way
   manager = TestDataManager(use_real_data=True)
   sequence = manager.get_test_sequence(length=1000)
   ```

2. **Update fixtures**:
   ```python
   # Old way
   @pytest.fixture
   def mock_ships_df():
       return pd.DataFrame({...})
   
   # New way
   def test_with_fixture(mock_ships_df):
       # mock_ships_df is provided by conftest_updated.py
   ```

3. **Add appropriate markers**:
   ```python
   @pytest.mark.real_data
   @pytest.mark.classification_pipeline
   def test_classification_workflow():
       # Test implementation
   ```

### Database Integration

The strategy integrates with existing database functions:

```python
from src.database.sql_manager import fetch_meta_data, fetch_ships

# These functions are used internally by TestDataManager
# when use_real_data=True
```

## Best Practices

### When to Use Real Data

- **BLAST Pipeline Tests**: Require real sequence data for accurate testing
- **Classification Workflow**: Need real data to test classification accuracy
- **Integration Tests**: Should use real data to test end-to-end workflows
- **Performance Tests**: Real data provides realistic performance metrics

### When to Use Mock Data

- **Unit Tests**: Fast, isolated testing of individual components
- **CI/CD Pipelines**: No database dependency
- **Development**: Quick feedback during development
- **Edge Cases**: Easy to generate specific test scenarios

### Error Handling

Always handle cases where real data is unavailable:

```python
@pytest.mark.real_data
def test_with_real_data(real_data_manager, database_available):
    if not database_available:
        pytest.skip("Database not available")
    
    try:
        data = real_data_manager.get_ships_data(limit=5)
        # Test with real data
    except Exception as e:
        pytest.skip(f"Real data not available: {e}")
```

## Troubleshooting

### Common Issues

1. **Database Connection Errors**
   - Check `DATABASE_URL` environment variable
   - Verify database is running and accessible
   - Use mock data as fallback

2. **Large Data Performance**
   - Use `limit` parameter to control data size
   - Consider using mock data for development
   - Use real data only when necessary

3. **Test Failures with Real Data**
   - Check if database schema has changed
   - Verify data availability in database
   - Add appropriate error handling

### Debugging

Enable debug logging to see which data source is being used:

```python
import logging
logging.basicConfig(level=logging.DEBUG)

# TestDataManager will log which data source it's using
```

## Future Enhancements

1. **Data Caching**: Cache real data to reduce database calls
2. **Data Validation**: Validate real data structure and content
3. **Performance Metrics**: Track test performance with different data sources
4. **Data Versioning**: Handle database schema changes gracefully
5. **Automated Data Refresh**: Periodically update test data from database

## Conclusion

The mixed test data strategy provides a robust, flexible approach to testing that balances the need for real data accuracy with the practical requirements of efficient development and CI/CD workflows. By using real data for critical pipeline tests and mock data for unit tests, we achieve comprehensive test coverage while maintaining fast, reliable test execution.
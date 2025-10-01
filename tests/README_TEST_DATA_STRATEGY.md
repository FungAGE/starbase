# Test Data Strategy - Mixed Approach

This document describes the mixed strategy for managing test data in the Starbase project.

## Overview

The mixed strategy addresses the challenge of testing with large sequence files while maintaining efficient git workflows:

1. **Real Data from Database**: Use `fetch_meta_data()` and `fetch_ships()` for BLAST/classification pipeline tests
2. **Mock Data**: Use generated mock data for unit tests and other functionality
3. **No Large Files in Git**: Large sequence files are not tracked by git

## Key Components

### 1. TestDataManager (`tests/utils/test_data_manager.py`)

The central class for managing test data:

```python
from tests.utils.test_data_manager import TestDataManager

# For mock data (default)
manager = TestDataManager(use_real_data=False)
ships_df = manager.get_ships_data(with_sequence=True)

# For real data from database
manager = TestDataManager(use_real_data=True)
ships_df = manager.get_ships_data(with_sequence=True, limit=10)
```

### 2. Pytest Fixtures

Ready-to-use fixtures for common test scenarios:

```python
def test_with_mock_data(mock_ships_df, mock_meta_df):
    # Uses mock data automatically
    pass

def test_with_real_data(real_ships_df, real_meta_df):
    # Uses real data from database
    pass
```

### 3. Test Configuration (`tests/test_config.py`)

Environment-based configuration:

```bash
# Enable real data for tests
export STARBASE_TEST_USE_REAL_DATA=true

# Set maximum records for real data
export STARBASE_TEST_MAX_REAL_DATA_LIMIT=50

# Set custom database URL
export STARBASE_TEST_DATABASE_URL=sqlite:///test.db
```

## Usage Patterns

### Unit Tests (Mock Data)

For testing individual functions in isolation:

```python
def test_assign_accession_new_sequence(test_data_manager):
    """Test accession assignment with mock data."""
    ships_df = test_data_manager.get_ships_data(with_sequence=True)
    new_sequence = "ATCGATCGATCGATCGATCG"
    
    accession, needs_review = assign_accession(new_sequence, ships_df)
    
    assert accession.startswith("SBS")
    assert needs_review is False
```

### Integration Tests (Real Data)

For testing BLAST/classification pipeline with real sequences:

```python
@use_real_data
def test_classification_workflow_with_real_data(real_data_manager):
    """Test full classification workflow with real data."""
    ships_df = real_data_manager.get_ships_data(with_sequence=True, limit=5)
    
    if ships_df.empty:
        pytest.skip("No real ships data available")
    
    # Test with real sequences from database
    # ...
```

### Mixed Strategy Tests

For tests that can use either approach:

```python
@pytest.mark.parametrize("use_real_data", [False, True])
def test_data_strategy_parametrized(use_real_data):
    """Test that works with both mock and real data."""
    manager = TestDataManager(use_real_data=use_real_data)
    ships_df = manager.get_ships_data(with_sequence=True, limit=1)
    
    if use_real_data and ships_df.empty:
        pytest.skip("No real data available")
    
    # Test implementation
```

## File Organization

```
tests/
├── utils/
│   └── test_data_manager.py    # Core test data management
├── test_data/
│   ├── large/                  # Large files (gitignored)
│   │   └── *.fasta
│   └── mock/                   # Small mock files (tracked)
│       └── small_*.fasta
├── test_config.py              # Test configuration
├── test_classification_mixed.py # Example mixed strategy tests
└── README_TEST_DATA_STRATEGY.md # This file
```

## Git Configuration

The `.gitignore` file is configured to:

- **Ignore**: Large sequence files (`tests/test_data/*.fasta`)
- **Track**: Small mock files (`tests/test_data/small_*`, `tests/test_data/mock_*`)

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `STARBASE_TEST_USE_REAL_DATA` | `false` | Enable real data from database |
| `STARBASE_TEST_MAX_REAL_DATA_LIMIT` | `100` | Max records for real data |
| `STARBASE_TEST_DATABASE_URL` | (auto) | Custom database URL |

## Best Practices

### 1. Choose the Right Data Type

- **Mock Data**: Unit tests, function testing, CI/CD
- **Real Data**: Integration tests, BLAST pipeline, classification workflow

### 2. Limit Real Data Usage

```python
# Good: Limit real data to reasonable amounts
ships_df = manager.get_ships_data(with_sequence=True, limit=10)

# Avoid: Fetching all data
ships_df = manager.get_ships_data(with_sequence=True)  # Could be huge
```

### 3. Handle Missing Data Gracefully

```python
def test_with_real_data(real_data_manager):
    ships_df = real_data_manager.get_ships_data(with_sequence=True, limit=5)
    
    if ships_df.empty:
        pytest.skip("No real ships data available in database")
    
    # Test implementation
```

### 4. Use Appropriate Markers

```python
@pytest.mark.real_data
def test_blast_pipeline():
    """Test that requires real data."""
    pass

@pytest.mark.mock_data
def test_unit_function():
    """Test that uses mock data."""
    pass
```

## Migration Guide

### From Old Test Files

1. **Replace hardcoded mock data**:
   ```python
   # Old
   mock_ships = pd.DataFrame({...})
   
   # New
   ships_df = test_data_manager.get_ships_data(with_sequence=True)
   ```

2. **Use fixtures instead of manual setup**:
   ```python
   # Old
   def test_something():
       ships_df = create_mock_ships()
   
   # New
   def test_something(mock_ships_df):
       ships_df = mock_ships_df
   ```

3. **Add real data tests for BLAST/classification**:
   ```python
   @use_real_data
   def test_classification_with_real_data(real_data_manager):
       # Test with real sequences
   ```

### For New Tests

1. **Start with mock data** for unit tests
2. **Add real data tests** for integration/BLAST pipeline
3. **Use appropriate fixtures** and markers
4. **Handle missing data gracefully**

## Troubleshooting

### Database Connection Issues

If real data tests fail due to database issues:

```python
def test_with_real_data(config, database_available):
    if not database_available:
        pytest.skip("Database not available")
    
    # Test implementation
```

### Large File Issues

If you need to add large test files:

1. **Don't commit them to git**
2. **Use the database instead**:
   ```python
   # Instead of large files, use real data
   ships_df = real_data_manager.get_ships_data(with_sequence=True)
   ```
3. **Or use temporary files**:
   ```python
   temp_file = test_data_manager.create_temp_fasta_file(sequences)
   # Use temp_file
   # Cleanup is automatic
   ```

## Examples

See `tests/test_classification_mixed.py` for comprehensive examples of:
- Unit tests with mock data
- Integration tests with real data
- Mixed strategy tests
- Proper error handling
- Cleanup patterns

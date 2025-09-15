"""
Test configuration for managing test data strategies.

This module provides configuration options for controlling how tests use data.
"""

import os
from typing import Dict, Any


class TestConfig:
    """Configuration for test data management."""
    
    # Environment variables that control test behavior
    USE_REAL_DATA_ENV = "STARBASE_TEST_USE_REAL_DATA"
    DATABASE_URL_ENV = "STARBASE_TEST_DATABASE_URL"
    MAX_REAL_DATA_LIMIT_ENV = "STARBASE_TEST_MAX_REAL_DATA_LIMIT"
    
    def __init__(self):
        self.use_real_data = self._get_bool_env(self.USE_REAL_DATA_ENV, False)
        self.database_url = os.getenv(self.DATABASE_URL_ENV)
        self.max_real_data_limit = int(os.getenv(self.MAX_REAL_DATA_LIMIT_ENV, "100"))
        
        # Test data file paths
        self.test_data_dir = os.path.join(os.path.dirname(__file__), "test_data")
        self.large_data_dir = os.path.join(self.test_data_dir, "large")
        self.mock_data_dir = os.path.join(self.test_data_dir, "mock")
    
    def _get_bool_env(self, env_var: str, default: bool) -> bool:
        """Get boolean value from environment variable."""
        value = os.getenv(env_var, str(default)).lower()
        return value in ('true', '1', 'yes', 'on')
    
    def should_use_real_data(self, test_name: str = None) -> bool:
        """
        Determine if real data should be used for a test.
        
        Args:
            test_name: Name of the test (for future filtering logic)
            
        Returns:
            True if real data should be used
        """
        return self.use_real_data
    
    def get_real_data_limit(self, test_type: str = "default") -> int:
        """
        Get the limit for real data based on test type.
        
        Args:
            test_type: Type of test ("blast", "classification", "unit", etc.)
            
        Returns:
            Maximum number of records to fetch
        """
        limits = {
            "blast": min(self.max_real_data_limit, 50),
            "classification": min(self.max_real_data_limit, 100),
            "unit": min(self.max_real_data_limit, 10),
            "integration": min(self.max_real_data_limit, 20),
            "default": min(self.max_real_data_limit, 50)
        }
        return limits.get(test_type, limits["default"])
    
    def get_test_data_path(self, filename: str, use_large: bool = False) -> str:
        """
        Get path to test data file.
        
        Args:
            filename: Name of the test data file
            use_large: If True, look in large data directory
            
        Returns:
            Path to test data file
        """
        if use_large:
            return os.path.join(self.large_data_dir, filename)
        else:
            return os.path.join(self.mock_data_dir, filename)
    
    def is_database_available(self) -> bool:
        """
        Check if database is available for real data tests.
        
        Returns:
            True if database is available
        """
        try:
            from src.database.sql_manager import StarbaseSession
            session = StarbaseSession()
            session.close()
            return True
        except Exception:
            return False


# Global test configuration instance
test_config = TestConfig()


# Pytest markers for different test types
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "real_data: mark test as requiring real data from database"
    )
    config.addinivalue_line(
        "markers", "mock_data: mark test as using mock data only"
    )
    config.addinivalue_line(
        "markers", "blast_pipeline: mark test as part of BLAST pipeline"
    )
    config.addinivalue_line(
        "markers", "classification_pipeline: mark test as part of classification pipeline"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )
    config.addinivalue_line(
        "markers", "unit: mark test as unit test"
    )


# Pytest fixtures for configuration
@pytest.fixture
def config():
    """Fixture providing test configuration."""
    return test_config


@pytest.fixture
def use_real_data_flag(config):
    """Fixture providing real data usage flag."""
    return config.should_use_real_data()


@pytest.fixture
def database_available(config):
    """Fixture indicating if database is available."""
    return config.is_database_available()


# Utility functions for tests
def skip_if_no_real_data(config, test_name: str = None):
    """Skip test if real data is not available."""
    if not config.should_use_real_data(test_name):
        pytest.skip("Real data not enabled")
    
    if not config.is_database_available():
        pytest.skip("Database not available")


def skip_if_no_mock_data():
    """Skip test if mock data is not available (should rarely happen)."""
    # Mock data should always be available, but this is here for completeness
    pass


# Example usage in tests:
"""
@pytest.mark.real_data
def test_with_real_data(config, database_available):
    skip_if_no_real_data(config, "test_with_real_data")
    # Test implementation here

@pytest.mark.mock_data  
def test_with_mock_data():
    # Test implementation here
"""

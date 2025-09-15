"""
Test Data Manager - Mixed Strategy for Test Data

This module implements a mixed strategy for test data:
1. Use real data from database for BLAST/classification pipeline tests
2. Use mock data for other tests
3. Avoid tracking large sequence files in git
"""

import os
import tempfile
import pandas as pd
from typing import Optional, Dict, Any, List
from unittest.mock import patch, MagicMock
import hashlib

from src.database.sql_manager import fetch_meta_data, fetch_ships


class TestDataManager:
    """Manages test data using a mixed strategy of real and mock data."""
    
    def __init__(self, use_real_data: bool = False):
        """
        Initialize the test data manager.
        
        Args:
            use_real_data: If True, use real data from database for BLAST/classification tests
        """
        self.use_real_data = use_real_data
        self._real_data_cache = {}
    
    def get_ships_data(self, curated: bool = False, with_sequence: bool = False, 
                      limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get ships data for testing.
        
        Args:
            curated: If True, only return curated entries
            with_sequence: If True, include sequence data
            limit: Limit number of records returned
            
        Returns:
            DataFrame with ships data
        """
        if self.use_real_data and with_sequence:
            # Use real data from database for BLAST/classification tests
            cache_key = f"ships_curated_{curated}_with_seq_{with_sequence}"
            if cache_key not in self._real_data_cache:
                try:
                    df = fetch_ships(curated=curated, with_sequence=with_sequence, dereplicate=True)
                    if limit:
                        df = df.head(limit)
                    self._real_data_cache[cache_key] = df
                except Exception as e:
                    print(f"Warning: Could not fetch real ships data: {e}")
                    # Fall back to mock data
                    return self._get_mock_ships_data(with_sequence=with_sequence, limit=limit)
            return self._real_data_cache[cache_key]
        else:
            # Use mock data for other tests
            return self._get_mock_ships_data(with_sequence=with_sequence, limit=limit)
    
    def get_meta_data(self, curated: bool = False, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Get metadata for testing.
        
        Args:
            curated: If True, only return curated entries
            limit: Limit number of records returned
            
        Returns:
            DataFrame with metadata
        """
        if self.use_real_data:
            # Use real data from database
            cache_key = f"meta_curated_{curated}"
            if cache_key not in self._real_data_cache:
                try:
                    df = fetch_meta_data(curated=curated)
                    if limit:
                        df = df.head(limit)
                    self._real_data_cache[cache_key] = df
                except Exception as e:
                    print(f"Warning: Could not fetch real metadata: {e}")
                    # Fall back to mock data
                    return self._get_mock_meta_data(limit=limit)
            return self._real_data_cache[cache_key]
        else:
            # Use mock data
            return self._get_mock_meta_data(limit=limit)
    
    def get_test_sequence(self, sequence_type: str = "nucl", length: int = 1000) -> str:
        """
        Get a test sequence for testing.
        
        Args:
            sequence_type: "nucl" for nucleotide, "prot" for protein
            length: Length of sequence to generate
            
        Returns:
            Test sequence string
        """
        if sequence_type == "nucl":
            # Generate a simple nucleotide sequence
            bases = "ATCG"
            import random
            return ''.join(random.choice(bases) for _ in range(length))
        else:
            # Generate a simple protein sequence
            amino_acids = "ACDEFGHIKLMNPQRSTVWY"
            import random
            return ''.join(random.choice(amino_acids) for _ in range(length))
    
    def create_temp_fasta_file(self, sequences: Dict[str, str], 
                              prefix: str = "test_sequence") -> str:
        """
        Create a temporary FASTA file with test sequences.
        
        Args:
            sequences: Dict mapping sequence IDs to sequences
            prefix: Prefix for temporary file
            
        Returns:
            Path to temporary FASTA file
        """
        fd, temp_path = tempfile.mkstemp(suffix='.fasta', prefix=prefix)
        try:
            with os.fdopen(fd, 'w') as f:
                for seq_id, sequence in sequences.items():
                    f.write(f">{seq_id}\n{sequence}\n")
            return temp_path
        except Exception:
            os.close(fd)
            raise
    
    def _get_mock_ships_data(self, with_sequence: bool = False, limit: Optional[int] = None) -> pd.DataFrame:
        """Generate mock ships data for testing."""
        mock_data = {
            "accession_tag": ["SBS000001", "SBS000002", "SBS000003", "SBS000004", "SBS000005"],
            "accession_display": ["SBS000001", "SBS000002", "SBS000003", "SBS000004", "SBS000005"],
            "curated_status": ["curated", "curated", "uncurated", "curated", "uncurated"],
            "familyName": ["Voyager", "Voyager", "Voyager", "Voyager", "Voyager"],
            "name": ["Test Organism 1", "Test Organism 2", "Test Organism 3", "Test Organism 4", "Test Organism 5"],
            "elementBegin": [100, 200, 300, 400, 500],
            "elementEnd": [1100, 1200, 1300, 1400, 1500],
            "contigID": ["contig1", "contig2", "contig3", "contig4", "contig5"]
        }
        
        if with_sequence:
            # Add sequence data
            sequences = [
                "ATGCATGCATGC" * 50,  # 600 bp
                "TTTTTTTTTTTT" * 50,  # 600 bp
                "GGGGGGGGGGGG" * 50,  # 600 bp
                "CCCCCCCCCCCC" * 50,  # 600 bp
                "AAAAAAAAAAAA" * 50,  # 600 bp
            ]
            mock_data["sequence"] = sequences
            mock_data["md5"] = [hashlib.md5(seq.encode()).hexdigest() for seq in sequences]
            mock_data["rev_comp_md5"] = [hashlib.md5(seq[::-1].encode()).hexdigest() for seq in sequences]
        
        df = pd.DataFrame(mock_data)
        if limit:
            df = df.head(limit)
        return df
    
    def _get_mock_meta_data(self, limit: Optional[int] = None) -> pd.DataFrame:
        """Generate mock metadata for testing."""
        mock_data = {
            "accession_tag": ["SBS000001", "SBS000002", "SBS000003", "SBS000004", "SBS000005"],
            "accession_display": ["SBS000001", "SBS000002", "SBS000003", "SBS000004", "SBS000005"],
            "curated_status": ["curated", "curated", "uncurated", "curated", "uncurated"],
            "familyName": ["Voyager", "Voyager", "Voyager", "Voyager", "Voyager"],
            "name": ["Test Organism 1", "Test Organism 2", "Test Organism 3", "Test Organism 4", "Test Organism 5"],
            "taxID": [12345, 12346, 12347, 12348, 12349],
            "strain": ["strain1", "strain2", "strain3", "strain4", "strain5"],
            "order": ["TestOrder", "TestOrder", "TestOrder", "TestOrder", "TestOrder"],
            "family": ["TestFamily", "TestFamily", "TestFamily", "TestFamily", "TestFamily"],
            "elementLength": [1000, 1000, 1000, 1000, 1000],
            "upDR": [100, 100, 100, 100, 100],
            "downDR": [100, 100, 100, 100, 100],
            "contigID": ["contig1", "contig2", "contig3", "contig4", "contig5"],
            "captainID": ["captain1", "captain2", "captain3", "captain4", "captain5"],
            "elementBegin": [100, 200, 300, 400, 500],
            "elementEnd": [1100, 1200, 1300, 1400, 1500],
            "type_element_reference": ["reference1", "reference2", "reference3", "reference4", "reference5"],
            "navis_name": ["navis1", "navis2", "navis3", "navis4", "navis5"],
            "haplotype_name": ["haplo1", "haplo2", "haplo3", "haplo4", "haplo5"],
            "ome": ["ome1", "ome2", "ome3", "ome4", "ome5"],
            "version": ["1.0", "1.0", "1.0", "1.0", "1.0"],
            "genomeSource": ["source1", "source2", "source3", "source4", "source5"],
            "citation": ["citation1", "citation2", "citation3", "citation4", "citation5"],
            "assembly_accession": ["assembly1", "assembly2", "assembly3", "assembly4", "assembly5"]
        }
        
        df = pd.DataFrame(mock_data)
        if limit:
            df = df.head(limit)
        return df


# Pytest fixtures for easy use in tests
@pytest.fixture
def test_data_manager():
    """Fixture providing a test data manager with mock data."""
    return TestDataManager(use_real_data=False)


@pytest.fixture
def real_data_manager():
    """Fixture providing a test data manager with real data from database."""
    return TestDataManager(use_real_data=True)


@pytest.fixture
def mock_ships_df(test_data_manager):
    """Fixture providing mock ships data."""
    return test_data_manager.get_ships_data(with_sequence=True)


@pytest.fixture
def real_ships_df(real_data_manager):
    """Fixture providing real ships data from database."""
    return real_data_manager.get_ships_data(with_sequence=True, limit=10)


@pytest.fixture
def mock_meta_df(test_data_manager):
    """Fixture providing mock metadata."""
    return test_data_manager.get_meta_data()


@pytest.fixture
def real_meta_df(real_data_manager):
    """Fixture providing real metadata from database."""
    return real_data_manager.get_meta_data(limit=10)


@pytest.fixture
def test_sequence(test_data_manager):
    """Fixture providing a test sequence."""
    return test_data_manager.get_test_sequence(length=1000)


@pytest.fixture
def temp_fasta_file(test_data_manager):
    """Fixture providing a temporary FASTA file with test sequences."""
    sequences = {
        "test_seq_1": test_data_manager.get_test_sequence(length=500),
        "test_seq_2": test_data_manager.get_test_sequence(length=750),
        "test_seq_3": test_data_manager.get_test_sequence(length=1000),
    }
    temp_file = test_data_manager.create_temp_fasta_file(sequences)
    yield temp_file
    # Cleanup
    try:
        os.unlink(temp_file)
    except OSError:
        pass


# Context manager for using real data in specific tests
class UseRealData:
    """Context manager for temporarily using real data in tests."""
    
    def __init__(self, manager: TestDataManager):
        self.manager = manager
        self.original_use_real_data = manager.use_real_data
    
    def __enter__(self):
        self.manager.use_real_data = True
        return self.manager
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.manager.use_real_data = self.original_use_real_data


# Decorator for marking tests that should use real data
def use_real_data(func):
    """Decorator to mark tests that should use real data from database."""
    func._use_real_data = True
    return func


# Utility function to check if test should use real data
def should_use_real_data(test_func) -> bool:
    """Check if a test function should use real data."""
    return hasattr(test_func, '_use_real_data') and test_func._use_real_data

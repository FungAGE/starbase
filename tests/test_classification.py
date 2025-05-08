import pytest
import pandas as pd
from src.utils.classification_utils import (
    assign_accession,
    check_contained_match,
    check_similar_match,
    generate_new_accession,
)
from unittest.mock import patch


# Mock data for testing
@pytest.fixture
def mock_ships_df():
    return pd.DataFrame(
        {
            "accession_tag": ["SBS000001", "SBS000002", "SBS000003"],
            "sequence": [
                "ATGCATGCATGC",  # Simple sequence for exact match
                "ATGCATGCATGCATGC",  # Longer sequence for contained match
                "ATGCATGCATTT",  # Similar sequence for similarity match
            ],
        }
    )


@pytest.fixture
def mock_sequence():
    return "ATGCATGCATGC"  # Same as SBS000001


@pytest.fixture
def mock_contained_sequence():
    return "ATGCATGC"  # Contained within SBS000002


@pytest.fixture
def mock_similar_sequence():
    return "ATGCATGCATTA"  # Similar to SBS000003


@pytest.fixture
def mock_similarities():
    """Mock similarity calculation results."""
    return {
        "query_sequence": {
            "SBS000003": 0.85  # Similar to SBS000003
        }
    }


def test_check_contained_match(mock_ships_df, mock_contained_sequence):
    """Test contained sequence matching."""
    # Test contained match - should return SBS000002 as it's longer
    result = check_contained_match(mock_contained_sequence, mock_ships_df)
    assert result == "SBS000002"  # Now matches the longer sequence

    # Test no match
    result = check_contained_match("TTTT", mock_ships_df)
    assert result is None


@patch("src.utils.classification.calculate_similarities")
def test_check_similar_match(
    mock_calc_sim, mock_ships_df, mock_similar_sequence, mock_similarities
):
    """Test similar sequence matching using k-mer comparison."""
    # Set up mock return value
    mock_calc_sim.return_value = mock_similarities

    # Test similar match
    result = check_similar_match(mock_similar_sequence, mock_ships_df, threshold=0.8)
    assert result == "SBS000003"

    # Verify mock was called correctly
    mock_calc_sim.assert_called_once()


def test_generate_new_accession(mock_ships_df):
    """Test new accession number generation."""
    result = generate_new_accession(mock_ships_df)
    assert result == "SBS000004"  # Next number after existing ones

    # Test with empty DataFrame
    empty_df = pd.DataFrame(columns=["accession_tag"])
    result = generate_new_accession(empty_df)
    assert result == "SBS000001"  # First number when no existing accessions


def test_assign_accession_exact_match(mock_ships_df, mock_sequence):
    """Test full workflow with exact match."""
    accession, needs_review = assign_accession(
        mock_sequence, existing_ships=mock_ships_df
    )
    assert accession == "SBS000001"
    assert needs_review is False


def test_assign_accession_contained_match(mock_ships_df, mock_contained_sequence):
    """Test full workflow with contained match."""
    accession, needs_review = assign_accession(
        mock_contained_sequence, existing_ships=mock_ships_df
    )
    assert accession == "SBS000002"
    assert needs_review is True


@patch("src.utils.classification.calculate_similarities")
def test_assign_accession_similar_match(
    mock_calc_sim, mock_ships_df, mock_similar_sequence
):
    """Test full workflow with similar match."""
    mock_calc_sim.return_value = {"query_sequence": {"SBS000003": 0.85}}

    accession, needs_review = assign_accession(
        mock_similar_sequence, existing_ships=mock_ships_df, threshold=0.8
    )
    assert accession == "SBS000003"
    assert needs_review is True


@patch("src.utils.classification.calculate_similarities")
def test_assign_accession_new(mock_calc_sim, mock_ships_df):
    """Test full workflow with new sequence."""
    mock_calc_sim.return_value = {
        "query_sequence": {}  # No similarities
    }

    accession, needs_review = assign_accession("TTTTTTTT", existing_ships=mock_ships_df)
    assert accession == "SBS000004"
    assert needs_review is False


@patch("src.utils.classification.calculate_similarities")
def test_assign_accession_empty_db(mock_calc_sim):
    """Test workflow with empty database."""
    empty_df = pd.DataFrame(columns=["accession_tag", "sequence"])

    # Mock empty similarities result
    mock_calc_sim.return_value = {"query_sequence": {}}

    accession, needs_review = assign_accession("ATGC", existing_ships=empty_df)
    assert accession == "SBS000001"
    assert needs_review is False


def test_assign_accession_invalid_input():
    """Test error handling for invalid input."""
    with pytest.raises(Exception):
        assign_accession(None)

    with pytest.raises(Exception):
        assign_accession("")

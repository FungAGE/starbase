import pytest
from unittest.mock import patch
import pandas as pd

from src.database.sql_manager import fetch_meta_data
from src.components.data import (
    create_ship_accession_modal_data,
    create_accession_modal_data,
)


def test_ship_accession_modal_data_creation(single_accession_meta_data):
    """Test successful creation of ship accession modal data"""
    with patch("src.components.data.fetch_meta_data") as mock_fetch:
        mock_fetch.return_value = pd.DataFrame([single_accession_meta_data])

        with patch("src.components.data.get_quality_tags") as mock_quality:
            mock_quality.return_value = [
                {"tag_type": "missing_direct_repeats", "tag_value": None}
            ]

            result = create_ship_accession_modal_data(
                single_accession_meta_data["ship_accession_tag"]
            )

            assert isinstance(result, dict)
            assert result["title"] == single_accession_meta_data["ship_accession_tag"]
            assert result["familyName"] == "Prometheus"
            assert result["curated_status"] == "curated"


def test_accession_modal_data_creation():
    """Test successful creation of accession modal data"""
    modal_df = pd.DataFrame(
        {
            "accession_tag": ["SSA002851"],
            "familyName": ["Prometheus"],
        }
    )
    with patch("src.components.data.fetch_meta_data") as mock_fetch:
        mock_fetch.return_value = modal_df

        with patch("src.components.data.get_quality_tags") as mock_quality:
            mock_quality.return_value = []

            result = create_accession_modal_data("SSA002851.1")

            assert isinstance(result, dict)
            assert result["title"] == "Group Accession: SSA002851"
            assert result["familyName"] == "Prometheus"


@pytest.mark.integration
def test_single_accession_data(single_accession_meta_data):
    meta_df = fetch_meta_data(
        accessions=single_accession_meta_data["ship_accession_tag"]
    )
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(single_accession_meta_data))


@pytest.mark.integration
def test_multiple_accession_data(multiple_accession_meta_data):
    meta_df = fetch_meta_data(
        accessions=multiple_accession_meta_data["ship_accession_tag"]
    )
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(multiple_accession_meta_data))


@pytest.mark.integration
def test_single_accession_multiple_ships(single_accession_multiple_ships_meta_data):
    meta_df = fetch_meta_data(
        accessions=single_accession_multiple_ships_meta_data["ship_accession_tag"]
    )
    assert not meta_df.empty
    assert meta_df.equals(pd.DataFrame(single_accession_multiple_ships_meta_data))

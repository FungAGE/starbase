import pytest
from unittest.mock import patch
import pandas as pd

from src.components.data import (
    create_ship_accession_modal_data,
    create_accession_modal_data,
)


@pytest.fixture
def ship_accession_modal_data(single_accession_meta_data):
    with (
        patch("src.components.data.fetch_meta_data") as mock_fetch,
        patch("src.components.data.get_quality_tags", return_value=[]),
    ):
        mock_fetch.return_value = pd.DataFrame([single_accession_meta_data])
        return create_ship_accession_modal_data(
            single_accession_meta_data["ship_accession_tag"]
        )


@pytest.fixture
def accession_modal_data(multiple_accession_meta_data):
    with (
        patch("src.components.data.fetch_meta_data") as mock_fetch,
        patch("src.components.data.get_quality_tags", return_value=[]),
    ):
        mock_fetch.return_value = pd.DataFrame(multiple_accession_meta_data)
        return create_accession_modal_data("SSA002851.1")

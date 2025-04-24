from unittest.mock import patch
import pandas as pd
from src.components.callbacks import create_accession_modal


def test_create_accession_modal(mock_accession_data):
    """Test the /api/accession/<accession_id> route."""
    accession_id = "ABC123"

    with patch("src.components.callbacks.fetch_meta_data") as mock_fetch_meta_data:
        mock_fetch_meta_data.return_value = mock_accession_data

        modal_content, modal_title = create_accession_modal(accession_id)

        assert "Ship Accession: ABC123" in str(modal_title)
        assert "Starship Family" in str(modal_content)
        assert "Family1" in str(modal_content)

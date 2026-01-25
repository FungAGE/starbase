import pytest
from unittest.mock import patch


def test_get_ship_accession_details(test_client):
    """Test the /api/ship_accession_details/<ship_accession_id> route."""
    ship_accession_id = "SSA002851"

    with patch("src.api.routes.create_ship_accession_modal_data") as mock_create_modal_data:
        mock_modal_data = {
            "title": ship_accession_id,
            "familyName": "TestFamily",
            "curated_status": "curated"
        }
        mock_create_modal_data.return_value = mock_modal_data

        response = test_client.get(f"/api/ship_accession_details/{ship_accession_id}")

        assert response.status_code == 200
        assert response.json == mock_modal_data
        mock_create_modal_data.assert_called_once_with(ship_accession_id)


def test_get_accession_details(test_client):
    """Test the /api/accession_details/<accession_id> route."""
    accession_id = "SSA002851"

    with patch("src.api.routes.create_accession_modal_data") as mock_create_modal_data:
        mock_modal_data = {
            "title": f"Starship Accession: {accession_id}",
            "familyName": "TestFamily",
            "genomes_present": "1"
        }
        mock_create_modal_data.return_value = mock_modal_data

        response = test_client.get(f"/api/accession_details/{accession_id}")

        assert response.status_code == 200
        assert response.json == mock_modal_data
        mock_create_modal_data.assert_called_once_with(accession_id)


@pytest.mark.parametrize("endpoint", ["/api/blast/blast-submit", "/api/cache/refresh"])
def test_handle_429(test_client, endpoint):
    """Test that rate limits are enforced with a 429 Too Many Requests response."""
    for _ in range(11):  # Exceed rate limit
        test_client.post(endpoint)
    response = test_client.post(endpoint)
    assert response.status_code == 429
    assert response.json == {
        "error": "Too Many Requests",
        "message": "Please wait before making more requests.",
    }

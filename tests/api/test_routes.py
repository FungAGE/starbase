import pytest
from unittest.mock import patch


def test_get_accession_details(test_client):
    """Test the /api/accession/<accession_id> route."""
    accession_id = "ABC123"

    with patch("src.api.routes.create_accession_modal") as mock_create_accession_modal:
        mock_create_accession_modal.return_value = ("bla", "test")

        response = test_client.get(f"/api/accession/{accession_id}")

        assert response.status_code == 200
        assert response.json == {
            "content": "bla",
            "title": "test",
        }

        mock_create_accession_modal.assert_called_once_with(accession_id)


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

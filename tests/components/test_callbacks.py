from src.components.ui import create_modal


def test_create_modal_ship_accession(test_ship_accession_modal_data):
    """Test the create_modal function with ship accession modal data."""
    modal_content, modal_title = create_modal(
        test_ship_accession_modal_data, mode="ship_accession"
    )

    # Check that modal content and title are created
    assert modal_content is not None
    assert modal_title is not None

    # Check title contains accession
    assert test_ship_accession_modal_data["title"] in str(modal_title)

    # Check content contains family name if present
    if test_ship_accession_modal_data.get("familyName"):
        assert test_ship_accession_modal_data["familyName"] in str(modal_content)


def test_create_modal_accession(test_accession_modal_data):
    """Test the create_modal function with accession modal data."""
    modal_content, modal_title = create_modal(
        test_accession_modal_data, mode="accession"
    )

    # Check that modal content and title are created
    assert modal_content is not None
    assert modal_title is not None

    # Check title contains "Starship Accession"
    assert "Starship Accession" in str(modal_title)

    # Check content contains family name if present
    if test_accession_modal_data.get("familyName"):
        assert test_accession_modal_data["familyName"] in str(modal_content)


def test_create_modal_error_handling():
    """Test the create_modal function with error data."""
    error_data = {"title": "Error Test", "error": "Test error message"}

    modal_content, modal_title = create_modal(error_data)

    # Check error handling
    assert modal_content is not None
    assert modal_title is not None
    assert "Error" in str(modal_title)
    assert "Test error message" in str(modal_content)

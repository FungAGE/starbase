from src.components.callbacks import create_modal


def test_create_modal(test_meta_data):
    """Test the create_modal function."""
    accession_id = test_meta_data["accession_tag"][0]
    family_name = test_meta_data["familyName"][0]
    modal_content = create_modal(accession_id)

    assert f"Starship Accession: {accession_id}" in str(modal_content)
    assert "Starship Family" in str(modal_content)
    assert family_name in str(modal_content)
